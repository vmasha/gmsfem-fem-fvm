#include "load.h"
#include "quadmesh.h"
#include "./space/PressureP0.h"
#include "./space/PressureP1.h"
#include "./space/PressureDG1.h"

using namespace dolfin;

// FVM: 
// ./solver F 0 20 0.0003 80 6400 ../data/modelF/ ../data/out/ 1 ./ err.txt
// ./solver C 0 20 0.0003 80 6400 ../data/modelF/ ../data/out/ 1936 ../data/modelMs/R100 err.txt

// FEM perf: 
// ./solver F 1 20 0.003 ../data/mesh/mesh-p 6101 ../data/modelF/ ../data/out/ 1 ./ err.txt
// ./solver C 1 20 0.003 ../data/mesh/mesh-p 6101 ../data/modelF/ ../data/out/ 1936 ../data/modelMs/R100 err.txt
// FEM heter:
// ./solver F 1 20 0.03 ../data/mesh/mesh-h 6373 ../data/modelF/ ../data/out/ 1 ./ err.txt
// ./solver C 1 20 0.03 ../data/mesh/mesh-h 6373 ../data/modelF/ ../data/out/ 1936 ../data/modelMs/R100 err.txt

// FEM DG1 perf: 
// ./solver F 2 20 1.0e-3 ../data/mesh/mesh-p 34824 ../data/modelF/ ../data/out/ 1 ./ err.txt
// FEM DG1 heter: 
// ./solver F 2 20 1.0e-3 ../data/mesh/mesh-h 37512 ../data/modelF/ ../data/out/ 1 ./ err.txt
int main(int argc, char** argv){
    parameters["reorder_dofs_serial"] = false;
    parameters["allow_extrapolation"] = true;

    Timer timer("solverFC");
    timer.start();

    /*
    * types:
    *   F - fine
    *   C - coarse
    * stypes:
    *   0 - fvm (cell, P0)
    *   1 - fem (vertex, P1)
    *   2 - fem (DG1)
    */
    std::string type = argv[1];
    int stype = std::atoi(argv[2]);

    int timeIt = std::atoi(argv[3]);
    double ta = std::atof(argv[4]);

    int n = std::atoi(argv[6]);
    std::string indir = argv[7];
    std::string outdir = argv[8];

    int Nc = std::atoi(argv[9]); 
    std::string infileR = argv[10];
    std::string errfilename = argv[11];

    // fine mesh and space for saving
    std::shared_ptr<FunctionSpace> W;
    std::shared_ptr<Mesh> mesh;
    if(stype == 1 || stype == 2){ // FEM
        std::string mesh_name = argv[5];
        mesh = std::make_shared<Mesh>(mesh_name + ".xml");
        if(stype == 1)
            W = std::make_shared<PressureP1::FunctionSpace>(mesh);
        else
            W = std::make_shared<PressureDG1::FunctionSpace>(mesh);
    }else{ // FVM
        int NN = std::atoi(argv[5]);
        Mesh mesh2; 
        build(mesh2, NN, NN, 0, 1.0, 0, 1.0);
        mesh = std::make_shared<Mesh>(mesh2);
        W = std::make_shared<PressureP0::FunctionSpace>(mesh);
    }
    Function up(W), up2(W), uerr(W);
    File filep(outdir + "p.pvd"), filep2(outdir + "p2.pvd"), fileerr(outdir + "err.pvd");
    mesh->init(); info(*mesh);


    // system
    std::string infileA = indir + "mat-K.txt";
    std::string infileM = indir + "mat-M.txt";
    std::string infileRhs = indir + "rhs.txt";
    std::string infileDbc = indir + "dbc.txt";
    std::string infileDbcg = indir + "dbc-g.txt";
    
    // time
    double timeMax = ta*timeIt; 
    double dt = timeMax/timeIt;
    int tdel = timeIt/20;
    info("Time = %f and dt = %f", timeMax, dt);
    info("FINE SIZE =  %d", n);
    
    // system parameters
    double uinit = 0.0;

    // ----- load A -----
    MyMatrix mT(infileA, n, n);
    Mat PT = mT.get();
    info("loaded K %d-%d", n, n);
    // ----- load M -----
    MyMatrix mM(infileM, n, n);
    Mat PM = mM.get();
    MatScale(PM, 1.0/dt);
    info("loaded M %d-%d", n, n);
    // ----- load Rhs -----
    MyVector mq(infileRhs, n);
    PETScVector q(mq.get());
    info("loaded Rhs q in (%g, %g)", q.min(), q.max());
    // ----- load DBC -----
    MyVector mb(infileDbc, n);
    std::set<int> bounddofs;
    for(int i = 0; i < n; i++){
        if(mb.getValue(i) > 0.1)
            bounddofs.insert(i);
    }
    MyVector mg(infileDbcg, n);
    PETScVector g(mg.get());
    info("loaded DBC g in (%g, %g)", g.min(), g.max());
    info("Loaded DBC, DOFs on boundary %d from %d vertices", bounddofs.size(), n);

    // init solver
    LinearSolver solverF("default", "default");  
    LinearSolver solverC("default", "default");  
    PETScVector x(PETSC_COMM_SELF, n);
    PETScVector xms(PETSC_COMM_SELF, n);
    PETScVector xerr(PETSC_COMM_SELF, n);
    PETScVector b(PETSC_COMM_SELF, n);
    PETScVector b2(PETSC_COMM_SELF, n);
    Vec pxms, pbc, Pb, Pb2, Aerr;
    PetscScalar vdot, vdot2;
    VecCreateSeq(PETSC_COMM_SELF, n, &pxms);
    VecCreateSeq(PETSC_COMM_SELF, Nc, &pbc);
    VecCreateSeq(PETSC_COMM_SELF, n, &Pb);
    VecCreateSeq(PETSC_COMM_SELF, n, &Pb2);
    VecCreateSeq(PETSC_COMM_SELF, n, &Aerr);
    info("System init %d", n);
    PETScVector bc(PETSC_COMM_SELF, Nc);
    PETScVector xc(PETSC_COMM_SELF, Nc);
    info("System init %d", Nc);

    // A = K + M
    Mat PA;
    MatDuplicate(PT, MAT_COPY_VALUES, &PA);
    MatAXPY(PA,  1., PM, DIFFERENT_NONZERO_PATTERN);

    // Apply DBC
    PetscInt ncols;
    const PetscInt *cols;
    const PetscReal *vals;
    double val, fval;     
    for (std::set<int>::iterator it = bounddofs.begin(); it != bounddofs.end(); ++it) {
        int I = *it;
        MatGetRow(PT, I, &ncols, &cols, &vals);
        for(int ej = 0; ej < ncols; ej++){
            int J = cols[ej];
            double val = (I == J)?1.0:0.0;
            MatSetValues(PA, 1, &I, 1, &J ,&val, INSERT_VALUES);
        }
    }
    MatAssemblyBegin(PA, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(PA, MAT_FINAL_ASSEMBLY);
    info("Apply dbc %d", bounddofs.size());

    PETScMatrix *A = new PETScMatrix(PA);  

    Mat PR, PRT;
    Mat PRA, PRART;
    PETScMatrix *A2;
    if(type == "C"){
        // load R
        MyMatrix mR(infileR, Nc, n);
        PR = mR.get();
        MatTranspose(PR, MAT_INITIAL_MATRIX, &PRT);
        info("COARSE: loaded R %d-%d", Nc, n);
        // mult RART
        MatMatMult(PR, PA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &PRA);
        MatMatMult(PRA, PRT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &PRART);
        A2 = new PETScMatrix(PRART);  
        info("COARSE: matrix generated");
    }

    // initial conditions
    for(int I = 0; I < n; I++){
        x.setitem(I, uinit);
        xms.setitem(I, uinit);
    }
    x.apply("insert"); 
    xms.apply("insert");
    
    // info("TIME Pinit %g", timer1.stop());
    // timer1.start();

    remove(errfilename.data()); 
    std::ofstream errfile(errfilename, std::ios_base::app);

    double t = 0;
    int tcounter = 0;
    while(t <= timeMax && tcounter <= timeIt){
        info("\n Time[%d] %g, dt = %g", tcounter, t, dt);
        
        // ----- PRESSURE FINE -----   
        // rhs
        MatMult(PM, x.vec(),   Pb);// y = Mx
        MatMult(PM, xms.vec(), Pb2);// y = Mx
        for(int I = 0; I < n; I++){
            // f
            VecGetValues(Pb, 1, &I, &val);
            b.setitem(I, val + q.getitem(I));
            // ms
            VecGetValues(Pb2, 1, &I, &val);
            b2.setitem(I, val + q.getitem(I));
        }
        b.apply("insert");
        b2.apply("insert");
        // Apply DBC
        for (std::set<int>::iterator it = bounddofs.begin(); it != bounddofs.end(); ++it) {
            int I = *it;
            b.setitem(I, g.getitem(I));
            b2.setitem(I, g.getitem(I));
        }
        b.apply("insert");
        b2.apply("insert");
        info("system [%d, %d] and rhs in (%g, %g)", n, n, b.min(), b.max());
        // solve
        solverF.solve(*A, x, b);
        info("solve: p in (%g, %g)", x.min(), x.max());

        // ----- PRESSURE COARSE -----   
        if(type == "C"){
            // rhs
            MatMult(PR, b2.vec(), pbc);
            // solve
            for(int I = 0; I < Nc; I++){
                VecGetValues(pbc, 1, &I, &val);
                bc.setitem(I, val);
            }
            bc.apply("insert");
            info("COARSE: system [%d, %d] and rhs in (%g, %g)", Nc, Nc, bc.min(), bc.max());
            solverC.solve(*A2, xc, bc);
            // project 
            MatMult(PRT, xc.vec(), pxms);
            for(int I = 0; I < n; I++){
                VecGetValues(pxms, 1, &I, &val);
                xms.setitem(I, val);
            }
            xms.apply("insert");
            info("COARSE: solve p in (%g, %g)", xms.min(), xms.max());
        }
        
        // error and save
        double x1, x2;
        for(int I = 0; I < n; I++){
            // f
            x1 = x.getitem(I);
            up.vector().get()->setitem(I, x1);
            // ms
            if(type == "C"){
                x2 = xms.getitem(I);
                up2.vector().get()->setitem(I, x2);
                // error
                xerr.setitem(I, x1 - x2);
                uerr.vector().get()->setitem(I, x1-x2);
            }
        }
        up.vector().get()->apply("insert");  
        if(tcounter%tdel == 0)
            filep << up;
        if(type == "C"){
            up2.vector().get()->apply("insert"); 
            uerr.vector().get()->apply("insert"); 
            xerr.apply("insert");
            if(tcounter%tdel == 0){
                filep2 << up2;
                fileerr << uerr;
            }
            // error
            MatMult(PA, xerr.vec(), Aerr);
            VecDot(Aerr, xerr.vec(), &vdot2);
            VecDot(xerr.vec(), xerr.vec(), &vdot);
            double err1 = vdot;
            double err2 = vdot2;
            MatMult(PA, x.vec(), Aerr);
            VecDot(Aerr, x.vec(), &vdot2);
            VecDot(x.vec(), x.vec(), &vdot);
            err1 = std::sqrt(err1/vdot)*100;
            err2 = std::sqrt(err2/vdot2)*100;
            info("ERROR L2 %f, energy %f", err1, err2);
            errfile << tcounter << " " << err1 << " " << err2 << "\n";
        }
        
        t += dt;
        tcounter++; 
    }
    errfile.close();
    info("SUCCESS %d time steps", tcounter);
    // PetscFinalize();
    return 0;
}
