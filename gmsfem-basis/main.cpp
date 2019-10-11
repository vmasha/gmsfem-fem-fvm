#include "load.h"
#include "quadmesh.h"
#include "dofs.h"
#include "./space/PressureP0.h"
#include "./space/PressureP1.h"
#include "./space/PressureDG1.h"

class ShapeFunction : public Expression{
    double vx[8];
    double vy[8];
    double vz[8];
    int mainind;
    double hx, hy, hz;
    int _dim;
 public:
    
    ShapeFunction(int dim, std::string filename):_dim(dim){
        // load cell borders 
        std::ifstream inlines(filename.data());
        std::string line;
        getline(inlines, line);
        std::vector<double> borders;
        readDoubleNumbers(line, borders);
        inlines.close();

        double minX = borders[0];
        double maxX = borders[1];
        hx = maxX - minX;
        double minY = borders[2];
        double maxY = borders[3];
        hy = maxY - minY;
        double minZ, maxZ;
        if(_dim == 3){
            minZ = borders[4];
            maxZ = borders[5];      
            hz = maxZ - minZ;
        }
        // z
        vx[0] = minX;  vy[0] = minY;  vz[0] = minZ;
        vx[1] = maxX;  vy[1] = minY;  vz[1] = minZ;
        vx[2] = minX;  vy[2] = maxY;  vz[2] = minZ;
        vx[3] = maxX;  vy[3] = maxY;  vz[3] = minZ;
        // z
        vx[4] = minX;  vy[4] = minY;  vz[4] = maxZ;
        vx[5] = maxX;  vy[5] = minY;  vz[5] = maxZ;
        vx[6] = minX;  vy[6] = maxY;  vz[6] = maxZ;
        vx[7] = maxX;  vy[7] = maxY;  vz[7] = maxZ;
    }

    void setMain(std::vector<double> vcoords, double meps){
        if(_dim == 2){
            for(int i = 0; i < 4; i++){
                // info("%d: %g, %g", i, vx[i], vy[i]);
                if( std::abs(vx[i] - vcoords[0]) < meps && std::abs(vy[i] - vcoords[1]) < meps )
                    mainind = i;
            }
            // info("%d: %g, %g", mainind, vcoords[0], vcoords[1]);
        }else{
            for(int i = 0; i < 8; i++){
                // info("%d: %g, %g, %g", i, vx[i], vy[i], vz[i]);
                if( std::abs(vx[i] - vcoords[0]) < meps && std::abs(vy[i] - vcoords[1]) < meps && std::abs(vz[i] - vcoords[2]) < meps )
                    mainind = i;
            }
            // info("%d: %g, %g, %g", mainind, vcoords[0], vcoords[1], vcoords[2]);
        }
        info("MAIN VERTEX FOR SHAPE FUNCTIOB %d (%g, %g, %g)", mainind, vx[mainind], vy[mainind], vz[mainind]);       
    }

    void eval(Array<double>& values, const Array<double>& x) const{
        values[0] = 0.0;
        bool inCell;
        if(_dim == 2)
            inCell = (x[0] >= vx[0] && x[0] <= vx[1]) && (x[1] >= vy[0] && x[1] <= vy[2]);
        else
            inCell = (x[0] >= vx[0] && x[0] <= vx[1]) && (x[1] >= vy[0] && x[1] <= vy[2]) && (x[2] >= vz[0] && x[2] <= vz[4]);
        if(!inCell)
            return;
        if(mainind == 0 || mainind == 4){
            values[0] =  (x[0] - vx[1])/hx * (x[1] - vy[3])/hy;
            if(_dim == 3 && mainind == 0)
                values[0] = -(x[0] - vx[1])/hx * (x[1] - vy[3])/hy*(x[2] - vz[4])/hz;
            if(_dim == 3 && mainind == 4)
                values[0] =  (x[0] - vx[1])/hx * (x[1] - vy[3])/hy*(x[2] - vz[0])/hz;
        }
        if(mainind == 1 || mainind == 5){
            values[0] = -(x[0] - vx[0])/hx * (x[1] - vy[2])/hy;
            if(_dim == 3 && mainind == 1)
                values[0] =  (x[0] - vx[0])/hx * (x[1] - vy[2])/hy*(x[2] - vz[5])/hz;
            if(_dim == 3 && mainind == 5)
                values[0] = -(x[0] - vx[0])/hx * (x[1] - vy[2])/hy*(x[2] - vz[1])/hz;
        }
        if(mainind == 2 || mainind == 6){
            values[0] = -(x[0] - vx[3])/hx * (x[1] - vy[1])/hy;
            if(_dim == 3 && mainind == 2)
                values[0] =  (x[0] - vx[3])/hx * (x[1] - vy[1])/hy*(x[2] - vz[6])/hz;
            if(_dim == 3 && mainind == 6)
                values[0] = -(x[0] - vx[3])/hx * (x[1] - vy[1])/hy*(x[2] - vz[2])/hz;
        }
        if(mainind == 3 || mainind == 7){
            values[0] =  (x[0] - vx[2])/hx * (x[1] - vy[0])/hy;
            if(_dim == 3 && mainind == 3)
                values[0] = -(x[0] - vx[2])/hx * (x[1] - vy[0])/hy*(x[2] - vz[7])/hz;
            if(_dim == 3 && mainind == 7)
                values[0] =  (x[0] - vx[2])/hx * (x[1] - vy[0])/hy*(x[2] - vz[3])/hz;
        }
    }
};

class QuadCell : public SubDomain {
 public:
    QuadCell(std::string filename, int dim, double eps): _eps(eps), _dim(dim){
        // load cell borders 
        std::ifstream inlines(filename.data());
        std::string line;
        getline(inlines, line);
        std::vector<double> borders;
        readDoubleNumbers(line, borders);
        inlines.close();

        _minX = borders[0];
        _maxX = borders[1];
        _minY = borders[2];
        _maxY = borders[3];        
        if(_dim == 3){
            _minZ = borders[4];
            _maxZ = borders[5];        
        }
        // info("Quad cell [%g, %g] - [%g, %g]", _minX, _maxX, _minY, _maxY);
    }

    bool inside(const Array<double> &x, bool on_boundary) const{
        if(_dim == 2){
            if (x[0] < (_minX - _eps) || x[0] > (_maxX + _eps) || 
                x[1] < (_minY - _eps) || x[1] > (_maxY + _eps))
                return false;
        }else{
            if (x[0] < (_minX - _eps) || x[0] > (_maxX + _eps) || 
                x[1] < (_minY - _eps) || x[1] > (_maxY + _eps) || 
                x[2] < (_minZ - _eps) || x[2] > (_maxZ + _eps))
                return false;
        }
        return true;
    }

 private:
    double _minX, _maxX, _minY, _maxY, _minZ, _maxZ;
    double _eps;
    int _dim;
};

int main(int argc, char** argv){
    PetscInitialize(&argc, &argv, (char*)0, "");
    parameters["reorder_dofs_serial"] = false;
    parameters["allow_extrapolation"] = true;
    
    /*
    * stypes:
    *   0 - fvm (cell, P0)
    *   1 - fem (vertex, P1)
    */
    int stype = std::atoi(argv[1]);
    std::string outDir = argv[2];
    std::string omDir = argv[3];
    int cind = std::atoi(argv[4]);

    std::string modelDir = argv[6]; 

    int save0 = atoi(argv[7]);
    bool saveit = (save0 == 1)?true:false;

    int ecount =  std::atoi(argv[8]);

    // fine mesh
    std::shared_ptr<Mesh> fineMesh;
    int NL;
    if(stype == 1 || stype == 2){ // FEM
        std::string fine_mesh_name = argv[5];
        fineMesh = std::make_shared<Mesh>(fine_mesh_name + ".xml");
    }else{ // FVM 
        NL = std::atoi(argv[5]); // fine mesh
        Mesh mesh2; 
        build(mesh2, NL, NL, 0.0, 1.0, 0.0, 1.0);
        fineMesh = std::make_shared<Mesh>(mesh2);
    }
    int mdim = fineMesh->geometry().dim();
    double myeps = fineMesh->hmin()*1.0e-5;
    
    std::string resultdir = outDir + "results/";
    std::string dofFile   = outDir + "dof/dof-" + std::to_string(cind);// + ".txt";
    std::string eigenFile = outDir + "eigen/" + std::to_string(cind) + "-m";
    std::string countFile = outDir + "eigen/basisCount-" + std::to_string(cind) + ".txt";
    std::string volFile   = outDir + "eigen/basisVol-" + std::to_string(cind) + ".txt";
    remove(dofFile.data());
    remove(eigenFile.data());
    remove(countFile.data());
    remove(volFile.data());

    std::string filenameM   = modelDir + "mat-S.txt";
    std::string filenameK   = modelDir + "mat-K.txt";
    
    // -----------------
    // ----- OMEGA -----
    // -----------------
    std::string omfile = omDir + "w" + std::to_string(cind) + ".txt";
    std::ifstream inlines(omfile.data());
    std::string line;
    getline(inlines, line);
    std::vector<int> omcells;
    readIntNumbers(line, omcells);
    inlines.close();
    // get coords
    std::string omfile2 = omDir + "v" + std::to_string(cind) + ".txt";
    std::ifstream inlines2(omfile2.data());
    std::string line2;
    getline(inlines2, line2);
    std::vector<double> vcoords;
    readDoubleNumbers(line2, vcoords);
    inlines2.close();
    info("omega with %d cells", omcells.size());
    // fine mesh markers
    MeshFunction<std::size_t> fineSundomForOmega(fineMesh, mdim);
    fineSundomForOmega = 0;
    for (int ii = 0; ii < omcells.size(); ii++){
        int ci = omcells[ii];
        std::string cfile = omDir + "c" + std::to_string(ci) + ".txt";
        QuadCell qcell(cfile, mdim, myeps);
        qcell.mark(fineSundomForOmega, 1);
    }
    

    // -----------------------    
    // ----- LOCAL MESH ------
    // -----------------------
    std::shared_ptr<Mesh> mesh;
    std::shared_ptr<FunctionSpace> W;
    std::shared_ptr<FunctionSpace> gW;
    if(stype == 1 || stype == 2){
        mesh = std::make_shared<SubMesh>(*fineMesh, fineSundomForOmega, 1);
        if(stype == 1){
            W = std::make_shared<PressureP1::FunctionSpace>(mesh);
            gW = std::make_shared<PressureP1::FunctionSpace>(fineMesh);
        }
        if(stype == 2){
            W = std::make_shared<PressureDG1::FunctionSpace>(mesh);
            gW = std::make_shared<PressureDG1::FunctionSpace>(fineMesh);
        }
    }else{
        // all cells
        double minX = vcoords[0], minY = vcoords[1];
        double maxX = vcoords[0], maxY = vcoords[1];
        for (int ii = 0; ii < omcells.size(); ii++){
            int ci = omcells[ii];
            std::string cfile = omDir + "c" + std::to_string(ci) + ".txt";
            // load cell coords
            std::ifstream inlines(cfile.data());
            std::string line;
            getline(inlines, line);
            std::vector<double> borders;
            readDoubleNumbers(line, borders);
            inlines.close();
            // find local domain coord for QUAD cells
            double minXi = borders[0], minYi = borders[2];
            double maxXi = borders[1], maxYi = borders[3];
            minX = (minXi < minX)?minXi:minX;
            minY = (minYi < minY)?minYi:minY;
            maxX = (maxXi > maxX)?maxXi:maxX;
            maxY = (maxYi > maxY)?maxYi:maxY;
        }
        double hh = 1.0/NL;
        int NSx = (maxX - minX + 1.0e-10)/hh;
        int NSy = (maxY - minY+ 1.0e-10)/hh;
        info("LOCAL MESH %d x %d", NSx, NSy);
        Mesh mesh2; 
        build(mesh2, NSx, NSy, minX, maxX, minY, maxY);
        mesh = std::make_shared<Mesh>(mesh2);
        W = std::make_shared<PressureP0::FunctionSpace>(mesh);
        gW = std::make_shared<PressureP0::FunctionSpace>(fineMesh);
    }
    mesh->init(); info(*mesh);
    int Ncells = mesh->num_cells();
    if(saveit){
        File filePvD1(resultdir + "mesh-" + std::to_string(cind) + ".pvd");
        filePvD1 << *mesh;
        info("mesh saved");
    }
    // function
    Function p(W);
    Function gp(gW);
    // size
    int sizeA = p.vector().get()->size();
    int Ngl = gp.vector().get()->size();


    // -----------------------
    // --------- POU ---------
    // -----------------------
    Function pu(W);
    Function pu2(W);
    // all cells
    for (int ii = 0; ii < omcells.size(); ii++){
        int ci = omcells[ii];
        std::string cfile = omDir + "c" + std::to_string(ci) + ".txt";
        ShapeFunction fpou(mdim, cfile);    
        fpou.setMain(vcoords, myeps);

        pu2.vector().get()->zero();
        pu2 = fpou;
        for(int di = 0; di < pu2.vector().get()->size(); di++){
            double val = pu2.vector().get()->getitem(di);
            if(val > 1.0e-10)
                pu.vector().get()->setitem(di, val);
        }
        pu.vector().get()->apply("insert");
    }
    info("pou (%g, %g)", pu.vector().get()->min(), pu.vector().get()->max());
    if(saveit){
        File fileMsbU(resultdir + "pou-" + std::to_string(cind) + ".pvd");
        fileMsbU << pu;
    }


    // ----------------  
    // --- DOF MAP ----
    // ----------------
    std::vector<int> dofmap;
    if(stype == 1 || stype == 2){ // FEM CG1
        // stype == 1
        // std::vector<std::size_t>& localParentMap = mesh->data().array("parent_vertex_indices", mdim);
        // for(int i = 0; i < mesh->num_vertices(); i++){
        //     lgdofmap.push_back(localParentMap[i]);
        // }  
        std::vector<std::size_t>& localParentMap = mesh->data().array("parent_cell_indices", mdim);
        std::vector<int> omcellmap;
        for(int i = 0; i < mesh->num_cells(); i++){
            omcellmap.push_back(localParentMap[i]);
        } 
        std::map<int, Point> coords = generateCoords(*mesh, *W);
        std::map<int, Point> fcoords = generateCoords(*fineMesh, *gW);
        info("fff DOFs size global = %d, local = %d", Ngl, sizeA);
        std::map<int, int> fmap = findLGMap(    *mesh,  *W,  coords,  0, sizeA, 
                                            *fineMesh, *gW, fcoords,  0, Ngl, omcellmap);
        for(std::map<int, int>::iterator it = fmap.begin(); it != fmap.end(); it++)
            dofmap.push_back(it->second);
    }else{ // FVM 
        dofmap.resize(Ncells);
        for(int i = 0; i < mesh->num_cells(); i++){
            Cell lcell(*mesh, i);
            Point lpoint = lcell.midpoint();
            for(int j = 0; j < fineMesh->num_cells(); j++){
                if(fineSundomForOmega[j] == 1){
                    Cell gcell(*fineMesh, j);
                    if(lpoint.distance(gcell.midpoint()) < myeps)
                        dofmap[i] = j;
                }
            }
        }
    } 
    info("DOFs size global = %d, local = %d", Ngl, sizeA);
    
    // ------------------------------
    // ------ LOAD GLOBAL A, M ------
    // ------------------------------
    MyMatrix matK(filenameK, Ngl, Ngl, 20);
    Mat GA = matK.get();
    MyMatrix matM(filenameM, Ngl, Ngl, 20);
    Mat GS = matM.get();
    info("LOAD GLOBAL DOF and Matrices %d", Ngl);

    // ------------------------
    // ------ LOCAL A, M ------
    // ------------------------
    Mat Ap, Mp;
    MatCreateSeqAIJ(PETSC_COMM_SELF, sizeA, sizeA, sizeA, 0, &Ap);
    MatCreateSeqAIJ(PETSC_COMM_SELF, sizeA, sizeA, sizeA, 0, &Mp);
    PetscInt ncols, ncols2;
    const PetscInt *cols, *cols2;
    const PetscReal *vals, *vals2;
    double val, sumval, dval; 
    std::vector<int>::iterator dmit;    
    for(int li = 0; li < sizeA; li++){
        int gi = dofmap[li];
        // A
        MatGetRow(GA, gi, &ncols, &cols, &vals);
        sumval = 0; dval = 0;
        for(int j = 0; j < ncols; j++){
            int gj = cols[j];
            dmit = std::find(dofmap.begin(), dofmap.end(), gj);
            if(gi != gj && dmit != dofmap.end()){
                int lj = std::distance(dofmap.begin(), dmit);
                val = vals[j];
                MatSetValues(Ap, 1, &li, 1, &lj ,&val, INSERT_VALUES);
                sumval += -val;
            }
            if(gi == gj)
                dval = vals[j];
            
        } 
        MatSetValues(Ap, 1, &li, 1, &li ,&sumval, INSERT_VALUES);
        // S
        MatGetRow(GS, gi, &ncols2, &cols2, &vals2);
        sumval = 0;
        for(int j = 0; j < ncols2; j++){
            int gj = cols2[j];
            dmit = std::find(dofmap.begin(), dofmap.end(), gj);
            if(gi != gj && dmit != dofmap.end()){
                int lj = std::distance(dofmap.begin(), dmit);
                val = vals2[j];
                MatSetValues(Mp, 1, &li, 1, &lj ,&val, INSERT_VALUES);
                sumval += val;
            }
            if(gi == gj)
                dval = vals2[j];
        } 
        if(stype == 1)
            MatSetValues(Mp, 1, &li, 1, &li ,&sumval, INSERT_VALUES);
        else
            MatSetValues(Mp, 1, &li, 1, &li ,&dval, INSERT_VALUES);
    }
    MatAssemblyBegin(Ap, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Ap, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Mp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Mp, MAT_FINAL_ASSEMBLY);
    info("create A, M[%d, %d]", sizeA, sizeA);

    // ----------------------------
    // --------- SPECTRAL BASIS ---
    // ----------------------------
    auto As = std::make_shared<PETScMatrix>(Ap);
    auto Ms = std::make_shared<PETScMatrix>(Mp);
    Mat PR, PRT;
    PETScVector uv(PETSC_COMM_SELF, sizeA);

    // Eigen SOLVER
    SLEPcEigenSolver *esolver = new SLEPcEigenSolver(MPI_COMM_SELF, As, Ms);
    esolver->parameters["solver"] = "krylov-schur";
    // esolver->parameters["problem_type"] = "gen_hermitian";
    esolver->parameters["spectral_transform"] = "shift-and-invert";
    esolver->parameters["spectral_shift"] = 1.0e-5;
    // esolver->parameters["spectral_shift"] = 1.0e-3;//1.0e-5;
    // esolver->parameters["spectrum"] = "smallest magnitude";
     esolver->parameters["spectrum"] = "smallest magnitude";// magnitude real
    esolver->solve(ecount);
    int ccount = esolver->get_number_converged();
    ecount = (ecount > ccount)?ccount:ecount;
    info("ecount %d ", ecount);

    double lmd, v;
    PETScVector rx, vx;   
    // int eigcounter = 0;
    // for(int i = 0; i < ecount; i++){
    //     esolver->get_eigenpair(lmd, v, rx, vx, i);
    //     info("lambda[%d] = %g, z in (%g, %g)", i, lmd, rx.min(), rx.max());
    //     // if(lmd > -1.0e-3)
    //     // if(lmd < 1.0e-1 && lmd > -1.0e-3)
    //         eigcounter++;
    // }
    // info("ECOUNT %d", ecount);
    // ecount = eigcounter;
    // info("ECOUNT %d", ecount);
    
    // save dof
    PetscInt arraydof[sizeA];
    for(int i = 0; i < sizeA; i++)
        arraydof[i] = dofmap[i];
    int fdd;
    PetscViewer viewerDof;
    PetscViewerBinaryOpen(PETSC_COMM_SELF, dofFile.data(), FILE_MODE_WRITE, &viewerDof);
    PetscViewerBinaryGetDescriptor(viewerDof, &fdd);
    PetscBinaryWrite(fdd, &arraydof, sizeA, PETSC_INT, PETSC_FALSE);
    PetscViewerDestroy(&viewerDof);
    // save bin
    int fd1;
    PetscViewer viewerSol1;
    PetscViewerBinaryOpen(PETSC_COMM_SELF,  eigenFile.data(), FILE_MODE_WRITE, &viewerSol1);
    PetscViewerBinaryGetDescriptor(viewerSol1, &fd1);
    // vol
    std::ofstream ofileVol(volFile, std::ios_base::app);   

    PetscScalar *array;
    File filep0(resultdir + "ev-m.pvd");
    int ccounter = 0;
    for(int i = 0; i < ecount; i++){
        double vol = 1.0;
        esolver->get_eigenpair(lmd, v, rx, vx, i);
        // if(lmd < 1.0e-1 && lmd > -1.0e-3){
        if(lmd > -1.0e-3){
            ofileVol << ccounter << " " << vol << "\n";
            info("lambda[%d] = %g, z in (%g, %g)", i, lmd, rx.min(), rx.max());

            p.vector().get()->zero();
            for(int jj = 0; jj < sizeA; jj++){
                double mval = rx.getitem(jj);
                p.vector().get()->setitem(jj, mval);
                double vpou = pu.vector().get()->getitem(jj);
                double rval =  vpou * mval;
                rx.setitem(jj, rval);
                // ebasisFile << rval << std::endl;
            }
            rx.apply("insert");
            p.vector().get()->apply("insert");

            VecGetArray(rx.vec(), &array);
            PetscBinaryWrite(fd1, array, sizeA, PETSC_DOUBLE, PETSC_FALSE);

            if(saveit)
                filep0 << p;
            ccounter++;
        }
    }
    // ebasisFile.close();
    PetscViewerDestroy(&viewerSol1);    
    ofileVol.close();

    // count
    std::ofstream ofile(countFile, std::ios_base::app);   
    ofile << ccounter << std::endl << sizeA << "\n";
    ofile.close();
    
 
    info("SUCCESS %d", ccounter);
//    PetscFinalize();
    return 0;
}
