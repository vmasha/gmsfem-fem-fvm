#ifndef GMSFEMPROBLEM_H
#define GMSFEMPROBLEM_H

#include <math.h>
#include <dolfin.h>

using namespace dolfin;

class LocalSolution{
 public:
    int ecount, esize;
    std::vector<int> dofs;
    std::vector<double> vols;
    std::vector<int> cmarkers;
    std::vector<std::vector<double>* > bfunctions;

    std::pair<int, int> init(std::string fileNameCount){
        std::ifstream efile(fileNameCount.data());
        efile >> ecount >> esize;
        efile.close();
        std::pair<int, int> pp(ecount, esize);
        return pp;
    }

    int loadDOF(std::string fileNameDOF){
        dofs.clear();
        int fd;
        PetscViewer viewerDof;
        PetscViewerBinaryOpen(PETSC_COMM_SELF, fileNameDOF.data(), FILE_MODE_READ, &viewerDof);
        PetscViewerBinaryGetDescriptor(viewerDof, &fd);

        PetscInt avec[esize];
        PetscBinaryRead(fd, avec, esize, PETSC_INT);
        for(int li = 0; li < esize; li++)
            dofs.push_back(avec[li]);

        PetscViewerDestroy(&viewerDof);
        return dofs.size();
    }

    void loadVol(std::string vfilename, int ecount){
        std::ifstream efile(vfilename.data());
        int ccc = 1;
        double val = 0;
        int mc = 0;
        for(int ei = 0; ei < ecount; ei++){
            efile >> mc >> val;
            mc++;
            cmarkers.push_back(mc);
            vols.push_back(val);
        }
        efile.close();
    }

    void loadBasis(std::string bfilename, int bcc){
        int fd;
        PetscViewer viewerDof;
        PetscViewerBinaryOpen(PETSC_COMM_SELF, bfilename.data(), FILE_MODE_READ, &viewerDof);
        PetscViewerBinaryGetDescriptor(viewerDof, &fd);

        PetscScalar avec[esize];
        for(int ei = 0; ei < bcc; ei++){
            PetscBinaryRead(fd, avec, esize, PETSC_SCALAR);
            std::vector<double> *ebasis = new std::vector<double>();
            for(int li = 0; li < esize; li++)
                ebasis->push_back(avec[li]);
            bfunctions.push_back(ebasis);
        }

        PetscViewerDestroy(&viewerDof);
    }

    std::vector<double>& getBasis(int ei){
        return *bfunctions[ei];
    }

};

class MsRGenerator {
public:
    MsRGenerator(std::string indir) {
        timer = new Timer("solve");
        _indir = indir;
    };

    // classic and save in txt
    void generateR(int Nom, int ecount, int Nf, std::string outfileR, std::string outfileDOF, std::string mort){
        timer->start();
        // LOAD
        int maxNZ = 0;
        int countfracdomains = 0;
        std::vector<LocalSolution*> localsols;
        std::vector<int> localsols_ecount;
        int ecounter = 0;
        for(int vi = 0; vi < Nom; vi++){
            std::string filenameDOF = _indir + "dof/dof-" + std::to_string(vi);// + ".txt";
            std::string bfilename = _indir + "eigen/" + std::to_string(vi) + "-" + mort;
            std::string vfilename = _indir + "eigen/basisVol-" + std::to_string(vi) + ".txt";
            std::string cfilename = _indir + "eigen/basisCount-" + std::to_string(vi) + ".txt";
            // pore matrix
            LocalSolution *ld = new LocalSolution();
            std::pair<int, int> bc = ld->init(cfilename);

            int bcf = (bc.first < ecount)?bc.first:ecount;
            int sizem = bc.second;

            ld->loadDOF(filenameDOF);
            ld->loadBasis(bfilename, bcf);
            ld->loadVol(vfilename, sizem);
            
            // if(vi%(Nom/10) == 0)
                info("m %d: dof %d  and ecount %d", vi, sizem, bcf);
            if(sizem != 0){
                localsols.push_back(ld);
                localsols_ecount.push_back(bcf);
                ecounter += bcf;
            }else{
                info("om %d NO BASIS", vi);
            }
            maxNZ = (maxNZ < sizem)?sizem:maxNZ;
        }
        info("Nom %d, Nc %d", localsols.size(), ecounter);

        // put to R
        int colCount = Nf;
        int rowCount = ecounter;
        MatCreateSeqAIJ(PETSC_COMM_SELF, rowCount, colCount, maxNZ, 0, &PR);
        info("R[%d, %d] with maxnnz = %d", rowCount, colCount, maxNZ);        

        // for saving
        std::string bufferDof = "";
        std::string bufferMat = "";
        
        int counter = 0;
        double p0 = 0.0;
        for(int omi = 0; omi < localsols.size(); omi++){
            // info("m %d", omi);
            int ecountm = localsols_ecount[omi];
            LocalSolution *localsolution = localsols[omi];
            std::vector<int> lddofs = localsolution->dofs;
            for(int ei = 0; ei < ecountm; ei++){
                // save dof 
                int mcnt = localsolution->cmarkers[ei];
                double vol = localsolution->vols[ei];//(ei == 0)?0.0025:0.05;//0.00125:0.0353;//0.0025:0.05;
                bufferDof += std::to_string(counter) + " " + std::to_string(omi) + " " + std::to_string(mcnt) + " " 
                          +  std::to_string(p0) + " " + std::to_string(vol) + "\n";
                // info("i %d, mc %d, vol %g", omi, mcnt, vol);
                // set R
                std::vector<double> &mbasis = localsolution->getBasis(ei);
                for(int li = 0; li < lddofs.size(); li++){
                    int gi = lddofs[li];
                    double val = mbasis[li];
                    // info("f: %d %g", li, val);
                    MatSetValues(PR, 1, &counter, 1, &gi ,&val, INSERT_VALUES);
                    bufferMat += std::to_string(counter) + " " + std::to_string(gi) + " " + std::to_string(val) + "\n";
                }
                counter++;
            }
        }
        MatAssemblyBegin(PR, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(PR, MAT_FINAL_ASSEMBLY);

        // dof file
        remove(outfileDOF.c_str());
        std::ofstream doffile(outfileDOF.c_str(), std::ios_base::app);
        doffile << bufferDof;
        doffile.close();
        
        // mat file
        remove(outfileR.c_str());
        std::ofstream matfile(outfileR.c_str(), std::ios_base::app);
        matfile << bufferMat;
        matfile.close();

        info("generate global R[%d, %d], time=%g", rowCount, colCount, timer->stop());
    }

    
protected:
    Timer *timer;
    std::string _indir;
    std::string _indir2;
    std::vector<int> celloffcets_start;
    std::vector<int> celloffcets_end;
    std::map<int, std::vector<int>> mapvertexcells;

    Mat PR;

    void saveMatrixBin(PETScMatrix& S, std::string fname){
        Mat PA1 = S.mat();
        saveMatrixBin(PA1, fname);
    }
    
    void saveMatrixBin(Mat& PA1, std::string fname){
        int fd, fd1, fd2;
        std::string mname1 = fname + "1.dat";
        std::string mname2 = fname + "2.dat";
        std::string mnameo = fname + "o.dat";
        info("START: save %s", fname.data());

        PetscInt dim1, dim2;
        MatGetSize(PA1, &dim1, &dim2);
        
        PetscInt ncols;
        const PetscInt *cols;
        const PetscReal *vals;
        PetscInt olen[dim1];
        int maxnz = 0;
        for(int row = 0; row < dim1; row++){ // rows
            MatGetRow(PA1, row, &ncols, &cols, &vals);
            maxnz = (maxnz < ncols)?ncols:maxnz;
            olen[row] = ncols;
        }
        PetscViewer viewer1o;
        PetscViewerBinaryOpen(PETSC_COMM_SELF, mnameo.data(), FILE_MODE_WRITE, &viewer1o);
        PetscViewerBinaryGetDescriptor(viewer1o, &fd);
        PetscBinaryWrite(fd, &olen, dim1, PETSC_INT, PETSC_FALSE);
        PetscViewerDestroy(&viewer1o);
                 
        PetscViewer viewer11, viewer12;
        PetscViewerBinaryOpen(PETSC_COMM_SELF, mname1.data(), FILE_MODE_WRITE, &viewer11);
        PetscViewerBinaryOpen(PETSC_COMM_SELF, mname2.data(), FILE_MODE_WRITE, &viewer12);
        PetscViewerBinaryGetDescriptor(viewer11, &fd1);
        PetscViewerBinaryGetDescriptor(viewer12, &fd2);
        PetscReal vals2[maxnz];
        PetscInt cols2[maxnz];
        for(int row = 0; row < dim1; row++){ // rows
            MatGetRow(PA1, row, &ncols, &cols, &vals);
            for(int j = 0; j < ncols; j++){
                vals2[j] = vals[j];
                cols2[j] = cols[j];
            }
            PetscBinaryWrite(fd1, vals2, ncols, PETSC_SCALAR, PETSC_FALSE);//PETSC_SCALAR, PETSC_FALSE);
            PetscBinaryWrite(fd2, cols2, ncols, PETSC_INT, PETSC_FALSE);//PETSC_SCALAR, PETSC_FALSE);
            // if(row == 25){
            //     info("--- ncols %d ---", ncols);
            //     for(int j = 0; j < ncols; j++){
            //         //    info("c %d, v %g", cols[j], vals[j]);
            //         if(std::abs(vals[j]) > 1.0e-10)
            //             info("%d with nz %d: c %d, v %g", row, ncols, cols[j], vals[j]);
            //     }
            // }
        }
        //MatView(PA1, viewer1);
        info("END: save %s [%d, %d] with nz %d", fname.data(), dim1, dim2, maxnz);
        PetscViewerDestroy(&viewer11);
        PetscViewerDestroy(&viewer12);
    }
    
};

#endif