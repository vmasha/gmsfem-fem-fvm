
#include <fstream>
#include "MsRGenerator.h"

using namespace dolfin; 

// ./rgen 1 6400 121 ../data/modelMs/ 16
// ./rgen 1 6101 121 ../data/modelMs/ 16
// ./rgen 1 6373 121 ../data/modelMs/ 16
int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, (char*)0, "");
    parameters["reorder_dofs_serial"] = false;
    parameters["allow_extrapolation"] = true;
    
    // int type = atoi(argv[1]);

    int Nf = atoi(argv[2]);
    int Nomega = atoi(argv[3]);
    std::string msbasisDir = argv[4];
    int ecount = atoi(argv[5]); 

    std::string mort = "m";
    std::string outFileName = msbasisDir + "R100";
    std::string outFileNameDOF = msbasisDir + "dof100";
     
    MsRGenerator generator(msbasisDir);
    generator.generateR(Nomega, ecount, Nf, outFileName, outFileNameDOF, mort);

    PetscFinalize();
    return 0;
}
