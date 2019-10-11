#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>

int getIndex(int i, int j, int Nx){
    return (j * Nx + i); 
    // return (i * Nx + j); 
}

class MyCell{
 public:
    int ci;
    std::vector<int> vertices;
    std::vector<double> borders;
    MyCell(int i, int j, int Nx, double hx, double hy, double x0, double y0){
        ci = getIndex(i, j, Nx);
        vertices.push_back(getIndex(i,     j, Nx+1));
        vertices.push_back(getIndex(i+1,   j, Nx+1));
        vertices.push_back(getIndex(i,   j+1, Nx+1));
        vertices.push_back(getIndex(i+1, j+1, Nx+1));
        
        borders.push_back(x0 + i*hx);
        borders.push_back(x0 + (i+1)*hx);
        borders.push_back(y0 + j*hy);
        borders.push_back(y0 + (j+1)*hy);
    }

    void info(){
        printf("c%d: ", ci);
        for(int jj = 0; jj < borders.size(); jj++){
            printf("%g, ", borders[jj]);
        }
        printf("\n");
    }

    bool containVertex(int item){
        bool cont = false;
        for(int i = 0; i < vertices.size(); i++){
            if(vertices[i] == item )
                cont = true;
        }
        return cont;//(std::find(vertices.begin(), vertices.end(), item) != vertices.end());
    }
};

// FEM: ./omegas 2 ../data/modelMs/omega10-fem/ 10 0.2 -1.0 10 0.2 -1.0
// FVM: ./omegas 2 ../data/modelMs/omega10-fvm/ 10 0.1 0.0 10 0.1 0.0
int main(int argc, char** args) {
    /*
    * type: 
    *       2 - CG 2d
    */
    int type =  std::atoi(args[1]);
    std::string outfile = args[2];

    int Nx = std::atoi(args[3]);
    double hx = std::atof(args[4]);
    double x0 = std::atof(args[5]);
    
    int Ny = std::atoi(args[6]);
    double hy = std::atof(args[7]);
    double y0 = std::atof(args[8]);
            
    if(type == 2){// CG 2d
        // generate cells
        std::vector<MyCell> cells;
        for(int i = 0; i < Nx; i++)
            for(int j = 0; j < Ny; j++)
                cells.push_back( MyCell(i, j, Ny, hx, hy, x0, y0) );
        for(int i = 0; i < cells.size(); i++){
            cells[i].info();
            std::ofstream outdomain(outfile + "c" + std::to_string(cells[i].ci) + ".txt");
            std::vector<double> borders = cells[i].borders;
            for(int jj = 0; jj < borders.size(); jj++){
                outdomain << borders[jj] << " ";
            }
        }
        // printf("cells %d\n", cells.size()); 
        // find vertices cells and save to file
        for(int i = 0; i < Nx+1; i++){// vertex index
            for(int j = 0; j < Ny+1; j++){// vertex index
                int vi = getIndex(i, j, Nx+1);
                double xx = x0 + i*hx;
                double yy = y0 + j*hy;
                printf("\n omega %d:", vi);
                // vertex coord
                std::ofstream outdomain2(outfile + "v" + std::to_string(vi) + ".txt");
                outdomain2 << xx << " " << yy;
                outdomain2.close();
                // cells
                std::ofstream outdomain(outfile + "w" + std::to_string(vi) + ".txt");
                for(int ci = 0; ci < cells.size(); ci++){
                    if(cells[ci].containVertex(vi)){
                        outdomain << cells[ci].ci << " ";
                        printf(" %d", cells[ci].ci );
                    }
                }
                outdomain.close();
            }
        }
    }

    printf("\n generated \n");
    return 0;
}
