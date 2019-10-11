#ifndef LOAD_H
#define LOAD_H

#include <dolfin.h>
#include <fstream>

using namespace dolfin;

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

double readDoubleNumbers(const std::string &s, std::vector<double> &v){
    std::istringstream is(s);
    double n;
    while( is >> n ) 
        v.push_back( n );
    return v.size();
}

int readIntNumbers(const std::string &s, std::vector<int> &v){
    std::istringstream is(s);
    int n;
    while( is >> n ) 
        v.push_back( n );
    return v.size();
} 

void loadVector(std::string fileVec, std::vector<double> &vals){
    std::ifstream infile(fileVec.data());
    std::string line;
    while( getline(infile,line) ){
        std::vector<std::string> vec = split(line, ' ');
        int dof = std::atoi(vec[0].data());
        double val = std::atof(vec[1].data());
        vals.push_back(val);
    }
    infile.close();
}

class MyVector {
  public:
    MyVector(std::string infile, int n){
        info(" ====== LOAD VECTOR %d from %s ====== ", n, infile.data());
        std::ifstream inRhs(infile.data());
        VecCreateSeq(PETSC_COMM_SELF, n, &vp);
        std::string line;
        while( getline(inRhs, line) ){
            std::vector<std::string> vec = split(line, ' ');
            int i = std::atoi(vec[0].data());
            double val = std::atof(vec[1].data());
            VecSetValues(vp, 1, &i, &val, INSERT_VALUES);
        }
        VecAssemblyBegin(vp);
        VecAssemblyEnd(vp);
        inRhs.close();
        info(" ====== LOADED ====== ");
    }

    double getValue(int i){
        double val;
        VecGetValues(vp, 1, &i, &val );
        return val;
    }

    Vec& get(){
        return vp;
    }

 private:
    Vec vp;
};

class MyMatrix {
  public:
    MyMatrix(std::string infile, int n1, int n2, int nnz){
        info(" ====== LOAD MATRIX [%d, %d] from %s ====== ", n1, n2, infile.data());
        MatCreateSeqAIJ(PETSC_COMM_SELF, n1, n2, nnz, 0, &Ap);
        std::string line;
        std::ifstream inMat(infile.data());
        while( getline(inMat, line) ){
            std::vector<std::string> vec = split(line, ' ');
            int i = std::atoi(vec[0].data());
            int j = std::atoi(vec[1].data());
            double val = std::atof(vec[2].data());
            // set value
            MatSetValues(Ap, 1, &i, 1, &j ,&val, INSERT_VALUES);
        }
        inMat.close();
        MatAssemblyBegin(Ap, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Ap, MAT_FINAL_ASSEMBLY);
        info(" ====== LOADED ====== ");
    }

    Mat& get(){
        return Ap;
    }

 private:
    Mat Ap;
};

#endif