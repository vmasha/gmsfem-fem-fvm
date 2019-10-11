#ifndef DOFS_H
#define	DOFS_H

#include <dolfin.h>
#include <cmath>

using namespace dolfin;

std::map<int, Point> generateCoords(Mesh& mesh, FunctionSpace& V) {
    // Timer timer("aa"); 
    // timer.start();
    int dim = mesh.geometry().dim();
    const GenericDofMap* dofMap = V.dofmap().get();
    std::vector<dolfin::la_index> dofs = dofMap->dofs();
    std::vector<double> coords = V.tabulate_dof_coordinates();
    std::map<int, Point> coordMap;
    for(int i = 0; i < dofs.size(); i++){
        int dof = dofs[i];
        if(dim == 3){
            Point p(coords[dim*i+0], coords[dim*i+1], coords[dim*i+2]);
            coordMap.insert(std::pair<int, Point>(dof, p));
        }else{
            Point p(coords[dim*i+0], coords[dim*i+1]);
            coordMap.insert(std::pair<int, Point>(dof, p));
        }
    }
    // info("generateCoords time %f, coordMap %d", timer.stop(), coordMap.size());
    info("generateCoords %d for dofs %d", coordMap.size(), dofs.size());
    return coordMap;
}

std::map<int, int> findLGMap(Mesh& lMesh, FunctionSpace& lV, std::map<int, Point> &lCoordMap, int cN0, int cN1,
Mesh& fMesh, FunctionSpace& fV, std::map<int, Point> &fCoordMap, int fN0, int fN1, 
std::vector<int>& localParentMap){
    // Timer timer("aa"); 
    // timer.start();
    double heps = lMesh.hmin()/10;
    int dim = lMesh.geometry().dim();
    const GenericDofMap* lDofMap = lV.dofmap().get();
    const GenericDofMap* fDofMap = fV.dofmap().get();
    
    std::map<int, int> localDofToFineDof;
    for(int lci = 0; lci < lMesh.num_cells(); lci++){
        // local mesh dofs
        ArrayView<const dolfin::la_index> lInds = lDofMap->cell_dofs(lci);
        // parent mesh dofs
        int pci = localParentMap[lci];
        ArrayView<const dolfin::la_index> pInds = fDofMap->cell_dofs(pci);
        for(int lp = 0; lp < lInds.size(); lp++){
            int lDof = lInds[lp];
            if(lDof < cN1 && lDof >= cN0){
                Point lPoint = lCoordMap[lDof];
                for(int pp = 0; pp < pInds.size(); pp++){
                    int pDof = pInds[pp];
                    if(pDof < fN1 && pDof >= fN0){
                        Point pPoint = fCoordMap[pDof];
                        double d = lPoint.distance(pPoint);
                        if(d < heps)
                            localDofToFineDof.insert(std::pair<int,int>(lDof, pDof));
                    }// check for current part of space dof
                }
            }
        }
    }
    // info("findLocalDofsToFineDofsUniv time %f", timer.stop() );
    return localDofToFineDof;
}

#endif
