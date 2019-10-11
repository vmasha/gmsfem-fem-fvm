# gmsfem-fem-fvm
Implementation of the Generalized Multiscale Finite Element Method (GMsFEM) for solution problems in heterogeneous or/and perforated media based on the system representation (fvm or fem fine grid approximation)

My papers about multiscale method:

* Chung ET, Efendiev Y, Li G, Vasilyeva M. Generalized multiscale finite element methods for problems in perforated heterogeneous domains [PDF](https://arxiv.org/abs/1501.03536)
* Vasilyeva M, Chung ET, Efendiev Y, Tyrylgin A. A three-level multi-continua upscaling method for flow problems in fractured porous media [PDF](https://arxiv.org/abs/1810.01581)

*I would be grateful if you add citations to my relevant publications when you use it in your research.*


Implementation of the method contains:

* **fine grid system generation** (./systemFEM.py and ./systemFVM.py) - create mass and stiffness matrices and right - hand side vector (create fine grid system)
* **local domains (coarse grid) generation** (./local-domain/) - create files with coarse cells coordinates and cell indices in local domains
* **multiscale basis function calculation** (./gmsfem-basis/) - solve local spectral problems to generate and save multiscale basis functions (use system generated from first step)
* **projection matrix generation** (./ms-rgen/) - load multiscale basis functions and create projection matrix R
* **fine scale and multiscale solver** (./solver/) - solve fine grid system or/and multiscale solver (use system generated from first step)

Implementation based on the [FEniCS](https://fenicsproject.org) (geometry objects, functions for saving and visualization) and [PETSc](https://www.mcs.anl.gov/petsc/) (matrices, vectors and solvers).

## How to use

Fine grid simulations:

1. run fenics container
  > ```docker run -ti -v $(pwd):/home/fenics/shared quay.io/fenicsproject/stable```
  
2. create folders ./data/out/, ./data/modelF/

3. generate fine grid system 
  > ```python systemFEM.py```
  
  or
  
  > ```python systemFVM.py```
  
4. run fine grid solver in ./solver/
  > ```./solver F 0 20 0.0003 80 6400 ../data/modelF/ ../data/out/ 1 ./ err.txt```

  or 

  > ```./solver F 1 20 0.003 ../data/mesh/mesh-p 6101 ../data/modelF/ ../data/out/ 1 ./ err.txt```

  or

  > ```./solver F 1 20 0.03 ../data/mesh/mesh-h 6373 ../data/modelF/ ../data/out/ 1 ./ err.txt```


Multiscale simulations:

1. create folders ./data/modelMs/omega10fem/, ./data/modelMs/eigen/, ./data/modelMs/dof/

2. local domains generations (coarse grid) in ./local-domain/
  > ```./omegas 2 ../data/omega10/ 10 0.1 0.0 10 0.1 0.0```

  or
  
  > ```./omegas 2 ../data/modelMs/omega10/ 10 0.2 -1.0 10 0.2 -1.0```
  
3. multiscale basis generation in ./gmsfem-basis/
  > ``` ./run```
  
4. generate R in ./ms-rgen/
  > ```./rgen 1 6400 121 ../data/modelMs/ 16```
  
  or
  
  > ```./rgen 1 6101 121 ../data/modelMs/ 16```
  
  or 
  
  > ```./rgen 1 6373 121 ../data/modelMs/ 16```
  
5. solve multiscale in ./solver/
  > ```./solver C 0 20 0.0003 80 6400 ../data/modelF/ ../data/out/ 1936 ../data/modelMs/R100 err.txt```
  
  or
  
  > ```./solver C 1 20 0.003 ../data/mesh/mesh-p 6101 ../data/modelF/ ../data/out/ 1936 ../data/modelMs/R100 err.txt```
  
  or
  
  > ```./solver C 1 20 0.03 ../data/mesh/mesh-h 6373 ../data/modelF/ ../data/out/ 1936 ../data/modelMs/R100 err.txt```
  
