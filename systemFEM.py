import numpy as np
import sys, math, time
from dolfin import *
import petsc4py
from petsc4py import PETSc
petsc4py.init(sys.argv)
print("FEM system generator")

parameters["reorder_dofs_serial"] = False;
parameters["allow_extrapolation"] = True;

outdir = './data/modelF/'

# perforated
mesh_name  = './data/mesh/mesh-p'
mheat = 5

# heterogeneous
# mesh_name  = './data/mesh/mesh-h'
# mheat = 2

c1 = Constant(1.0)
c2 = Constant(0.01)
k1 = Constant(1.0)
k2 = Constant(100.0)
f = Constant(0.0)

# mesh
mesh = Mesh(mesh_name + ".xml")
subdomains = MeshFunction("size_t", mesh, mesh_name + "_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, mesh_name + "_facet_region.xml")

# define system 
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

ds = ds(subdomain_data=boundaries)
dx = dx(subdomain_data=subdomains)

# forms
a = inner(k1*grad(u), grad(v))*dx(1) + inner(k2*grad(u), grad(v))*dx(2)
m = c1*u*v*dx(1) + c2*u*v*dx(2)
s = k1*u*v*dx(1) + k2*u*v*dx(2)
L = f*v*dx

M = PETScMatrix()
assemble(m, tensor=M)
A = PETScMatrix()
assemble(a, tensor=A)
S = PETScMatrix()
assemble(s, tensor=S)
b = PETScVector()
assemble(L, tensor=b)

n = b.size()
print('system init, SIZE = ' + str(n))

# DBC
mg = PETSc.Vec().createSeq(n) 
g = PETSc.Vec().createSeq(n) 
for fi in range(boundaries.size()):
    mbc = boundaries[fi]
    if(mbc == mheat):#mbc == 1 or mbc == 2 or mbc == 3 or mbc == 4 or mbc == 5):
        facet = Facet(mesh, fi)
        mg.setValue(facet.entities(0)[0], 1.0)
        mg.setValue(facet.entities(0)[1], 1.0)
        g.setValue(facet.entities(0)[0], 0.0)
        g.setValue(facet.entities(0)[1], 0.0)
        if(mbc == mheat):
            g.setValue(facet.entities(0)[0], 1.0)
            g.setValue(facet.entities(0)[1], 1.0)
print('generate Dbc') 

# save DBC
outRhs = outdir + 'dbc.txt'
fileRhs = open(outRhs, "w")
bufferRhs = ''
for I in range(n):
    bufferRhs += str(I) + ' ' + str(mg.getValue(I)) + '\n'
fileRhs.write(bufferRhs)
fileRhs.close()
print('save rhs into ' + outRhs)      

# save DBC
outRhs = outdir + 'dbc-g.txt'
fileRhs = open(outRhs, "w")
bufferRhs = ''
for I in range(n):
    bufferRhs += str(I) + ' ' + str(g.getValue(I)) + '\n'
fileRhs.write(bufferRhs)
fileRhs.close()
print('save rhs into ' + outRhs)       


# save A and S and RHS
PA = A.mat()
PM = M.mat()
PS = S.mat()
Pb = b.vec()

# save A
filenameT = outdir + 'mat-K.txt'
fileT = open(filenameT, "w")
bufferT = ''
for I in range(n):
    cols,vals = PA.getRow(I)
    for cj in range(len(cols)):
        bufferT += str(I) + ' ' + str(cols[cj]) + ' ' + str(vals[cj]) + '\n'
fileT.write(bufferT)
fileT.close()
print('save mat A into ' + filenameT)

# save M
filenameM = outdir + 'mat-M.txt'
fileM = open(filenameM, "w")
bufferM = ''
for I in range(n):
    cols,vals = PM.getRow(I)
    for cj in range(len(cols)):
        bufferM += str(I) + ' ' + str(cols[cj]) + ' ' + str(vals[cj]) + '\n'
fileM.write(bufferM)
fileM.close()
print('save mat M into ' + filenameM)

# save S
filenameS = outdir + 'mat-S.txt'
fileS = open(filenameS, "w")
bufferS = ''
for I in range(n):
    cols,vals = PS.getRow(I)
    for cj in range(len(cols)):
        bufferS += str(I) + ' ' + str(cols[cj]) + ' ' + str(vals[cj]) + '\n'
fileS.write(bufferS)
fileS.close()
print('save mat S into ' + filenameS)

# save Rhs
outRhs = outdir + 'rhs.txt'
fileRhs = open(outRhs, "w")
bufferRhs = ''
for I in range(n):
    sval = Pb.getValue(I)
    bufferRhs += str(I) + ' ' + str(sval) + '\n'
fileRhs.write(bufferRhs)
fileRhs.close()
print('save rhs into ' + outRhs)        


