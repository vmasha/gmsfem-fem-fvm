import numpy as np
import random, time
print("FVM system generator")

outdir = './data/modelF/'
kdir = './data/perm/'
cval = 1.0
fval = 0.0

dd = 8

# load k(x)
NS = 10
NN = 64
Nf = NN*NS
for ii in range(NS):
    locN = ii*NS + 0
    infile0 = kdir+'/k'+str(locN)+'.txt'
    with open(infile0) as f:
        lines_list = f.readlines()
        my_data0 = [float(val) for val in lines_list[1::2]]# with scipping
    arr0 = np.reshape(my_data0, (NN, NN))
    for jj in range(1, NS):
        locN = ii*NS + jj
        infileij = kdir+'/k'+str(locN)+'.txt'
        with open(infileij) as f:
            lines_list = f.readlines()
            my_dataij = [float(val) for val in lines_list[1::2]]# with scipping
        arrij = np.reshape(my_dataij, (NN, NN))
        arr0 = np.concatenate((arr0, arrij), axis=1)# concatinate col-wise
    if ii==0:
        arrK = arr0 # concatinate row-wise
    else:
        arrK = np.concatenate((arrK, arr0))# concatinate col-wise
print(arrK.min(), arrK.max())
Nx, Ny = arrK.shape
print(Nx, Ny)


# rescale k(x)
NN = int(NN/dd)
arrK = arrK[::dd, ::dd]
Nx, Ny = arrK.shape
print(Nx, Ny, NN)


# generate T and S
import sys, math
import petsc4py
from petsc4py import PETSc
petsc4py.init(sys.argv)

hh = 1.0/Nx; volK = hh*hh
n = Nx*Ny

T = PETSc.Mat().createAIJ([n, n], nnz=5)
M = PETSc.Mat().createAIJ([n, n], nnz=1)
S = PETSc.Mat().createAIJ([n, n], nnz=1)
for i in range(Nx):
    for j in range(Ny):
        I = i*Ny+j
        diagval = 0
        if j!=0:
            val = 2.0/(1.0/arrK[i,j-1]+1.0/arrK[i,j]); diagval += val
            T.setValue(I, I-1, -val)
        if j!=(Ny-1):
            val = 2.0/(1.0/arrK[i,j+1]+1.0/arrK[i,j]); diagval += val
            T.setValue(I, I+1, -val)
        if i!=0:
            val = 2.0/(1.0/arrK[i-1,j]+1.0/arrK[i,j]); diagval += val
            T.setValue(I, I-Ny, -val)
        if (i!=(Nx-1)):
            val = 2.0/(1.0/arrK[i+1,j]+1.0/arrK[i,j]); diagval += val
            T.setValue(I, I+Ny, -val)
        T.setValue(I, I, diagval)
        # M
        mval = cval*volK
        M.setValue(I, I, mval)
        # S for GMsFEM
        sval = arrK[i,j]*volK
        S.setValue(I, I, sval)
T.assemblyBegin()
T.assemblyEnd()
M.assemblyBegin()
M.assemblyEnd()
S.assemblyBegin()
S.assemblyEnd()
print('generate T, S and M')

q = PETSc.Vec().createSeq(n) 
for i in range(Nx):
    for j in range(Ny):
        I = i*Ny+j
        q.setValue(I, fval*volK)
print('generate Rhs') 


# save T and S for pres and RHS

# save T
filenameT = outdir + 'mat-K.txt'
fileT = open(filenameT, "w")
bufferT = ''
for I in range(n):
    cols,vals = T.getRow(I)
    for cj in range(len(cols)):
        bufferT += str(I) + ' ' + str(cols[cj]) + ' ' + str(vals[cj]) + '\n'
fileT.write(bufferT)
fileT.close()
print('save mat T into ' + filenameT)

# save M
filenameT = outdir + 'mat-M.txt'
fileT = open(filenameT, "w")
bufferT = ''
for I in range(n):
    cols,vals = M.getRow(I)
    for cj in range(len(cols)):
        bufferT += str(I) + ' ' + str(cols[cj]) + ' ' + str(vals[cj]) + '\n'
fileT.write(bufferT)
fileT.close()
print('save mat M into ' + filenameT)

# save S
filenameT = outdir + 'mat-S.txt'
fileT = open(filenameT, "w")
bufferT = ''
for I in range(n):
    cols,vals = T.getRow(I)
    for cj in range(len(cols)):
        bufferT += str(I) + ' ' + str(cols[cj]) + ' ' + str(vals[cj]) + '\n'
fileT.write(bufferT)
fileT.close()
print('save mat S into ' + filenameT)

# save Rhs
outRhs = outdir + 'rhs.txt'
fileRhs = open(outRhs, "w")
bufferRhs = ''
for I in range(n):
    sval = q.getValue(I)
    bufferRhs += str(I) + ' ' + str(sval) + '\n'
fileRhs.write(bufferRhs)
fileRhs.close()
print('save rhs into ' + outRhs)    

# DOF of DBC
hh = 1.0/Nx
mg = PETSc.Vec().createSeq(n) 
g = PETSc.Vec().createSeq(n) 
for i in range(Nx):
    for j in range(Ny):
        I = i*Ny+j
#         if (i*hh < hh/2 or j*hh < hh/2 or i*hh > 1.0-hh-hh/2 or j*hh > 1.0-hh-hh/2):
        if (j*hh < hh/2):
            mg.setValue(I, 1.0)
            g.setValue(I, 1.0)
        else:
            mg.setValue(I, 0.0)
            g.setValue(I, 0.0)
print('generate Rhs') 

# save DBC
outRhs = outdir + 'dbc.txt'
fileRhs = open(outRhs, "w")
bufferRhs = ''
for I in range(n):
    bufferRhs += str(I) + ' ' + str(mg.getValue(I)) + '\n'
fileRhs.write(bufferRhs)
fileRhs.close()
print('save dbc into ' + outRhs)      

outRhs = outdir + 'dbc-g.txt'
fileRhs = open(outRhs, "w")
bufferRhs = ''
for I in range(n):
    bufferRhs += str(I) + ' ' + str(g.getValue(I)) + '\n'
fileRhs.write(bufferRhs)
fileRhs.close()
print('save g into ' + outRhs)        


# save figs
from dolfin import *
import math

meshc = UnitSquareMesh(NS, NS)
Vc = FunctionSpace(meshc, 'DG', 0)
uc = Function(Vc)
uarrc = uc.vector().array()

mesh = UnitSquareMesh(Nx, Ny)
V = FunctionSpace(mesh, 'DG', 0)
u = Function(V)
uarr = u.vector().array()
print("functions fenics")

for i in range(Nx):
    for j in range(Ny):
        I = i*Ny + j
        val = arrK[i,j]
        uarr[2*I] = val
        uarr[2*I+1] = uarr[2*I]

u.vector().set_local(uarr)
filef = File(outdir+"results/k.pvd")
filef << u 
print('k saved')