using PyPlot
using LinearAlgebra
using SparseArrays
using Arpack
using ForwardDiff
using BenchmarkTools
using Formatting
using FileIO




#readfile
println("opening file... \n")
a = open("/home/billykon/programing/2d-elasticity-quadrilateral/plaka5x1.unv","r")
include("readfile.jl") #returns array of: nodes, edge, elem, bc
#include("mesh.jl")
include("elastic_functions.jl") #defines basic functions

#start building global array
n = length(elem[:,1])
dofs = 2*length(nodes[:,1])
global E = Ematrix(220e9,0.5,"stress")
global K = sparse([],[],Float64[],dofs,dofs);
#global M = sparse([],[],Float64[],dofs,dofs);

println("calculating stiffness/mass matrix...\n")
for i =1:n
    i1 = elem[i,:]'
    xs = nodes[i1,1]
    ys = nodes[i1,2]
    i1 = [i1 (i1 .+ Int64(0.5*dofs))]
    kre = ke(E,xs,ys)
    #mre = 7e3 .* me(xs,ys)
    K[i1,i1] = K[i1,i1][1,:,1,:] + kre
    #M[i1,i1] = M[i1,i1][1,:,1,:] + mre
end

uv = [.5 .5]




#start boundary conditions
println("calculating right-hand-side...\n")
F = zeros(dofs,1)

global f(x) = [[1e6 0]',[0 0]']

(nbc,ni,ix) = IDs(bc) #nbc has number of conditions
                   #ix has edges nodes ids for each condition

global bc,edge,elem,nodes,i,F
i=1

while true
   global i,ix,elem
   if i>nbc
       break
   end
   rows = [j for k=1:length(ni[i]) for j=1:length(elem[:,1]) if issubset(ix[i][k,1],elem[j,:]) &
                                                                issubset(ix[i][k,2],elem[j,:]) &
                                                                issubset(ix[i][k,3],elem[j,:])] #contains superset id element
   xedge = nodes[elem[rows,:],:][:,:,1]; yedge = nodes[elem[rows,:],:][:,:,2] #contains x,y of edges

   superelem = elem[rows,:]
   mapx = [sum(dims=2,superelem[k,:] .== ix[i][k,:]') for k=1:length(superelem[:,1])]#zero the non edge indexes of element
   for k=1:length(mapx)
       g(x) = f(x)[i]
       a = BC(g,xedge[k,:]',yedge[k,:]',mapx[k])
       id = [ix[i][k,:]' (ix[i][k,:]' .+ Int64(0.5*dofs))]
       id = sort(id,dims = 2)

       idlocal = [mapx[k]' mapx[k]'] #has the global degrees of freedoms
       idlocal = [idlocal[i] .==1 for i=1:18]
       F[id] .= F[id] .+ a[idlocal]'
   end
   i = i+1
end
#gib traction nodes where (u,v) are known
#reminder, ix contains boundary condition nodes INDEXES
println("solving for U....\n")
freedofs = zeros(dofs,1)
fixedofs = ix[2][:]'#take first's b.c. nodes

fixedofs = [fixedofs 17]
freedofs[fixedofs] .= 1
freedofs = [i for i=1:dofs if freedofs[i]==0]






xs = nodes[:,1]'
ys = nodes[:,2]'


U = zeros(18,1)

Ufree = K[freedofs,freedofs]\F[freedofs]
U[freedofs] .= Ufree

σ = E*LN([.5 .5],xs,ys)*U

figure(11)
for j=0:0.1:.5
    for i=-1:0.01:1

    uv = [i j]
    (x,y)=quadra(uv,xs,ys)
    s = E*LN(uv,xs,ys)*U
    scatter(i,s[1])
end
end
