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
a = open("../programing/efelkismos-simple.unv","r")
#include("readfile.jl") #returns array of: nodes, edge, elem, bc
include("example-mesh.jl")
include("elastic_functions.jl") #defines basic functions

#start building global array
n = length(elem[:,1])
dofs = 2*length(nodes[:,1])
global E = Ematrix(220e9,0.5,"stress")
global K = sparse([],[],Float64[],dofs,dofs);
global M = sparse([],[],Float64[],dofs,dofs);

println("calculating stiffness/mass matrix...\n")
for i =1:n
    i1 = elem[i,:]'
    xs = nodes[i1,1]
    ys = nodes[i1,2]
    i1 = [i1 (i1 .+ Int64(0.5*dofs))]
    kre = ke(E,xs,ys)
    mre = 7e3 .* me(xs,ys)
    K[i1,i1] = K[i1,i1][1,:,1,:] + kre
    M[i1,i1] = M[i1,i1][1,:,1,:] + mre
end

include("rhs.jl") #calculate F

#gib traction nodes where (u,v) are known
#reminder, ix contains boundary condition nodes INDEXES
println("solving for U....\n")
freedofs = zeros(dofs,1)
fixedofs = ix[1][:]'#take first's b.c. nodes

fixedofs = [fixedofs[1] (fixedofs .+Int64(0.5*dofs))]
freedofs[fixedofs] .= 1
freedofs = [i for i=1:dofs if freedofs[i]==0]

U = zeros(dofs,1)
rhs = F
rhs = rhs[freedofs] - K[freedofs,fixedofs'][:,:,1]*U[fixedofs']
Kfree = K[freedofs,freedofs]
Mfree = M[freedofs,freedofs]

Ufree = U[freedofs]
U[freedofs] = Kfree\rhs

n = length(elem[:,1])

println("calculating element average stress....\n")
#calculate element stress
include("stress.jl")
#println(stress')
println("DONE! \n")

i1 = elem[2,:]'
xs = nodes[i1,1]
ys = nodes[i1,2]

println(E*LN([0 0 1],xs,ys)*U[[ i1 (i1 .+ Int64(0.5*dofs))]]')
display(stress[:,1:3:18])
xs = [1 0 0 0.5 0 0.5]
ys = [0 -.5 .5 -.25 0 .25]
l = [1 0 0]

#=
F = zeros(12,1)
F[1] = 0.5*1e6
freedofs = ones(1,12)
fixedofs = [2 5 3 11]
freedofs[fixedofs] .= 0
freedofs = freedofs .==1
freedofs = freedofs[:]
kre = ke(E,xs,ys)

U  = zeros(12,1)
Ufree = kre[freedofs,freedofs]\F[freedofs]
U[freedofs] .= Ufree
println(E*LN([2 1 1]./4,xs,ys)*U)
scatter()
for l1 = 0.:0.01:1.
    l2 = 0.5*(1-l1)
    l3 = l2
    sigma = E*LN([l1 l2 l3],xs,ys)*U
    scatter!([l1],[sigma[1]/1e6],legend=false,label="sx")
    #scatter!([l1],[sigma[2]/1e6],legend=false,label="sy")
    #scatter!([l1],[sigma[3]/1e6],legend=false,label="txy")
end
scatter!()
=#
