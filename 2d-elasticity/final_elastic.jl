using Plots
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
include("readfile.jl") #returns array of: nodes, edge, elem, bc

include("elastic_functions.jl") #defines basic functions

#start building global array
n = length(elem[:,1])
dofs = 2*length(nodes[:,1])
global E = Ematrix(220e9,0.5,"stress")
global K = sparse([],[],Float64[],dofs,dofs);
global M = sparse([],[],Float64[],dofs,dofs);

println("calculating stiffness/mass matrix...\n")
h =1
@time begin
for i =1:n
    i1 = elem[i,:]'
    xs = nodes[i1,1]
    ys = nodes[i1,2]
    i1 = [2 .*i1 (i1 .* 2 .-1)]
    i2 = sort(i1,dims=2)
    @show i2
    kre = h .* ke(E,xs,ys)
    mre = h .* 7e3 .* me(xs,ys)
    K[i1,i1] = K[i1,i1][1,:,1,:] + kre
    M[i1,i1] = M[i1,i1][1,:,1,:] + mre
end
end

include("rhs.jl") #calculate F

#gib traction nodes where (u,v) are known
#reminder, ix contains boundary condition nodes INDEXES
println("solving for U....\n")
freedofs = zeros(dofs,1)
fixedofs = edge[bc[ix[1]],:][:]'

freedofs[[2 .*fixedofs (2 .* fixedofs .-1)]] .= 1
freedofs = [i for i=1:dofs if freedofs[i]==0]
fixedofs = [2 .* fixedofs (2 .* fixedofs .-1)]
U = zeros(dofs,1)
rhs = F
rhs = rhs[freedofs] - K[freedofs,fixedofs'][:,:,1]*U[fixedofs']
Kfree = K[freedofs,freedofs]
Mfree = M[freedofs,freedofs]

Ufree = U[freedofs]
U[freedofs] = Kfree\rhs

n = length(elem[:,1])
for i=1:n
    i1 = elem[i,:]'
    xs = nodes[i1,1]
    ys = nodes[i1,2]
    i1 = [2 .*i1 (i1 .* 2 .-1)]
end

println("calculating element average stress....\n")
#calculate element stress
include("stress.jl")
#println(stress')
println("DONE! \n")

i1 = elem[1,:]'
xs = nodes[i1,1]
ys = nodes[i1,2]

println(E*b([1 1 1]./3,xs,ys)*U[[2 .* i1 (2 .* i1 .- 1)]]')

xs = [1 0 0 0.5 0 0.5]
ys = [0 -.5 .5 -.25 0 .25]
l = [1 0 0]
