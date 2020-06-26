using PyPlot
using LinearAlgebra
using SparseArrays
using Arpack
using ForwardDiff
using BenchmarkTools
using Formatting
using FileIO
using DelimitedFiles



#readfile
println("opening file... \n")
a = open(
    "/home/billykon/programing/2d-elasticity-quadrilateral/plaka15x1.unv",
    "r",
)
b = open(
    "/home/billykon/programing/2d-elasticity-quadrilateral/plaka15x1.dat",
    "r",
)
include("readfile.jl") #returns array of: nodes, edge, elem, bc
#include("mesh.jl")
include("elastic_functions.jl") #defines basic functions

#start building global array
n = length(elem[:, 1])
dofs = 2 * length(nodes[:, 1])
v = .3
Y = 200e9
global E = Ematrix(Y, v, "stress")

global K = sparse([], [], Float64[], dofs, dofs);
#global M = sparse([],[],Float64[],dofs,dofs);

println("calculating stiffness/mass matrix...\n")
for i = 1:n
    i1 = elem[i, :]'
    xs = nodes[i1, 1]
    ys = nodes[i1, 2]
    i1 = [i1 (i1 .+ Int64(0.5 * dofs))]
    kre = ke(E, xs, ys)
    #mre = 7e3 .* me(xs,ys)
    K[i1, i1] = K[i1, i1][1, :, 1, :] + kre
    #M[i1,i1] = M[i1,i1][1,:,1,:] + mre
end

uv = [0.5 0.5]




#start boundary conditions
println("calculating right-hand-side...\n")
F = zeros(dofs, 1)

#---------------------------------------------------
P=1e6
global f(x) = [[0 P*(x[2]-x[2]^2)*6]', [0 0]']
#---------------------------------------------------
(nbc, ni, ix) = IDs(bc) #nbc has number of conditions
#ix has edges nodes ids for each condition
global bc, edge, elem, nodes, i, F
i = 1

while true
    global i, ix, elem
    if i > nbc
        break
    end
    rows = [j for k=1:ni[i] for j=1:size(elem)[1] if issubset(ix[i][k, 1], elem[j, :]) &&
               issubset(ix[i][k, 2], elem[j, :]) &&
               issubset(ix[i][k, 3], elem[j, :])] #contains superset id element
    xedge = nodes[elem[rows, :], :][:, :, 1]
    yedge = nodes[elem[rows, :], :][:, :, 2] #contains x,y of edges

    superelem = elem[rows, :]
    mapx = [sum(dims = 2, superelem[k, :] .== ix[i][k, :]')  for k = 1:length(superelem[:, 1])]#zero the non edge indexes of element
    for k = 1:length(mapx)
        g(x) = f(x)[i]
        a = BC(g, xedge[k, :]', yedge[k, :]', mapx[k])
        id = [ix[i][k, :]' (ix[i][k, :]' .+ Int64(0.5 * dofs))]
        id = sort(id, dims = 2)

        idlocal = [mapx[k]' mapx[k]'] #has the global degrees of freedoms
        idlocal = [idlocal[i] .== 1 for i = 1:18]
        F[id] .= F[id] .+ a[idlocal]'
    end
    i = i + 1
end
#gib traction nodes where (u,v) are known
#reminder, ix contains boundary condition nodes INDEXES
println("solving for U....\n")
freedofs = zeros(dofs, 1)
fixedofs = ix[2][:]'#take first's b.c. nodes

fixedofs = [fixedofs Int.(5 .+0.5*dofs)]
freedofs[fixedofs] .= 1
freedofs = [i for i = 1:dofs if freedofs[i] == 0]






xs = nodes[:, 1]'
ys = nodes[:, 2]'


U = zeros(dofs, 1)

Ufree = K[freedofs, freedofs] \ F[freedofs]
U[freedofs] .= Ufree

include("stress.jl")

writedlm("/home/billykon/programing/2d-elasticity-quadrilateral/stress.csv",stress,',')
writedlm("/home/billykon/programing/2d-elasticity-quadrilateral/stress.csv",fstress',',')
writedlm("/home/billykon/programing/2d-elasticity-quadrilateral/nodes.csv",nodes,',')
writedlm("/home/billykon/programing/2d-elasticity-quadrilateral/data.csv",["x" "y" "z" "sx" "sy" "txy" ; nodes fstress' ./ 1e6],',')
writedlm("/home/billykon/programing/2d-elasticity-quadrilateral/U.csv",U,',')

print(maximum(U))


#post processing of bending
#below analytical solution is compared to the found one
amplify = 50.
u = amplify .*U
nodex = nodes[:,1]
nodey = nodes[:,2]
nodedx = u[1:Int(0.5*dofs)].+nodex
nodedy = u[Int(0.5*dofs)+1:Int(dofs)].+nodey

figure(1)
scatter(nodedx,nodedy,marker="x")
axis("equal")

#analytical bending, I = 1/12, D=1, L=15
uexact = -(2*P/Y).*(nodey.-0.5).*(3. .*nodex.*(30. .-nodex).+2*v.*nodey.*(nodey.-1.))
vexact = (2*P/Y).*(nodex.^2 .*(45. .-nodex)+3*v*(15. .-nodex).*(nodey.-0.5).^2 +(4+5*v).*nodex./4.)
nodedxact = nodex .+ uexact.*amplify
nodedyact = nodey .+ vexact.*amplify

scatter(nodedxact,nodedyact,marker="o")
axis("equal")
