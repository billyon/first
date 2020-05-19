using Plots
using LinearAlgebra
#using Polynomials
using SparseArrays
using Arpack
using ForwardDiff
using BenchmarkTools
using Formatting
using FileIO
include("test.jl")
using lol
mutable struct L
    L1::Float64
    L2::Float64
    L3::Float64
end
mutable struct XYs
    Xs::Array
    Ys::Array
end

println("opening file... \n")
a = open("../programing/efelkismos-hard.unv","r")
line = readlines(a)
line = line[21:end]
nodes = Array{Float64}[]
#read nodes
println("reading nodes...\n")
global i=1
while true
    global i
    l = line[i]
    l = split(l)
    arr = [parse(Float64,k) for k in l]
    if length(arr)==1 & (arr[1]== -1.)
        break
    end
    i = i+2
    push!(nodes,arr)
end
line = line[i+2:end]
#read elements
println("reading elements....\n")
elem = Array{Int64}[]
edge = Array{Int64}[]

i = 1
while true

    global i
    l = line[i]
    l = split(l)
    arr = [parse(Float64,k) for k in l]
    if length(arr)==1 & (arr[1]==-1)
        break
    elseif arr[2]==22
        i = i+3
        l = line[i-1]
        l = split(l)
        arr = [parse(Float64,k) for k in l]
        push!(edge,arr)
    elseif arr[2]==42
        l = line[i+1]
        l = split(l)
        arr = [parse(Float64,k) for k in l]
        push!(elem,arr)
        i = i+2
    end

end
#read boundaries
println("reading boundaries...\n")

function splitToInt(l)
    l = split(l)
    l = [parse(Int64,k) for k in l]
    return l
end
i = i+3 #record 1 of 2467 format
n = splitToInt(line[i])[end] #shows how many elems are to read
line = line[i+2:end]
bc = Array{Int64}[]
global n = div(n,2)+mod(n,2)

while true
global i = 1
global line,n

println("splitting boundaries \n")
while true
    global i,line,n
    if (i>n)
        break
    end
    if tryparse(Float64, split(line[i])[1]) == nothing
        break
    end
    arr = splitToInt(line[i])
    push!(bc,[arr[2]])
    if length(arr)>5
        push!(bc,[arr[6]])
    end
    i = i+1
end
    push!(bc,[-1])
    arr = splitToInt(line[i])
    if arr[1]==-1
        break
    end
    n = arr[end]
    n = div(n,2)+mod(n,2)
    line = line[i+2:end]
end



#=
using RecursiveArrayTools
elem = VectorOfArray(elem); elem = convert(Array,elem)
edge = VectorOfArray(edge); edge = convert(Array,edge)
=#
println("making the arrays \n")
using Flux:batch
nodes = Transpose(batch(nodes))
edge = Transpose(batch(edge))
elem = Transpose(batch(elem))
elem = elem[:,[1 3 5 2 4 6]][:,1,:]
bc = Transpose(batch(bc))

xs = 10*[-1 0 0 -0.5 0 -0.5]
ys = 10*[0 0 1 0 0.5 0.5]
l = [1 0 0]

n = length(elem[:,1])
dofs = 2*length(nodes[:,1])
global E = Ematrix(220e9,0.33,"strain")
global K = sparse([],[],Float64[],dofs,dofs);
global M = sparse([],[],Float64[],dofs,dofs);

println("calculating stiffness/mass matrix...\n")
@time begin
for i =1:n
    i1 = elem[i,:]'
    xs = nodes[i1,1]
    ys = nodes[i1,2]
    i1 = [2 .*i1 (i1 .* 2 .-1)]
    kre = 0.1 .* ke(E,xs,ys)
    mre = 0.1 .* 7e3 .* me(xs,ys)
    K[i1,i1] = K[i1,i1][1,:,1,:] + kre
    M[i1,i1] = M[i1,i1][1,:,1,:] + mre
end
end
#start boundary conditions
println("calculating right-hand-side...\n")

F = zeros(dofs,1)
#find indexes of edge that are to be bounded
nbc = count(x->x==-1,bc)+1
ix = findall(x->x==-1,bc)
ix = [ix[i][1] for i=1:nbc-1]
prepend!(ix,0)
ix = [collect(ix[i-1]+1:ix[i]-1) for i=2:nbc]
#ix now contains all indexes, seperated so you can now iterate ix[i], for each i you get the list of edges that are on boundarycondition group i
global bc,edge,elem,nodes,i,F
i=1
#iterate through groups
f(x) = [[0 0]', [0 1e6]']
while true
    if i>nbc-1
        break
    end
    global bc,edge,elem,nodes,rows,xedge,yedge,mapx,i,F
    edix = bc[ix[i]] #contains edge id
    ed = edge[edix,:]#contains edge nodes id
    rows = [i for k=1:length(ix[i]) for i=1:length(elem[:,1]) if issubset(ed[k,1],elem[i,:]) & issubset(ed[k,2],elem[i,:]) & issubset(ed[k,3],elem[i,:])] #contains superset id element
    xedge = nodes[elem[rows,:],:][:,:,1]; yedge = nodes[elem[rows,:],:][:,:,2] #contains x,y of edges
    #zero the non edge indexes of element
    superelem = elem[rows,:]
    mapx = [sum(dims=2,superelem[i,:] .== ed[i,:]') for i=1:length(superelem[:,1])]
    for k=1:length(mapx)
        global bc,edge,elem,nodes,rows,xedge,yedge,mapx,i,F
        g(x) = f(x)[i]
        a = BC(g,xedge[k,:],yedge[k,:],mapx[k])
        id = [2 .*elem[rows[k],:] (2 .* elem[rows[k],:].-1)]
        F[id] = a
    end
    i = i+1
end
#gib traction nodes where (u,v) are known
#reminder, ix contains boundary condition nodes INDEXES
println("solving for U....\n")
freedofs = zeros(dofs,1)
fixedofs = edge[bc[ix[1]],:]
fixedofs = elem[[i for k=1:length(ix[1]) for i=1:length(elem[:,1]) if issubset(fixedofs[k,1],elem[i,:]) & issubset(fixedofs[k,2],elem[i,:]) & issubset(fixedofs[k,3],elem[i,:])],:][:] #contains superset id element

freedofs[[2 .*fixedofs' (2 .* fixedofs' .-1)]] .= 1
freedofs = [i for i=1:dofs if freedofs[i]==0]
fixedofs = [2 .* fixedofs' (2 .* fixedofs' .-1)]'
U = zeros(dofs,1)
rhs = F
rhs = rhs[freedofs] - K[freedofs,fixedofs][:,:,1]*U[fixedofs]
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
stress = zeros(n,3)
for i=1:length(elem[:,1])
    l =  [1 0 0;
          0 1 0;
          0 0 1;
          0.5 0.5 0;
          0 0.5 0.5;
          0.5 0 0.5]

    i1 = elem[i,:]'
    xs = nodes[i1,1]
    ys = nodes[i1,2]
    mean = zeros(3,1)
    for k =1:6
        epsilon = ε(l[k,:],xs,ys)
        sigma = E*sigma
        mean = mean + sigma/6;
     end
     stress[i,:] = mean
 end

#println(stress')
println("DONE! \n")
