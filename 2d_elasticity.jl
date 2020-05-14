using Plots
using LinearAlgebra
#using Polynomials
using SparseArrays
using Arpack
using ForwardDiff
using BenchmarkTools
using Formatting
using FileIO

mutable struct L
    L1::Float64
    L2::Float64
    L3::Float64
end
mutable struct XYs
    Xs::Array
    Ys::Array
end

<<<<<<< HEAD

a = open("../programing/efelkismos-hard.unv","r")
=======
a = open("../programing/efelkismos.unv","r")
>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d
line = readlines(a)
line = line[21:end]
nodes = Array{Float64}[]
#read nodes
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
elem = Array{Int64}[]
edge = Array{Int64}[]

i = 1
while true
<<<<<<< HEAD

=======
>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d
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
<<<<<<< HEAD
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

while true
    global i,line,n
    if (i>n)
        break
    end
    if tryparse(Float64, split(line[i])[1]) == nothing
        break
    end
    arr = splitToInt(line[i])
=======

i = i+3
n = split(line[i])
n = [parse(Int64,k) for k in n]
n = n[end]
line = line[i+2:end]
bc = Array{Int64}[]
global n = div(n,2)+mod(n,2)
while true
global i = 1
global line,n
while true
    global i,line,n
    if tryparse(Float64, split(line[i])[1]) == nothing
        break
    end
    l = line[i]
    l = split(l)
    arr = [parse(Int64,k) for k in l]
>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d
    push!(bc,[arr[2]])
    if length(arr)>5
        push!(bc,[arr[6]])
    end
    i = i+1
<<<<<<< HEAD
end
    push!(bc,[-1])
    arr = splitToInt(line[i])
=======
    if (i>=n)
        break
    end
end
    l = line[i]
    @show l
    l = split(l)
    arr = [parse(Int64,k) for k in l]
>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d
    if arr[1]==-1
        break
    end
    n = arr[end]
    n = div(n,2)+mod(n,2)
<<<<<<< HEAD
    line = line[i+2:end]
end



=======
    line = line[i+1:end]
end


>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d
#=
using RecursiveArrayTools
elem = VectorOfArray(elem); elem = convert(Array,elem)
edge = VectorOfArray(edge); edge = convert(Array,edge)
=#
using Flux:batch
nodes = Transpose(batch(nodes))
edge = Transpose(batch(edge))
elem = Transpose(batch(elem))
elem = elem[:,[1 3 5 2 4 6]][:,1,:]
<<<<<<< HEAD
bc = Transpose(batch(bc))
=======
>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d

xs = 10*[-1 0 0 -0.5 0 -0.5]
ys = 10*[0 0 1 0 0.5 0.5]
l = [1 0 0]


function tria6(l,xs,ys)
    (l1,l2,l3) = l
    N1 = (2*l1-1)*l1; N2 = (2*l2-1)*l2; N3 = (2*l3-1)*l3;
    N4 = 4*l1*l2; N5 = 4*l2*l3; N6 = 4*l3*l1;
    x = dot(xs,[N1 N2 N3 N4 N5 N6])
    y = dot(ys,[N1 N2 N3 N4 N5 N6])
    return x,y
end
function N(l)
    (l1,l2,l3) = l;
    N1 = (2*l1-1)*l1; N2 = (2*l2-1)*l2; N3 = (2*l3-1)*l3;
    N4 = 4*l1*l2; N5 = 4*l2*l3; N6 = 4*l3*l1;

    N = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0;
         0 N1 0 N2 0 N3 0 N4 0 N5 0 N6];
    return N
end
function J(l,xs,ys)
    (l1,l2,l3) = l
    Nl1 = [4*l1-1 0 0 4*l2 0 4*l3]
    Nl2 = [0 4*l2-1 0 4*l1 4*l3 0]
    Nl3 = [0 0 4*l3-1 0 4*l2 4*l1]
    Nl = [Nl1; Nl2; Nl3]
    J = Nl*[xs' ys']
end
function Nl(l)
    (l1,l2,l3) = l
    Nl1 = [4*l1-1 0 0 4*l2 0 4*l3]
    Nl2 = [0 4*l2-1 0 4*l1 4*l3 0]
    Nl3 = [0 0 4*l3-1 0 4*l2 4*l1]
    Nl = [Nl1; Nl2; Nl3]
end
function Nxy(l,xs,ys)
    (l1,l2,l3) = l
    j = J(l,xs,ys)
    nl = Nl(l)
    nxy = j\nl
    return nxy
end


function me(xs,ys)
    #neeeds multiplication by rho*h
    me  = zeros(12,12)
    a1 = 0.0597158717; b1 = 0.4701420641; a2 = 0.7974269853; b2 = 0.1012865073
    w1 = 0.2250000000; w2 =0.1323941527; w3=0.1323941527; w4=0.1323941527
    w5 =0.1259391805; w6 = 0.1259391805; w7 = 0.1259391805;
    g1 = [1 1 1]./3; g2 = [a1 b1 b1]; g3 = [b1 a1 b1]; g4 = [b1 b1 a1];
    g5 = [a2 b2 b2]; g6 = [b2 a2 b2]; g7 = [b2 b2 a2]
    w = [w1 w2 w3 w4 w5 w6 w7]
    g = [[g1] [g2] [g3] [g4] [g5] [g6] [g7]]

    for i=1:7
        n = N(g[i])
        me = me + w[i]*n'*n
    end
    return me
end
function Ematrix(e,v,cond)
    if cond=="stress"
        E = (e/(1-v^2))*[1 v 0;
                        v 1 0;
                        0 0 (1-v)/2]
    elseif cond=="strain"
        E = (e*(1-v)/((1+v)*(1-2*v)))*[1 v/(v-1) 0;
                                       v/(v-1) 1 0;
                                        0 0 (1-2*v)/(2-2*v)]
    else
        println("motherfucker no planestress or planestrain condition")
        exit()
    end
    return E
end
function ke(E,xs,ys)
    #neeeds multiplication by h
    me  = zeros(12,12)
    a1 = 0.0597158717; b1 = 0.4701420641; a2 = 0.7974269853; b2 = 0.1012865073
    w1 = 0.2250000000; w2 =0.1323941527; w3=0.1323941527; w4=0.1323941527
    w5 =0.1259391805; w6 = 0.1259391805; w7 = 0.1259391805;
    g1 = [1 1 1]./3; g2 = [a1 b1 b1]; g3 = [b1 a1 b1]; g4 = [b1 b1 a1];
    g5 = [a2 b2 b2]; g6 = [b2 a2 b2]; g7 = [b2 b2 a2]
    w = [w1 w2 w3 w4 w5 w6 w7]
    g = [[g1] [g2] [g3] [g4] [g5] [g6] [g7]]

    ke = zeros(12,12)
    for i=1:7
        nxy = Nxy(g[i],xs,ys)
        nx = nxy[1,:]; ny = nxy[2,:];
        nx = [1 0 0 0 0 0;0 0 0 0 0 0;0 1 0 0 0 0;0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0]*nx
        ny = [0 0 0 0 0 0;1 0 0 0 0 0;0 0 0 0 0 0;0 1 0 0 0 0;0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 0;0 0 0 0 0 1]*ny;
        b = [1 0;0 0;0 1]*[nx';ny'] + [0 0;0 1;1 0]*[nx';ny']

        ke = ke + w[i]*b'*E*b
    end
    return ke
end
<<<<<<< HEAD
function BC(f,xedge,yedge,map)
    #f is vertical array
    #needs to multiplied by h
    #(fx,fy)
    if map[1]==0
        l1 = zeros(7,1)
        l2 = [-0.949108  -0.741531  -0.405845  0.405845  0.741531  0.949108  0.0]
=======
function BC(f,xedge,yedge,mapx,mapy)
    #f is vertical array
    #needs to multiplied by h
    #(fx,fy)
    if mapx[1]==0
        l1 = zeros(7,1)
        l2 = [-0.949108  -0.741531  -0.405845  0.405845  0.741531  0.949108  0.0]
        l3 = 1 .- l2
>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d
        w = [0.129485  0.279705  0.38183  0.38183  0.279705  0.129485  0.417959]
        #transform to 0~1
        l2 = 0.5 .*l2 .+ 0.5
        l3 = 1 .- l2
        w = 0.5 .* w
<<<<<<< HEAD
    elseif map[2]==0
        l2 = zeros(7,1)
        l3 = [-0.949108  -0.741531  -0.405845  0.405845  0.741531  0.949108  0.0]
        w = [0.129485  0.279705  0.38183  0.38183  0.279705  0.129485  0.417959]
        #transform to 0~1
        l3 = 0.5 .*l3 .+ 0.5
        l1 = 1 .- l3
        w = 0.5 .* w
    elseif map[3]==0
        l3 = zeros(7,1)
        l1 = [-0.949108  -0.741531  -0.405845  0.405845  0.741531  0.949108  0.0]
=======
    elseif mapx[2]==0
        l2 = zeros(7,1)
        l3 = [-0.949108  -0.741531  -0.405845  0.405845  0.741531  0.949108  0.0]
        l3 = 1 .- l1
        w = [0.129485  0.279705  0.38183  0.38183  0.279705  0.129485  0.417959]
        #transform to 0~1
        l3 = 0.5 .*l3 .+ 0.5
        l1 = 1 .- l2
        w = 0.5 .* w
    elseif mapx[3]==0
        l3 = zeros(7,1)
        l1 = [-0.949108  -0.741531  -0.405845  0.405845  0.741531  0.949108  0.0]
        l2 = 1 .- l1
>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d
        w = [0.129485  0.279705  0.38183  0.38183  0.279705  0.129485  0.417959]
        #transform to 0~1
        l1 = 0.5 .*l1 .+ 0.5
        l2 = 1 .- l1
        w = 0.5 .* w
    else
        println("this is not a valid edge map")
        exit()
    end
    bc = zeros(12,1)
    for i=1:7
        l = (l1[i],l2[i],l3[i])
        (x,y) = tria6(l,xedge,yedge)
<<<<<<< HEAD
        bc = bc + N(l)'*f([x,y])
=======
        bc = bc + N(l)*f([x,y])
>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d
    end
    return bc
end



<<<<<<< HEAD


n = length(elem[:,1])
dofs = 2*length(nodes[:,1])
global E = Ematrix(220e9,0.33,"stress")
global K = sparse([],[],Float64[],dofs,dofs);
global M = sparse([],[],Float64[],dofs,dofs);

for i =1:n
    i1 = elem[i,:]'
    xs = nodes[i1,1]
    ys = nodes[i1,2]
    i1 = [2 .*i1 (i1 .* 2 .-1)]
    @show i1
=======
n = length(elem[:,1])
dofs = 12*n
global E = Ematrix(220e9,0.33,"stress")
global K = sparse([],[],Float64[],dofs,dofs);
global M = sparse([],[],Float64[],dofs,dofs);
@show @time begin
    for i =1:n

    i1 = elem[i,:]'
    xs = nodes[i1,1]
    ys = nodes[i1,2]
    i1 = [i1 (i1 .* 2)]
>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d
    kre = 0.1 .* ke(E,xs,ys)
    mre = 0.1 .* 7e3 .* me(xs,ys)
    K[i1,i1] = K[i1,i1][1,:,1,:] + kre
    M[i1,i1] = M[i1,i1][1,:,1,:] + mre
end
<<<<<<< HEAD

#start boundary conditions
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
    @show mapx[1]
    for k=1:length(mapx)
        global bc,edge,elem,nodes,rows,xedge,yedge,mapx,i,F
        g(x) = f(x)[i]
        a = BC(g,xedge[k,:],yedge[k,:],mapx[k])
        @show a
        id = [2 .*elem[rows[k],:] (2 .* elem[rows[k],:].-1)]
        F[id] = a
        @show F[id]
    end
    i = i+1
end
#gib traction nodes where (u,v) are known
#reminder, ix contains boundary condition nodes INDEXES
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


#calculate nodal stress
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
        nxy = Nxy(l,xs,ys)
        nx = nxy[1,:]; ny = nxy[2,:];
        nx = [1 0 0 0 0 0;0 0 0 0 0 0;0 1 0 0 0 0;0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0]*nx
        ny = [0 0 0 0 0 0;1 0 0 0 0 0;0 0 0 0 0 0;0 1 0 0 0 0;0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 0;0 0 0 0 0 1]*ny;
        nxy = [nx';ny']
        xy = [collect(zip(xs,ys))[i][k] for i=1:6 for k=1:2]
        epsilon = E*parag*xy
        mean = mean + epsilon/6;
     end
     stress[i,:] = mean
 end
=======
end
bottomedge1 = [4 44 14]; bottomedge1 = [0 4 14 0 44 0]
bottomedge2 = [14 45 2]; bottomedge2 = [14 2 0 45 0 0]
topedge1 = [24 56 3]; topedge1 = [3 24 0 0 0 56]
topedge2 = [1 57 24]; topedge2 = [0 1 24 0 57 0]

b = length(bc[:,1])
edges = [i for k=1:b for i=1:n if issubset(bc[k,1],elem[i,:]) & issubset(bc[k,2],elem[i,:]) & issubset(bc[k,3],elem[i,:])]
>>>>>>> 9326345e2427a87a121737451d7cb23681b8dc0d
