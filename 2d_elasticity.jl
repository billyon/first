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

a = open("../programing/efelkismos.unv","r")
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
    push!(bc,[arr[2]])
    if length(arr)>5
        push!(bc,[arr[6]])
    end
    i = i+1
    if (i>=n)
        break
    end
end
    l = line[i]
    @show l
    l = split(l)
    arr = [parse(Int64,k) for k in l]
    if arr[1]==-1
        break
    end
    n = arr[end]
    n = div(n,2)+mod(n,2)
    line = line[i+1:end]
end


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
function BC(f,xedge,yedge,mapx,mapy)
    #f is vertical array
    #needs to multiplied by h
    #(fx,fy)
    if mapx[1]==0
        l1 = zeros(7,1)
        l2 = [-0.949108  -0.741531  -0.405845  0.405845  0.741531  0.949108  0.0]
        l3 = 1 .- l2
        w = [0.129485  0.279705  0.38183  0.38183  0.279705  0.129485  0.417959]
        #transform to 0~1
        l2 = 0.5 .*l2 .+ 0.5
        l3 = 1 .- l2
        w = 0.5 .* w
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
        bc = bc + N(l)*f([x,y])
    end
    return bc
end



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
    kre = 0.1 .* ke(E,xs,ys)
    mre = 0.1 .* 7e3 .* me(xs,ys)
    K[i1,i1] = K[i1,i1][1,:,1,:] + kre
    M[i1,i1] = M[i1,i1][1,:,1,:] + mre
end
end
bottomedge1 = [4 44 14]; bottomedge1 = [0 4 14 0 44 0]
bottomedge2 = [14 45 2]; bottomedge2 = [14 2 0 45 0 0]
topedge1 = [24 56 3]; topedge1 = [3 24 0 0 0 56]
topedge2 = [1 57 24]; topedge2 = [0 1 24 0 57 0]

b = length(bc[:,1])
edges = [i for k=1:b for i=1:n if issubset(bc[k,1],elem[i,:]) & issubset(bc[k,2],elem[i,:]) & issubset(bc[k,3],elem[i,:])]
