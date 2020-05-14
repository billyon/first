using Plots
using LinearAlgebra

using SparseArrays
using Arpack
#komvos1 komvos2 komvos3
xy = [ 0 0;
       0 1]; #orismos komvon san xyz
nodes = size(xy)[2]; #arithmos komvon

#komvos1_start komvos1_end
numel = [1 2;


]; #orismos ton stixion os grammes me toys akraious komvous
N = size(numel)[1]; #arithmos stoixion
dofs = 3*nodes;
#orismos elastikotitas stixion
a = 10e-2;
E = 220e9.*ones(N,1); I = a^4 .*ones(N,1)./12; A = a^2 .*ones(N,1) #is square with axa
rho = 7.8e3;



#orismos ton desmeumenon komvon
BC = zeros(dofs,1)
#ignored = 1:3:dofs;#ignore turn of axis=ignore freedom
BC[1:3] .= 1; #BC[ignored].=1;
fixedofs = findall(x->x==1,BC[:])
freedofs = findall(x->x==0,BC[:])
#boundary conditions
U = zeros(dofs,1);
F = zeros(dofs,1); F[[5]] .= 10e3;
#euresi olon ton deikton pou tha prospelasoun ton K
alldofs = zeros(N,6);
#orismos olikou mitroou stivarotitas
global K = sparse([],[],Float64[],dofs,dofs);
global M = sparse([],[],Float64[],dofs,dofs);

#global K = zeros(dofs,dofs);
#start main loop to calculate stiffness matrix
for i=1:N
    i1 = numel[i,1]
    i2 = numel[i,2]
    #calclulate directions
    x1 = xy[1,i1]; x2 = xy[1,i2]; dx = x2-x1;
    y1 = xy[2,i1]; y2 = xy[2,i2]; dy = y2-y1;
    L = sqrt(dx^2+dy^2);
    nx = dx/L; ny = dy/L;
    #direction matrix
    c = nx; s = ny;
    Bbend = [-s c 0  0 0 0;
              0 0 1  0 0 0;
              0 0 0 -s c 0;
              0 0 0  0 0 1];
    Btens = [c s 0 0;
             0 0 c s];
    kbend = E[i]*I[i]/L^3;

    #insert tensile
    ktens = E[i]*A[i]/L
    global Ketens = ktens*[1 -1; -1 1]
    Ketens = Btens'*Ketens*Btens
    global Metens = rho*A[i]*L*[2 1; 1 2]./6;
    Metens = Btens'*Metens*Btens;


    Kebend = kbend*[12 6*L   -12   6*L;
                     0 4*L^2 -6*L 2*L^2;
                     0  0     12  -6*L;
                     0  0     0   4*L^2]
    Kebend = Symmetric(Kebend);
    Kebend = Bbend'*Kebend*Bbend;
    Mebend = rho*A[i]*(L/420)*[156   22*L  54     -13*L;
                                0   4*L^2 13*L   -3*L^2;
                                0     0   156     -22*L;
                                0     0    0     4*L^2];
    Mebend = Symmetric(Mebend);
    Mebend = Bbend'*Mebend*Bbend;

    tensiledofs = [1 2 4 5];
    Ke = Kebend .+0;
    Ke[tensiledofs,tensiledofs] = Ke[tensiledofs,tensiledofs][1,:,1,:]  + Ketens;
        #println(size(Ke),"\n")
    edofs = [3*i1-2 3*i1-1 3*i1  3*i2-2 3*i2-1 3*i2]
    K[edofs,edofs] = K[edofs,edofs][1,:,1,:] + Ke;
    #println(K,"\n")
    Me = Mebend .+0;
    Me[tensiledofs,tensiledofs] = Me[tensiledofs,tensiledofs][1,:,1,:]  + Metens;
    M[edofs,edofs] = M[edofs,edofs][1,:,1,:] + Me;
end
Kfinal = K[freedofs,freedofs]
Mfinal = M[freedofs,freedofs]
rhs = F
rhs = rhs[freedofs]
rhs = rhs - K[freedofs,fixedofs]*U[fixedofs]

Ufree = Kfinal\rhs
U[freedofs] = Ufree
U = reshape(U,(3,nodes))
show(U[:])
#println(varinfo())
omega = sqrt.(eigvals(inv(Matrix(Mfinal))*Matrix(Kfinal)))/(2*pi)
scatter(omega)

#calculate
