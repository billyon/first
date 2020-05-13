using Plots
using LinearAlgebra
using Polynomials

using SparseArrays
using Arpack

using ForwardDiff

using BenchmarkTools

#using SymPy
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
a = 20e-2;
E = 220e9.*ones(N,1); I = a^4 .*ones(N,1)./12; A = a^2 .*ones(N,1) #is square with axa
rho = 7.8e3;



#orismos ton desmeumenon komvon
BC = zeros(dofs,1)
#ignored = 1:3:dofs;#ignore turn of axis=ignore freedom
BC[1:3] .= 1; #BC[ignored].=1;
fixedofs = findall(x->x==1,BC[:])
freedofs = findall(x->x==0,BC[:])
#boundary conditions
U = zeros(dofs+3*N,1);
F = zeros(dofs+3*N,1); F[[4 5]] .= 10e3
#euresi olon ton deikton pou tha prospelasoun ton K
alldofs = zeros(N,6);
#orismos olikou mitroou stivarotitas
global K = sparse([],[],Float64[],dofs+3*N,dofs+3*N);
global M = sparse([],[],Float64[],dofs+3*N,dofs+3*N);
ignored = []

#build hermite interpolation
function Lagrange(x,xi,xj)
           p = poly(xj[:])
           p = p/p(xi)
           y = p(x)
           return y
end

function H(x,xi,xj)
    Lval(y) = Lagrange(y,xi,xj)
    Ldot = ForwardDiff.derivative(Lval,xi)
    H = (1-2*(x-xi)*Ldot)*Lagrange(x,xi,xj)^2
    return H
end
function Hbar(x,xi,xj)
    Hbar = (x-xi)*Lagrange(x,xi,xj)^2
    return Hbar
end
#define integration functions
function gauss(n)
    w = zeros(n,1)
    global coeffs = zeros(n+1,1);
    for k = Int64.(0:floor(n/2))
        coeffs[2*k+1] = (-1)^k*binomial(n,k)*binomial(2*n-2*k,n)/(2^n);
    end
    coeffs = Poly(reverse(coeffs,dims=1)[:])
    points = roots(coeffs);

    for i=1:n
        L = poly(points[[1:i-1; i+1:end]]);
        L = L/polyval(L,points[i]);
        L = polyint(L);
        w[i] = polyval(L,1)-polyval(L,-1);

    end
    #I = sum(w(i).*f(points));
    return w,points
end
gaussw,gausspoints = gauss(6) #used to integrate by gaussintegral

function gaussintegral(x1,x2,w,points,f)
    g(y) = 0.5*(x2-x1)*f.((x2-x1)*y/2 +(x2+x1)/2)
    I = sum(w.*g.(points))
    return I
end
#global K = zeros(dofs,dofs);
Nt(x) = [Lagrange(x,0,[0.5*L L]);
         Lagrange(x,0.5*L,[0 L]);
         Lagrange(x,L,[0 0.5*L])]
Ntd(x) = ForwardDiff.derivative(Nt,x)

Nb(x) = [H(x,0,[0.5*L L]);
         Hbar(x,0,[0.5*L L]);
         H(x,0.5*L,[0 L]);
         Hbar(x,0.5*L,[0 L]);
         H(x,L,[0 0.5*L]);
         Hbar(x,L,[0 0.5*L]) ]
Nbd(x) = ForwardDiff.derivative(Nb,x)
smolKt(x) = Ntd(x)*Ntd(x)'; smolKb(x) = Nb(x)*Nb(x)';
smolMt(x) = Nt(x)*Nt(x)'; smolMb(x) = Nbd(x)*Nbd(x)';

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
    Bbend = [-s c 0  0 0 0  0 0 0;
              0 0 1  0 0 0  0 0 0;
              0 0 0 -s c 0  0 0 0;
              0 0 0  0 0 1  0 0 0;
              0 0 0  0 0 0 -s c 0;
              0 0 0  0 0 0  0 0 1];
    Btens = [c s 0 0 0 0;
             0 0 c s 0 0;
             0 0 0 0 c s];
    kbend = E[i]*I[i]/L^3;

    #insert tensile
    #=
    ktens = E[i]*A[i]/L
    global Ketens = (ktens/3)*[ 7 -8  1;
                               -8 16 -8;
                                1 -8  7];
    Ketens = Btens'*Ketens*Btens
    global Metens = (rho*A[i]*L/15)*[2  1 -1/2;
                                     1  8  1;
                                    -1/2  1  2];
    Metens = Btens'*Metens*Btens;


    Kebend = (kbend/35/(12+13*L)^2).*[4*(25776+55560*L+30695*L^2) 2*L*(80592+173240*L+95025*L^2)   -512*(927+2040*L+1085*L^2) 512*L*(537+1175*L+665*L^2) -4*(56016+127800*L+75985*L^2) 2*L*(18384+42680*L+25585*L^2);
                     0  L^2*4*(118080+25400*L+13835*L^2) -256*L*(492+1040*L+555*L^2) 128*L^2*(354+770*L+435*L^2) -2*L*(17616+40120*L+23985*L^2) 2*L^2*(2832+6600*L+3995*L^2);
                     0  0     8192*(126+270*L+145*L^2)  1024*L*(12+40*L+25*L^2) -512*(1044+2280*L+1235*L^2) 256*L*(516+1120+605*L^2);
                     0  0     0   1024*L^2*(183+405*L+230*L^2) -512*L*(561+1255*L+715*L^2) 128*L^2*(378+850*L+485*L^2);
                     0  0     0     0   4*(189648+419640*L+234065*L^2) -2*L*(84432+186040*L+103025*L^2);
                     0  0     0     0  0   L^2*4*(12192+26680*L+14635*L^2)];
    Kebend = Symmetric(Kebend);
    Kebend = Bbend'*Kebend*Bbend;
    Mebend = rho*A[i]*(L/(12+13*L)^2)*[(337248+750276*L+417883*L^2)/15015   L*(57552+132764*L+76767*L^2)/45045  8*(18408+41602*L+23423*L^2)/30030     -2*L*(31848+60578*L+27535*L^2)/45045 (81168+181304*L+10061*L^2)/30030 -L*(48720+104128*L+54703*L^2)/45045;
                                       0   L^2*(4236+10167*L+6131*L^2)/45045 8*L*(5772+13573*L+7907*L^2)/45045 -2*L^2*(764+1381*L+556*L^2)/15015 L*(58704+139712*L+81263*L^2)/180180 -L^2*(1840+4208*L+2343*L^2)/60060;
                                       0     0   1024*(858+1846+993*L^2)/15015    128*L*(156+556*L+415*L^2)/45045 8*(15912+32706*L+16783*L^2)/15015 -8*L*(4524+9125*L+4587*L^2)/45045;
                                       0     0    0     128*L^2*(471+1029*L+568*L^2)/45045 2*L*(42456+98782*L+56745*L^2)/45045 -2*L^2*(1076+2515*L+1441*L^2)/15015;
                                       0     0    0    0     (316032+673868*L+359463*L^2)/15015 -L*(49440+103456*L+54197*L^2)/45045;
                                       0     0    0    0    0    L^2*(300+615*L+316*L^2)/4095];
    Mebend = Symmetric(Mebend);
    Mebend = Bbend'*Mebend*Bbend;
    =#
    ktens = E[i]*A[i]
    kbend = E[i]*I[i]

    Ketens = ktens*gaussintegral(0,L,gaussw,gausspoints,smolKt)
    Ketens = Btens'*Ketens*Btens
    Metens = rho[i]*A[i]*gaussintegral(0,L,gaussw,gausspoints,smolMt)
    Metens = Btens'*Metens*Btens

    Kebend = kbend.*gaussintegral(0,L,gaussw,gausspoints,smolKb)
    Kebend = Bbend'*Kebend*Bbend
    Mebend = rho[i]*A[i].*gaussintegral(0,L,gaussw,gausspoints,smolMb)
    Mebend = Bbend'*Mebend*Bbend



    tensiledofs = [1 2 4 5 7 8];
    global Ke = Kebend .+0;
    Ke[tensiledofs,tensiledofs] = Ke[tensiledofs,tensiledofs][1,:,1,:]  + Ketens;
        #println(size(Ke),"\n")
    ending = dofs+3*N
    #edofs = [3*i1-2 3*i1-1 3*i1 ending-3*i+3 ending-3*i+2 ending-3*i+1 3*i2-2 3*i2-1 3*i2]
    edofs = [1 2 3  7 8 9  4 5 6]
    push!(ignored,ending-3*i+3, ending-3*i+2, ending-3*i+1)
    K[edofs,edofs] = K[edofs,edofs][1,:,1,:] + Ke;
    #println(K,"\n")
    Me = Mebend .+0;
    Me[tensiledofs,tensiledofs] = Me[tensiledofs,tensiledofs][1,:,1,:]  + Metens;
    M[edofs,edofs] = M[edofs,edofs][1,:,1,:] + Me;
end
Kfinal = K[vcat(freedofs,ignored),vcat(freedofs,ignored)]
Mfinal = M[vcat(freedofs,ignored),vcat(freedofs,ignored)]
rhs = F
rhs = rhs[vcat(freedofs,ignored)]
rhs = rhs - K[vcat(freedofs,ignored),fixedofs]*U[fixedofs]

Ufree = Kfinal\rhs
U[freedofs] = Ufree[freedofs]
#Ufr = reshape(U,(3,nodes))
show(U[:])
#println(varinfo())
f = inv(Matrix(Mfinal))*Matrix(Kfinal)
omega = sqrt.(abs.(eigvals(f)))/(2*pi)
st = eigvecs(f)
plot(omega)

#calculate
