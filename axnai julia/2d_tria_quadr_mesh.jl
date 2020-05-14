using Plots
using LinearAlgebra
using Polynomials

using SparseArrays
using Arpack

using ForwardDiff

using BenchmarkTools

function tria6(l1,l2,l3,xs,ys)
    L1 = l1; L2=l2; L3 = l3;
    N1 = (2*L1-1)*L1; N2 = (2*L2-1)*L2; N3 = (2*L3-1)*L3;
    N4 = 4*L1*L2; N5 = 4*L2*L3; N6 = 4*L3*L1;
    x = dot(xs',[N1 N2 N3 N4 N5 N6])
    y = dot(ys',[N1 N2 N3 N4 N5 N6])
    return x,y
end
xs = 10*[-1 0.3 0 -.5 0.2 -.5]
ys = 10*[0 0 1 0.2 .5 .8]
scatter()
for l1 = 0:0.02:1
    for l2 = 0:0.02:1-l1
        scatter!(tria6(l1,l2,1-l1-l2,xs,ys))
    end
end
scatter!()

function tria6jacobian(l1,l2,l3,xs,ys)
    #calculates jacobian matrix dN1/dl1 ....
    L1 = l1; L2=l2; L3=l3;
    N1 = (2*L1-1)*L1; N2 = (2*L2-1)*L2; N3 = (2*L3-1)*L3;
    N4 = 4*L1*L2
    N5 = 4*L2*L3
    N6 = 4*L3*L1
    xy =[xs;ys];
    Ndot = [4*l1-1 0 0;
            0 4*l2-1 0;
            0 0 4*l3-1;
            4*l2 4*l1 0;
            0 4*l2 4*l3;
            4*l3 0 4*l1]
    jacobian = xy*Ndot;
    return jacobian
end
YM = 230e9; #YOUNGS MODULUS
ν = 0.33 # POISON RATIO
d11 = YM*(1-ν)/(1+ν)/(1-2*ν); d12 = ν*d11/(1-ν); d33 = YM/(2*(1+ν));
E = [d11 d12 0;
     d12 d11 0;
     0 0 d33]
function calcBEB(l1,l2,l3,xs,ys)
    L1 = l1; L2= l2; L3 = l3;
    jacobian = tria6jacobian(l1,l2,l3,xs,ys)
    JP = [1 1 1; jacobian] #according to p368 of Fellipa C. Colorado isop
    P = JP\[0 0; 1 0; 0 1]
    Ndot = [4*l1-1 0 0;
            0 4*l2-1 0;
            0 0 4*l3-1;
            4*l2 4*l1 0;
            0 4*l2 4*l3;
            4*l3 0 4*l1]
    Nderivxy = P'*Ndot'; # is [dN1/dx ... dN6/dx;
                        #     dN1/dy ... dN6/dy]
    N1 = (2*L1-1)*L1; N2 = (2*L2-1)*L2; N3 = (2*L3-1)*L3;
    N4 = 4*L1*L2
    N5 = 4*L2*L3
    N6 = 4*L3*L1
    N = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0;
         0 N1 0 N2 0 N3 0 N4 0 N5 0 N6]
    N1x = Nderivxy[1,1]; N2x = Nderivxy[1,2]; N3x = Nderivxy[1,3]; N4x = Nderivxy[1,4]; N5x = Nderivxy[1,5]; N6x = Nderivxy[1,6];
    N1y = Nderivxy[2,1]; N2y = Nderivxy[2,2]; N3y = Nderivxy[2,3]; N4y = Nderivxy[2,4]; N5y = Nderivxy[2,5]; N6y = Nderivxy[2,6];

    B = [N1x 0 N2x 0 N3x 0 N4x 0 N5x 0 N6x 0;
         0  N1y 0 N2y 0 N3y 0 N4y 0 N5y 0 N6y;
         N1y N1x N2y N2x N3y N3x N4y N4x N5y N5x N6y N6x]
    BEB = B'*E*B
    return BEB
end

function tria6Ke(xs,ys)
    Ke  = zeros(12,12)
    g1 = [0.108103, 0.445948, 0.445948]; w1 =0.223382;
    g2 = [0.445948, 0.108103, 0.445948]; w2 =0.223382;
    g3 = [0.445948, 0.445948, 0.108103]; w3 =0.223382;
    g4 = [0.816848, 0.0915762, 0.0915762]; w4 = 0.109952;
    g5 = [0.0915762, 0.816848, 0.0915762]; w5 = 0.109952;
    g6 = [0.0915762, 0.0915762, 0.816848]; w6 = 0.109952;
    g = [[g1] [g2] [g3] [g4] [g5] [g6]]; w = [w1 w2 w3 w4 w5 w6];

    for i=1:6
        Ke = Ke + w[i]*calcBEB(g[i][1],g[i][2],g[i][3],xs,ys)
    end
    return Ke
end

f = open("C:\\User\\billykon\\Desktop\\Mesh_1.unv","r")

function edgeToLs(xedge,yedge,TOPOLOGY):
