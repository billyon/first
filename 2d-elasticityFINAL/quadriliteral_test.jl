using PyPlot
using LinearAlgebra
using SparseArrays
using Arpack
using ForwardDiff
using BenchmarkTools
using Formatting
using FileIO



function Ematrix(e, v, cond)
    if cond == "stress"
        E = (e / (1 - v * v)) * [
            1 v 0
            v 1 0
            0 0 (1 - v) / 2]
    elseif cond == "strain"
        E =
            (e * (1 - v) / ((1 + v) * (1 - 2 * v))) * [
                1 v / (v - 1) 0
                v / (v - 1) 1 0
                0 0 (1 - 2 * v) / (2 - 2 * v)]
    else
        println("no planestress or planestrain condition")
        exit()
    end
    return E
end
function quadra(uv,xs,ys)
    (u,v) = uv
    N1 = u*v*(u-1)*(v-1)/4; N5 = v*(1-u*u)*(v-1)/2;
    N2 = u*v*(u+1)*(v-1)/4; N6 = u*(u+1)*(1-v*v)/2;
    N3 = u*v*(u+1)*(v+1)/4; N7 = v*(1-u*u)*(v+1)/2;
    N4 = u*v*(u-1)*(v+1)/4; N8 = u*(u-1)*(1-v*v)/2;
            N9 = (1-u*u)*(1-v*v);
    x = dot(xs,[N1 N2 N3 N4 N5 N6 N7 N8 N9])
    y = dot(ys,[N1 N2 N3 N4 N5 N6 N7 N8 N9])
    return x,y
end
function Nuv(uv)
    (u,v) = uv
    Nu1 = (2*u-1)*(v-1)*v/4;  Nv1 = (u-1)*u*(2*v-1)/4;
    Nu2 = (2*u+1)*(v-1)*v/4;  Nv2 = (u+1)*u*(2*v-1)/4;
    Nu3 = (2*u+1)*(v+1)*v/4;  Nv3 = (u+1)*u*(2*v+1)/4;
    Nu4 = (2*u-1)*(v+1)*v/4;  Nv4 = (u-1)*u*(2*v+1)/4;
    Nu5 = -u*(v-1)*v;         Nv5 = -(u*u-1)*(2*v-1)/2;
    Nu6 = -(v*v-1)*(2*u+1)/2; Nv6 = -u*(u+1)*v;
    Nu7 = -u*(v+1)*v;         Nv7 = -(u*u-1)*(2*v+1)/2;
    Nu8 = -(v*v-1)*(2*u-1)/2; Nv8 = -u*(u-1)*v;
    Nu9 = 2*u*(v*v-1);        Nv9 = 2*v*(u*u-1);

    Nu = [Nu1 Nu2 Nu3 Nu4 Nu5 Nu6 Nu7 Nu8 Nu9]
    Nv = [Nv1 Nv2 Nv3 Nv4 Nv5 Nv6 Nv7 Nv8 Nv9]

    return [Nu; Nv]
end
uv = [.5 .5]
function Jacobian(uv,xs,ys)
    J = Nuv(uv)*[xs' ys']
    return J
end
function Nxy(uv,xs,ys)
    J = Jacobian(uv,xs,ys)
    nuv = Nuv(uv)
    nxy = J\nuv
    return nxy
end
function N(uv)
    (u,v) = uv
    N1 = u*v*(u-1)*(v-1)/4; N5 = v*(1-u*u)*(v-1)/2;
    N2 = u*v*(u+1)*(v-1)/4; N6 = u*(u+1)*(1-v*v)/2;
    N3 = u*v*(u+1)*(v+1)/4; N7 = v*(1-u*u)*(v+1)/2;
    N4 = u*v*(u-1)*(v+1)/4; N8 = u*(u-1)*(1-v*v)/2;
            N9 = (1-u*u)*(1-v*v);

    N = [N1 N2 N3 N4 N5 N6 N7 N8 N9 zeros(1,9);
              zeros(1,9)           N1 N2 N3 N4 N5 N6 N7 N8 N9]
    return N
end
function LN(uv,xs,ys)
    nxy = Nxy(uv,xs,ys)
    nx = nxy[1,:]'
    ny = nxy[2,:]'
    nx = [nx zeros(1,9);
          zeros(1,9) nx]
    ny = [ny zeros(1,9);
          zeros(1,9) ny]

    LN = [1 0; 0 0; 0 1]*nx + [0 0; 0 1; 1 0]*ny

    return LN
end
function ke(E,xs,ys)
    g = [-0.949108 -0.741531 -0.405845 0.405845 0.741531 0.949108 0.0]
    w = [0.129485 0.279705 0.38183 0.38183 0.279705 0.129485 0.417959]
    kre = zeros(18,18)
    for i=1:7,j=1:7
        uv = [g[i] g[j]]
        ln = LN(uv,xs,ys)
        J = Jacobian(uv,xs,ys)
        kre = kre + w[i]*w[j]*ln'*E*ln*abs(det(J))
    end
    return kre
end

xs = [0 1 1 0 0.5 1 .5 0 .5].*100 .+ 1000
ys = [0 0 1 1 0 .5 1 .5 .5].*100 .+ 1000
xs = 2 .*[0 1 1 0 0.5 1 .5 0 .5] .+1000
ys = 2 .*[0 0 1 1 0 .5 1 .5 .5] .+ 1000
xs = 1 .*[0 1 1 0 0.5 1 .5 0 .5] .+ 0
ys = 1 .*[0 0 1 1 0 .5 1 .5 .5] .+ 0
E = Ematrix(220e9,0.5,"stress")
kre = ke(E,xs,ys)

freedofs = ones(1,18)
fixedofs = [1 8 4 10 17 13]
freedofs[fixedofs] .= 0
freedofs = [i for i=1:18 if freedofs[i] .== 1]'

global f = [1e6 0]';
F = zeros(1,18)
global g = [-0.949108 -0.741531 -0.405845 0.405845 0.741531 0.949108 0.0]
global w = [0.129485 0.279705 0.38183 0.38183 0.279705 0.129485 0.417959]
global rhs = zeros(18,1)
for i=1:7
    global rhs,f,g,w
    uv = [1 g[i]]
    nuv = Nuv(uv)
    xdot = nuv[2,:]'*xs'
    ydot = nuv[2,:]'*ys'
    road = (xdot[1]^2+ydot[1]^2)^(0.5)
    rhs = rhs + w[i]*N(uv)'*f*road
end
F[[2 3 6]]   .= rhs[[2 3 6]]

U = zeros(18,1)

Ufree = kre[freedofs,freedofs][1,:,1,:]\F[freedofs]'
U[freedofs] .= Ufree'

σ = E*LN([.5 .5],xs,ys)*U

figure(11)
for j=0:0.1:.5
    for i=-1:0.01:1

    uv = [i j]
    (x,y)=quadra(uv,xs,ys)
    s = E*LN(uv,xs,ys)*U
    scatter(x,s[1])
end
end
for k=1:1
    figure(k+300)
    for i=-1:0.1:1
        for j= -1:0.1:1
            n = Jacobian([i j],xs,ys)
            #f = n[1,k]
            f = det(n)
            #f = Nxy([i j],xs,ys)[1,k]

            scatter3D([i],[j],[f])
        end
    end
end
