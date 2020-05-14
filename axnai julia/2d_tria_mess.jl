using Plots
using LinearAlgebra
using Polynomials

using SparseArrays
using Arpack

using ForwardDiff

using BenchmarkTools



#define integration functions
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
#produces (w,p)
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
#results I
function gaussintegral(x1,x2,w,points,f)
    g(y) = 0.5*(x2-x1)*f.((x2-x1)*y/2 +(x2+x1)/2)
    I = sum(w.*g.(points))
    return I
end

function quad(ksi,eta,xs,ys)
    #x and y must be line vectors!!!!
    #returns x,y on selected ksi, eta points
    n = length(xs)
    xs = reshape(xs,(n,n))'
    ys = reshape(ys,(n,n))'
    nξ = n; nη = n;
    ξ = LinRange(-1,1,nξ); η = LinRange(-1,1,nη)
    Nij(i,j,ξi,ηi) = Lagrange(ξi,ξ[i],ξ[deleteat!(collect(1:nξ),i)])*Lagrange(ηi,η[j],η[deleteat!(collect(1:nη),j)])
    x = 0
    y = 0
    for i=1:n
        for j=1:n
            x = x + xs[i,j]*Nij(i,j,ksi,eta)
            y = y + ys[i,j]*Nij(i,j,ksi,eta)
        end
    end
    return x,y
end
mutable struct Foo2
    x::Int64
end

function setx(a::Foo2, v)
 a.x=v
end
#nodesξη =reverse(collect(Iterators.product(ξ, η)),dims=1)
#=
ξη
(-1.0, -1.0)       (-1.0, -0.333333)       (-1.0, 0.333333)       (-1.0, 1.0)
(-0.333333, -1.0)  (-0.333333, -0.333333)  (-0.333333, 0.333333)  (-0.333333, 1.0)
(0.333333, -1.0)   (0.333333, -0.333333)   (0.333333, 0.333333)   (0.333333, 1.0)
(1.0, -1.0)        (1.0, -0.333333)        (1.0, 0.333333)        (1.0, 1.0)
=#
#=
ij
-------------------------------
|(1, 1)  (1, 2)  (1, 3)  (1, 4)|
|(2, 1)  (2, 2)  (2, 3)  (2, 4)|
|(3, 1)  (3, 2)  (3, 3)  (3, 4)|
|(4, 1)  (4, 2)  (4, 3)  (4, 4)|
-------------------------------
=#
