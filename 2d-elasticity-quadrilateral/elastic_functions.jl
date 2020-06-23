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
function IDs(bc)
    #returns
    #nbc: number of boundary conditions
    #ix: ids of edges for each condition
    #ni: contains number of edges for each condition
    nbc = count(x->x==-1,bc)+1
    ix = findall(x->x==-1,bc)
    ix = [ix[i][1] for i=1:nbc-1]
    prepend!(ix,0)
    ix = [collect(ix[i-1]+1:ix[i]-1) for i=2:nbc]

    ni = [length(i) for i in ix]
    ix = [ [i for i in bc[ix[j]]] for j=1:nbc-1]
    ix = [edge[ix[i],:] for i=1:nbc-1]
    return Int64(nbc-1),ni,ix
end
function BC(f,xedge,yedge,map)
    uv = [0. 0.]
    if map[1]==1 & map[2]==1
        ind = 1
        uv[2] = -1
    elseif map[2]==1 & map[3]==1
        ind = 2
        uv[1] = 1
    elseif map[3]==1 & map[4]==1
        ind = 1
        uv[2] = 1
    else
        ind = 2
        uv[1] = -1
    end
    g = [-0.949108 -0.741531 -0.405845 0.405845 0.741531 0.949108 0.0]
    w = [0.129485 0.279705 0.38183 0.38183 0.279705 0.129485 0.417959]
    rhs = zeros(18,1)
    for i=1:7
        uv[ind] = g[i]
        (x,y) = quadra(uv,xedge,yedge)
        nuv = Nuv(uv)
        xdot = nuv[ind,:]'*xedge'
        ydot = nuv[ind,:]'*yedge'
        road = (xdot[1]^2+ydot[1]^2)^(0.5)
        rhs = rhs + w[i]*N(uv)'*f([x y])*road
    end
    return rhs
end
