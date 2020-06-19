#defining functions
function tria6(l,xs,ys)
    (l1,l2,l3) = l
    N1 = (2*l1-1)*l1; N2 = (2*l2-1)*l2; N3 = (2*l3-1)*l3;
    N4 = 4*l1*l2; N5 = 4*l2*l3; N6 = 4*l3*l1;
    x = dot(xs,[N1 N2 N3 N4 N5 N6])
    y = dot(ys,[N1 N2 N3 N4 N5 N6])
    return x,y
end
function BC(f, xedge, yedge, map)
    #f is vertical array
    #(fx,fy)
    if map[1] == 0
        l1 = zeros(7, 1)
        l2 = [-0.949108 -0.741531 -0.405845 0.405845 0.741531 0.949108 0.0]
        w = [0.129485 0.279705 0.38183 0.38183 0.279705 0.129485 0.417959]
        #transform to 0~1
        l2 = 0.5 .* l2 .+ 0.5
        l3 = 1 .- l2
        w = 0.5 .* w
        road_coord = 3
        map = 1
    elseif map[2] == 0
        l2 = zeros(7, 1)
        l3 = [-0.949108 -0.741531 -0.405845 0.405845 0.741531 0.949108 0.0]
        w = [0.129485 0.279705 0.38183 0.38183 0.279705 0.129485 0.417959]
        #transform to 0~1
        l3 = 0.5 .* l3 .+ 0.5
        l1 = 1 .- l3
        w = 0.5 .* w
        road_coord =3
        map = 2
    elseif map[3] == 0
        l3 = zeros(7, 1)
        l1 = [-0.949108 -0.741531 -0.405845 0.405845 0.741531 0.949108 0.0]
        w = [0.129485 0.279705 0.38183 0.38183 0.279705 0.129485 0.417959]
        #transform to 0~1
        l1 = 0.5 .* l1 .+ 0.5
        l2 = 1 .- l1
        w = 0.5 .* w
        road_coord  = 2
        map = 3
    else
        println("this is not a valid edge map")
        exit()
    end


    bc = zeros(12, 1)
    for i = 1:7
        l = (l1[i], l2[i], l3[i])
        (x, y) = tria6(l, xedge, yedge)

        road = NL1L2L3(l,xedge,yedge) #gives [xl1 yl1; xl2 yl2; xl3 yl3]
        if road_coord == 1
            road = road[1,:]
        elseif road_coord ==2
            road = road[2,:]
        else
            road = road[3,:]
        end

        bc = bc + w[i]*N(l)' * f([x, y])*sqrt(road[1]^2+road[2]^2)
    end
    return bc
end


#SHITS ABOUT TO GET REAL
function N(l)
    (l1, l2, l3) = l
    N1 = (2 * l1 - 1) * l1
    N2 = (2 * l2 - 1) * l2
    N3 = (2 * l3 - 1) * l3
    N4 = 4 * l1 * l2
    N5 = 4 * l2 * l3
    N6 = 4 * l3 * l1

    N = [N1 N2 N3 N4 N5 N6 0 0 0 0 0 0;
         0 0 0 0 0 0 N1 N2 N3 N4 N5 N6]
    return N
end
function NL1L2L3(l,xs,ys)
    J = JacobianL1L2L3(l,xs,ys)
    J = [1 1 1; J']
    P = J\[0 0;1 0; 0 1]
    nxy = Nxy(l,xs,ys)'

    nl1l2l3 = P'\nxy'
    xyl1l2l3 = nl1l2l3*[xs' ys']
    return xyl1l2l3
end
function Nl(l)
#    to be removed


    (l1, l2, l3) = l
    Nl1 = [(4 * l1 - 1) 0 0 (4 * l2) 0 (4 * l3)]
    Nl2 = [0 (4 * l2 - 1) 0 (4 * l1) (4 * l3) 0]
    Nl3 = [0 0 (4 * l3 - 1) 0 (4 * l2) (4 * l1)]
    Nl = [Nl1; Nl2; Nl3]
    #returns partials of Ni at l
    return Nl
end
function JacobianL1L2L3(l,xs,ys)

    J = Nl(l)*[xs' ys']
    return J
end
#=
function Nxyl1l2l3(l,xs,ys)
    #row 1 has x derivatives
    #row 2 has y derivatives
    nl = Nl(l)
    J = JacobianL1L2L3(l,xs,ys)
    nxy = J\nl
    return nxy
end
=#
function Nxy(l,xs,ys)
    #row 1 has x derivatives
    #row 2 has y derivatives
    Bo = Nl1l2(l)
    J = JacobianL1L2(l,xs,ys)
    nxy = J\Bo
    return nxy
end
function LN(l,xs,ys)
    nxy = Nxy(l,xs,ys)
    nx = nxy[1,:]'
    ny = nxy[2,:]'
    nx = [nx zeros(1,6);
          zeros(1,6) nx]
    ny = [ny zeros(1,6);
          zeros(1,6) ny]

    LN = [1 0; 0 0; 0 1]*nx + [0 0; 0 1; 1 0]*ny

    return LN
end

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
function Nl1l2(l)
    (l1,l2,l3) = l
    nl1 = [(4*l1-1) 0 (-4*(1-l1-l2)+1) 4*l2 -4*l2 4*(1-l2-2*l1)]
    nl2 = [0 (4*l2-1) (-4*(1-l1-l2)+1) 4*l1 4*(1-l1-2*l2) -4*l1]
    return [nl1 ; nl2]
end
function JacobianL1L2(l,xs,ys)
    (l1,l2,l3) = l
    nl1 = [(4*l1-1) 0 (-4*(1-l1-l2)+1) 4*l2 -4*l2 4*(1-l2-2*l1)]
    nl2 = [0 (4*l2-1) (-4*(1-l1-l2)+1) 4*l1 4*(1-l1-2*l2) -4*l1]
    J = [nl1 ; nl2]*[xs' ys']
    return J
end
function ke(E,xs,ys)
    a1 = 0.0597158717
    b1 = 0.4701420641
    a2 = 0.7974269853
    b2 = 0.1012865073
    w1 = 0.2250000000
    w2 = 0.1323941527
    w3 = 0.1323941527
    w4 = 0.1323941527
    w5 = 0.1259391805
    w6 = 0.1259391805
    w7 = 0.1259391805
    g1 = [1 1 1] ./ 3
    g2 = [a1 b1 b1]
    g3 = [b1 a1 b1]
    g4 = [b1 b1 a1]
    g5 = [a2 b2 b2]
    g6 = [b2 a2 b2]
    g7 = [b2 b2 a2]
    w = [w1 w2 w3 w4 w5 w6 w7]
    g = [[g1] [g2] [g3] [g4] [g5] [g6] [g7]]

    ke = zeros(12, 12)
    for i = 1:7
        li = g[i]
        wi = w[i]
        jl1l2 = JacobianL1L2(li,xs,ys)
        ln = LN(li,xs,ys)
        ke = ke + 0.5*wi*ln'*E*ln*abs(det(jl1l2))
    end
    return ke
end

function me(xs,ys)
    a1 = 0.0597158717
    b1 = 0.4701420641
    a2 = 0.7974269853
    b2 = 0.1012865073
    w1 = 0.2250000000
    w2 = 0.1323941527
    w3 = 0.1323941527
    w4 = 0.1323941527
    w5 = 0.1259391805
    w6 = 0.1259391805
    w7 = 0.1259391805
    g1 = [1 1 1] ./ 3
    g2 = [a1 b1 b1]
    g3 = [b1 a1 b1]
    g4 = [b1 b1 a1]
    g5 = [a2 b2 b2]
    g6 = [b2 a2 b2]
    g7 = [b2 b2 a2]
    w = [w1 w2 w3 w4 w5 w6 w7]
    g = [[g1] [g2] [g3] [g4] [g5] [g6] [g7]]

    me = zeros(12, 12)
    for i = 1:7
        li = g[i]
        wi = w[i]
        jl1l2 = JacobianL1L2(li,xs,ys)
        ni = N(li)
        me = me + 0.5*wi*ni'*ni*abs(det(jl1l2))
    end
    return me
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
#=
xs = [1 0 0 0.5 0 0.5]
ys = [0 -.5 .5 -.25 0 .25]
E = Ematrix(220e9,0.5,"stress")
figure(1)
for l1=0:0.1:1
    for l2=0:0.1:1
        l3 = 1 -l1-l2
        if l3<0 break end
        l = [l1 l2 l3]
        (x,y) = tria6(l,xs,ys)
        f  = sum(N(l),dims=2)
        scatter3D(x,y,f[1])
    end
end
=#
xs = [1 0 0 0.5 0 0.5]
ys = [0 -.5 .5 -.25 0 .25]
