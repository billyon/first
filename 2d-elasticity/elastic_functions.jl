#defining functions
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
    a = [xs[2]-xs[1] ; ys[2]-ys[1] ; 0]
    c = [xs[3]-xs[1] ; ys[3]-ys[1] ; 0]
    ae = 0.5*norm(cross(a,c))
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
    return ae.*me
end
function Ematrix(e,v,cond)
    if cond=="stress"
        E = (e/(1-v*v))*[1 v 0;
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
    a = [xs[2]-xs[1] ; ys[2]-ys[1] ; 0]
    c = [xs[3]-xs[1] ; ys[3]-ys[1] ; 0]
    ae = 0.5*norm(cross(a,c))

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
    return ae.*ke
end
function BC(f,xedge,yedge,map)
    #f is vertical array
    #needs to multiplied by h
    #(fx,fy)
    if map[1]==0
        l1 = zeros(7,1)
        l2 = [-0.949108  -0.741531  -0.405845  0.405845  0.741531  0.949108  0.0]
        w = [0.129485  0.279705  0.38183  0.38183  0.279705  0.129485  0.417959]
        #transform to 0~1
        l2 = 0.5 .*l2 .+ 0.5
        l3 = 1 .- l2
        w = 0.5 .* w
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
        bc = bc + N(l)'*f([x,y])
    end
    return bc
end
function b(l,xs,ys)
    nxy = Nxy(l,xs,ys)
    nx = nxy[1,:]; ny = nxy[2,:];
    nx = [1 0 0 0 0 0;0 0 0 0 0 0;0 1 0 0 0 0;0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0]*nx
    ny = [0 0 0 0 0 0;1 0 0 0 0 0;0 0 0 0 0 0;0 1 0 0 0 0;0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 0;0 0 0 0 0 1]*ny;
    nxy = [nx';ny']
    xy = [collect(zip(xs,ys))[i][k] for i=1:6 for k=1:2]
    b = [1 0;0 0;0 1]*[nx';ny'] + [0 0;0 1;1 0]*[nx';ny']
    epsilon = b
    return epsilon
end
