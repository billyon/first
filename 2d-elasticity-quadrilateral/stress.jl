stress = zeros(n,9,3)
for i=1:length(elem[:,1])
    uv = [
    -1 -1
    1 -1
    1 1
    -1 1
    0 -1
    1 0
    0 1
    -1 0
    0 0]

    i1 = elem[i,:]'
    xs = nodes[i1,1]
    ys = nodes[i1,2]

    u = U[[i1 (i1 .+ Int64(0.5*dofs))]]'
    stress[i,1,:] = E*LN(uv[1,:],xs,ys)*u
    stress[i,2,:] = E*LN(uv[2,:],xs,ys)*u
    stress[i,3,:] = E*LN(uv[3,:],xs,ys)*u
    stress[i,4,:] = E*LN(uv[4,:],xs,ys)*u
    stress[i,5,:] = E*LN(uv[5,:],xs,ys)*u
    stress[i,6,:] = E*LN(uv[6,:],xs,ys)*u
    stress[i,7,:] = E*LN(uv[7,:],xs,ys)*u
    stress[i,8,:] = E*LN(uv[8,:],xs,ys)*u
    stress[i,9,:] = E*LN(uv[9,:],xs,ys)*u
 end
