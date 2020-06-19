stress = zeros(n,3*9)
for i=1:length(elem[:,1])
    uv = [
    -1 -1
    -1 1
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
    mean = zeros(3,1)
    #u = reshape( U[[2 .* i1 (2 .* i1 .- 1)]],(2,6))
    u = U[[i1 (i1 .+ Int64(0.5*dofs))]]'
    stress[i,1:3] = E*LN(uv[1,:],xs,ys)*u
    stress[i,4:6] = E*LN(uv[2,:],xs,ys)*u
    stress[i,7:9] = E*LN(uv[3,:],xs,ys)*u
    stress[i,10:12] = E*LN(uv[4,:],xs,ys)*u
    stress[i,13:15] = E*LN(uv[5,:],xs,ys)*u
    stress[i,16:18] = E*LN(uv[6,:],xs,ys)*u
    stress[i,19:21] = E*LN(uv[7,:],xs,ys)*u
    stress[i,22:24] = E*LN(uv[8,:],xs,ys)*u
    stress[i,25:27] = E*LN(uv[9,:],xs,ys)*u
    for k =1:6
        B = LN(uv[k,:],xs,ys)

        sigma = E*B*u
        mean = mean + sigma/6;
     end
    #stress[i,:] = mean
 end
