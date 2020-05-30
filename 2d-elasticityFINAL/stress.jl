stress = zeros(n,3*6)
for i=1:length(elem[:,1])
    l =  [1 0 0;
          0 1 0;
          0 0 1;
          0.5 0.5 0;
          0 0.5 0.5;
          0.5 0 0.5]

    i1 = elem[i,:]'
    xs = nodes[i1,1]
    ys = nodes[i1,2]
    mean = zeros(3,1)
    #u = reshape( U[[2 .* i1 (2 .* i1 .- 1)]],(2,6))
    u = U[[i1 (i1 .+ Int64(0.5*dofs))]]'
    stress[i,1:3] = E*LN(l[1,:],xs,ys)*u
    stress[i,4:6] = E*LN(l[2,:],xs,ys)*u
    stress[i,7:9] = E*LN(l[3,:],xs,ys)*u
    stress[i,10:12] = E*LN(l[4,:],xs,ys)*u
    stress[i,13:15] = E*LN(l[5,:],xs,ys)*u
    stress[i,16:18] = E*LN(l[6,:],xs,ys)*u
    for k =1:6
        B = LN(l[k,:],xs,ys)

        sigma = E*B*u
        mean = mean + sigma/6;
     end
     #stress[i,:] = mean
 end
