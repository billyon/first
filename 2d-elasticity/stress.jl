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
    u = U[[2 .* i1 (2 .* i1 .- 1)]]'
    stress[i,1:3] = E*b(l[1,:],xs,ys)*u
    stress[i,4:6] = E*b(l[2,:],xs,ys)*u
    stress[i,7:9] = E*b(l[3,:],xs,ys)*u
    stress[i,10:12] = E*b(l[4,:],xs,ys)*u
    stress[i,13:15] = E*b(l[5,:],xs,ys)*u
    stress[i,16:18] = E*b(l[6,:],xs,ys)*u
    for k =1:6
        B = b(l[k,:],xs,ys)

        sigma = E*B*u
        mean = mean + sigma/6;
     end
     #stress[i,:] = mean
 end
