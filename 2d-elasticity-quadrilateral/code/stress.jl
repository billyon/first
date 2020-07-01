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
println("calculating nodal stress...")
fstress = zeros(3, size(nodes,1), 2)
fstress[:,:,2] .== 1.
for i = 1:n
    a = elem[i,:]' #global node number
    mapa = [findall(x -> x == a[k], elem[i, :])[1] for k = 1:length(a)] #lobal node number

    fstress[:, a, 1] = fstress[:, a, 1][:, 1, :] .+ stress[i, mapa, :]'
    fstress[:, a, 2] = fstress[:, a, 2][:, 1, :] .+ 1
end
println("averaging stress...")
num_indexed = fstress[:,:,2] #node sharing over elements
fstress = fstress[:,:,1] ./ fstress[:,:,2] #divide to obtain average
stress =1
elem = 1
edges = 1
bc = 1
ix = 1
F = 1
print("Averaging complete!")
