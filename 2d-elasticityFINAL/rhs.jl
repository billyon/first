
#start boundary conditions
println("calculating right-hand-side...\n")

F = zeros(dofs,1)

(nbc,ni,ix) = IDs(bc) #nbc has number of conditions
                   #ix has edges nodes ids for each condition
global bc,edge,elem,nodes,i,F
i=1

while true
    global i,ix,elem
    if i>nbc
        break
    end
    rows = [j for k=1:length(ni[i]) for j=1:length(elem[:,1]) if issubset(ix[i][k,1],elem[j,:]) &
                                                                 issubset(ix[i][k,2],elem[j,:]) &
                                                                 issubset(ix[i][k,3],elem[j,:])] #contains superset id element
    xedge = nodes[elem[rows,:],:][:,:,1]; yedge = nodes[elem[rows,:],:][:,:,2] #contains x,y of edges

    superelem = elem[rows,:]
    mapx = [sum(dims=2,superelem[k,:] .== ix[i][k,:]') for k=1:length(superelem[:,1])]#zero the non edge indexes of element
    for k=1:length(mapx)
        g(x) = f(x)[i]
        a = BC(g,xedge[k,:]',yedge[k,:]',mapx[k])
        id = [ix[i][k,:]' (ix[i][k,:]' .+ Int64(0.5*dofs))]
        id = sort(id,dims=2)
        
        idlocal = [mapx[k]' mapx[k]'] #has the global degrees of freedoms
        idlocal = [idlocal[i] .==1 for i=1:12]
        F[id] .= F[id] .+ a[idlocal]'
    end
    i = i+1
end
