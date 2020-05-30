
#start boundary conditions
println("calculating right-hand-side...\n")

F = zeros(dofs,1)
#find indexes of edge that are to be bounded
nbc = count(x->x==-1,bc)+1
ix = findall(x->x==-1,bc)
ix = [ix[i][1] for i=1:nbc-1]
prepend!(ix,0)
ix = [collect(ix[i-1]+1:ix[i]-1) for i=2:nbc]
#ix now contains all indexes, seperated so you can now iterate ix[i], for each i you get the list of edges that are on boundarycondition group i
global bc,edge,elem,nodes,i,F
i=1
#iterate through groups
f(x) = [[0 0]', [0 1e6]']
while true
    if i>nbc-1
        break
    end
    global bc,edge,elem,nodes,rows,xedge,yedge,mapx,i,F
    edix = bc[ix[i]] #contains edge id
    ed = edge[edix,:]#contains edge nodes id
    rows = [i for k=1:length(ix[i]) for i=1:length(elem[:,1]) if issubset(ed[k,1],elem[i,:]) & issubset(ed[k,2],elem[i,:]) & issubset(ed[k,3],elem[i,:])] #contains superset id element
    xedge = nodes[elem[rows,:],:][:,:,1]; yedge = nodes[elem[rows,:],:][:,:,2] #contains x,y of edges
    #zero the non edge indexes of element
    superelem = elem[rows,:]
    mapx = [sum(dims=2,superelem[i,:] .== ed[i,:]') for i=1:length(superelem[:,1])]
    for k=1:length(mapx)
        global bc,edge,elem,nodes,rows,xedge,yedge,mapx,i,F
        g(x) = f(x)[i]
        edix = bc[ix[i]] #contains edge id
        ed = edge[edix,:]#contains edge nodes id
        a = h .* BC(g,xedge[k,:],yedge[k,:],mapx[k])
        id = [2 .*ed (2 .* ed .-1)]
        idlocal = [collect(zip(mapx[k].==1,mapx[k].==1))[i][j] for i=1:6 for j=1:2]
        F[id] .= F[id] .+ a[idlocal]'
    end
    i = i+1
end
