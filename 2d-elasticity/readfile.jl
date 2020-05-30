
line = readlines(a)
line = line[21:end]
nodes = Array{Float64}[]
#read nodes
println("reading nodes...\n")
global i=1
while true
    global i
    l = line[i]
    l = split(l)
    arr = [parse(Float64,k) for k in l]
    if length(arr)==1 & (arr[1]== -1.)
        break
    end
    i = i+2
    push!(nodes,arr)
end
line = line[i+2:end]
#read elements
println("reading elements....\n")
elem = Array{Int64}[]
edge = Array{Int64}[]

i = 1
while true

    global i
    l = line[i]
    l = split(l)
    arr = [parse(Float64,k) for k in l]
    if length(arr)==1 & (arr[1]==-1)
        break
    elseif arr[2]==22
        i = i+3
        l = line[i-1]
        l = split(l)
        arr = [parse(Float64,k) for k in l]
        push!(edge,arr)
    elseif arr[2]==42
        l = line[i+1]
        l = split(l)
        arr = [parse(Float64,k) for k in l]
        push!(elem,arr)
        i = i+2
    end

end
#read boundaries
println("reading boundaries...\n")

function splitToInt(l)
    l = split(l)
    l = [parse(Int64,k) for k in l]
    return l
end
i = i+3 #record 1 of 2467 format
n = splitToInt(line[i])[end] #shows how many elems are to read
line = line[i+2:end]
bc = Array{Int64}[]
global n = div(n,2)+mod(n,2)

while true
global i = 1
global line,n

println("splitting boundaries \n")
while true
    global i,line,n
    if (i>n)
        break
    end
    if tryparse(Float64, split(line[i])[1]) == nothing
        break
    end
    arr = splitToInt(line[i])
    push!(bc,[arr[2]])
    if length(arr)>5
        push!(bc,[arr[6]])
    end
    i = i+1
end
    push!(bc,[-1])
    arr = splitToInt(line[i])
    if arr[1]==-1
        break
    end
    n = arr[end]
    n = div(n,2)+mod(n,2)
    line = line[i+2:end]
end

println("making the arrays \n")
using Flux:batch
nodes = Transpose(batch(nodes))
edge = Transpose(batch(edge))
elem = Transpose(batch(elem))
elem = elem[:,[1 3 5 2 4 6]][:,1,:]
bc = Transpose(batch(bc))
