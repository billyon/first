
linea = readlines(a)
lineb = readlines(b)
function splitToInt(l)
    l = split(l)
    l = [parse(Int64,k) for k in l]
    return l
end
function splitToFloat(l)
    l = split(l)
    l = [parse(Float64,k) for k in l]
    return l
end
(n, e) = splitToInt(lineb[1])
lineb = lineb[2:end]
nodes = zeros(n,3)

for i=1:n
    l = splitToFloat(lineb[i])
    nodes[i,:] = [l[2] l[3] l[4]]
end
lineb = lineb[n+1:end]

elem = Array{Int64}[]
edge = Array{Int64}[]

for i=1:e
    arr = splitToInt(lineb[i])

    if arr[2]==103
        push!(edge,arr[3:end])
    elseif arr[2]==209
        push!(elem,arr[3:end])
    elseif length(arr)<=1
        break
    end
end


#read boundaries
println("reading boundaries...\n")
linea = linea[21:end]
global i = 1
while true
    global i
    l = linea[i]
    l = splitToFloat(l)
    if l[1] == 2467 && length(l)==1
        break
    end
    i = i+1
end
i = i+1#record 1 of 2467 format
linea = linea[i:end]

n = splitToInt(linea[1])[end] #shows how many elems are to read

line = linea[3:end]
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
#nodes = Transpose(batch(nodes))
edge = Transpose(batch(edge))
elem = Transpose(batch(elem))
bc = Transpose(batch(bc))
elem = elem[:,[2 3 4 1 6 7 8 5 9]][:,1,:]
