
using LinearAlgebra
using Plots

m = 40 #parents
l = 3*m #childs
n = 10000 #DOFS
mr = 0.002*0 #mutation rate
par = 4 #parents
elites = 2
bl = -5.12 .*ones(1,n) #lower limits
bu = 5.12 .*ones(1,n) #upper limits
#1st generation
b = 10.24 .*rand(l,n) .-5.12
bnew = b
a = []
while true
    global b,bnew
    #evaluate
    #F = sum(b.^2,dims=2)
    F = 10*n .+ sum(b.^2 .-10 .*cos.(2*pi.*b),dims = 2)
    sort_ix = sortperm(F[:])
    sorted = F[sort_ix]
    #elitism
    i_elite = sort_ix[1:elites]
    #parent
    sorted = sorted./sum(F)
    for i = 2:l
        sorted[i] = sorted[i] + sorted[i-1] #make linear roulette
    end
    random = rand(m,1)
    selections = zeros(m,1)
    for i=1:m
        selections[i] = sort_ix[searchsortedfirst(sorted,random[i])] #select from roulette
    end
    selections = Int.(selections)
    # start crossover
    for i = 1:l
        indexOFthree = selections[rand(1:m,par)] #who will mate
        o = rand(par,n) #random interpolation constants
        pnew = o.*b[indexOFthree,:]./par
        pnew = sum(pnew,dims=1)
        if !(i in i_elite) #exclude elite
            bnew[i,:] = pnew #happy birthday
        end
    end
    #mutate
    for i=1:l
        if !(i in i_elite)
            for j=1:n
                if rand(1)[1]<mr
                    o = rand(1)[1]
                    bnew[i,j] = o*bl[j]+(1-o)*bu[j]
                end
            end
        end
    end
    b = bnew
    println(minimum(F))
    append!(a,minimum(F))
end
