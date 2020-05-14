using CUDAnative
using CuArrays

function applyfun(x::CuDeviceArray, F::CuDeviceArray, f)
    i = ( blockIdx().x - 1) * blockDim().x + threadIdx().x

    F[i] = f(x[i])
    return nothing
end

nthreads = 512
nblocks = 10
F = CuArrays.zeros(nthreads * nblocks)
x = CuArray( LinRange(0, 10, nthreads*nblocks) )

f(x, c) = c*x^2

#This workds
g1(x) = f(x,1)
@cuda blocks=nblocks threads=nthreads applyfun(x, F, g1 )

println("hey")