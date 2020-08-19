include("src/bispectrum_utilities.jl")
include("src/powerspectrum.jl")
include("src/bispectrum.jl")
include("test.jl")

grid_k = randn(ComplexF32, (33, 64, 64))
dk = 0.01
N = 10
L = 1000
kmax = 0.1
Bk = bispectrum(grid_k, dk, 20, L, kmax)
#write_powerspectrum(Pk, dk, "test.txt")
