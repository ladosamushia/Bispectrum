include("src/bispectrum.jl")
include("src/powerspectrum.jl")

grid_k = randn(ComplexF32, (257, 512, 512))
dk = 0.01
N = 5
L = 1000
kmax = 0.05
Bk = bispectrum(grid_k, dk, N, L, kmax)
#Pk = power_spectrum(grid_k, dk, 20, L)
