include("bispectrum.jl")

grid_k = randn(ComplexF32, (257, 512, 512))
dk = 0.02
N = 10
L = 1000
kmax = 0.2
bispectrum(grid_k, dk, N, L, kmax)