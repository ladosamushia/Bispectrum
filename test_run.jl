include("src/bispectrum_utilities.jl")
include("src/powerspectrum.jl")

grid_k = randn(ComplexF32, (33, 64, 64))
dk = 0.01
N = 10
L = 1000
kmax = 0.1
Pk = power_spectrum(grid_k, dk, 20, L)
write_powerspectrum(Pk, dk, "test.txt")