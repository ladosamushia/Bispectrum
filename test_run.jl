include("src/bispectrum_utilities.jl")
include("src/powerspectrum.jl")

<<<<<<< HEAD
grid_k = randn(ComplexF32, (257, 512, 512))
dk = 0.01
N = 10
L = 1000
kmax = 0.1
println(bispectrum_bins(10))
Bk = bispectrum(grid_k, dk, N, L, kmax)
write_bispectrum(Bk, dk, N, "test_output.txt")
#Pk = power_spectrum(grid_k, dk, 20, L)
