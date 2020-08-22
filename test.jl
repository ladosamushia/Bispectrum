include("src/bispectrum.jl")
include("src/powerspectrum.jl")
include("exact_bispectrum.jl")

Ngrid = 64

grid_k = rand(ComplexF32, div(Ngrid, 2) + 1, Ngrid, Ngrid)
dk = 0.01
N = 10
L = 1000
kmax = 0.1

Bk, Nk = bispectrum(grid_k, dk, N, L, kmax)
write_bispectrum(Bk, dk, N, "test.txt")
Pk = power_spectrum(grid_k, dk, N, L)
write_powerspectrum(Pk, dk, "test_pk.txt")