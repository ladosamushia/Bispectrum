using Test

include("../src/utilities.jl")

kx1, ky1, kz1 = Fourier_frequencies(10, 1)
kx2, ky2, kz2 = Fourier_frequencies(10, 0.1)
kx3, ky3, kz3 = Fourier_frequencies(100,1)

@test isapprox(kx1[2], 2*π)
@test kx2[2] > kx1[2]
@test kx3[end] > kx1[end]

@test wrap_grid(-1, 256) == 255
@test wrap_grid(257, 256) == 1
@test wrap_grid(120, 256) == 120

@test wrap_L(5, 10) == 5
@test wrap_L(-1.9, 100) == 98.1
@test wrap_L(1409, 1000) == 409
