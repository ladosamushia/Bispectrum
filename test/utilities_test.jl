using Test

include("../src/utilities.jl")

kx1, ky1, kz1 = Fourier_frequencies(10, 1)
kx2, ky2, kz2 = Fourier_frequencies(10, 0.1)
kx3, ky3, kz3 = Fourier_frequencies(100,1)

@test isapprox(kx1[2], 2*π)
@test kx2[2] > kx1[2]
@test kx3[end] > kx1[end]