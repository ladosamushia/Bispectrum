using Test

include("../src/powerspectrum.jl")

grid_k = ones(257, 512, 512)
dk = 0.1
Nkbins = 5
L = 1000

Pk = power_spectrum(grid_k, dk, Nkbins, L)
   
for P in Pk
    @test P == 1 * (L/512)^3 / 512^3
end