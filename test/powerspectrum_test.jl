using Test

include("../src/powerspectrum.jl")

grid_k = ones(257, 512, 512)
dk = 0.1
Nkbins = 5
L = 1000

Pk = power_spectrum(grid_k, dk, Nkbins, L)
   
for i in size(Pk)[2]
    @test Pk[1,i,1] == 1 * (L/512)^3 / 512^3
    @test -1e-6 < Pk[1,i,2] < 1e-6
    @test -1e-6 < Pk[1,i,3] < 1e-6
end