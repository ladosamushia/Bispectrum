using Test

include("../src/bispectrum.jl")

grid_k = ones(257, 512, 512)
dk = 0.01
N = 5
L = 1000
kmax = 0.05

Bk = bispectrum(grid_k, dk, N, L, kmax)

for B in Bk
    if isnan(B) == false
        @test B == 1 * (1000 / 512)^6 / 512^3
    end
end
