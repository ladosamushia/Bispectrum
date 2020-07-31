using Test

include("../src/bispectrum.jl")

grid_k = ones(257, 512, 512)
dk = 0.01
N = 10
L = 1000
kmax = 0.1
"""
Bk = bispectrum(grid_k, dk, N, L, kmax)

for B in Bk
    if isnan(B) == false
        @test B == 1 * (1000 / 512)^6 / 512^3
    end
end
"""
grid_k = zeros(ComplexF32, 257, 512, 512)
grid_k[1, 4, 2] = 1 + 2im
grid_k[3, 512 + 1 - 2, 3] = 2 - 3im
grid_k[2, 512 + 1 - 2, 512 + 1 - 5] = 1 - 1im

Bk = bispectrum(grid_k, dk, N, L, kmax)
for i in 1:length(Bk)
    if isnan(Bk[i])
        Bk[i] = 0
    end
end
@test sum(Bk) == real((1 + 2im)*(2 + 3im)*(1 - 1im))
