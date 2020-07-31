using Test

include("../src/bispectrum.jl")

Ngrid = 64

grid_k = ones(div(Ngrid, 2) + 1, Ngrid, Ngrid)
dk = 0.01
N = 10
L = 1000
kmax = 0.1

Bk = bispectrum(grid_k, dk, N, L, kmax)

for B in Bk
    if isnan(B) == false
        @test B == 1 * (1000 / Ngrid)^6 / Ngrid^3
    end
end

grid_k = zeros(ComplexF32, div(Ngrid, 2) + 1, Ngrid, Ngrid)

grid_k[1 + 1, 4 + 1, 2 + 1] = 1 + 2im
grid_k[3 + 1, -2 + Ngrid + 1, 3 + 1] = 2 - 3im
grid_k[2 + 1, -2 + Ngrid + 1, -5 + Ngrid + 1] = 1 - 1im

Bk = bispectrum(grid_k, dk, N, L, kmax)
for i in 1:length(Bk)
    if isnan(Bk[i])
        Bk[i] = 0
    end
end
@test sum(Bk) == real((1 + 2im)*(2 + 3im)*(1 - 1im)) * (1000 / Ngrid)^6 / Ngrid^3
