using Test

include("../src/bispectrum.jl")
include("../exact_bispectrum.jl")
Ngrid = 64

grid_k = ones(div(Ngrid, 2) + 1, Ngrid, Ngrid)
dk = 0.01
N = 10
L = 1000
kmax = 0.1

Bk = bispectrum(grid_k, dk, N, L, kmax)

for B in Bk
    if isnan(B) == false
        @test B == 1 * (L / Ngrid)^6 / Ngrid^3
    end
end

grid_k = rand(div(Ngrid, 2) + 1, Ngrid, Ngrid)
Bk = bispectrum(grid_k, dk, N, L, kmax)
Bk_exact = exact_bispectrum(grid_k, dk, N, L, kmax)
"""
for i in 1:length(Bk)
    if isnan(Bk[i])
        Bk[i] = 0
    end
    if isnan(Bk_exact[i])
        Bk_exact[i] = 0
    end
end
"""



for i in 1:length(Bk)
    if Bk[i] != Bk_exact[i]
        println(i," ",Bk[i]," ",Bk_exact[i])
    end
    #@test Bk[i] == Bk_exact[i]
end
