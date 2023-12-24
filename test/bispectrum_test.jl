using Test

include("../src/bispectrum.jl")
include("exact_bispectrum.jl")
Ngrid = 64

gk = ones(div(Ngrid, 2) + 1, Ngrid, Ngrid)
dk = 0.01
N = 10
L = 1000
kmax = 0.1

Bk, Nk = bispectrum(gk, dk, N, L, 0, kmax)

for B in Bk
    if isnan(B) == false
        @test B == 1 * (L / Ngrid)^6 / Ngrid^3
    end
end

gk = rand(ComplexF32, div(Ngrid, 2) + 1, Ngrid, Ngrid)
Bk, Nk = bispectrum(gk, dk, N, L, 0, kmax)
Bk_exact, Nk_exact = exact_bispectrum(gk, dk, N, L, 0, kmax)

for i in eachindex(Bk)
    if isnan(Bk[i]) == false
        @test isapprox(Bk[i], Bk_exact[i])
    end
end

dk = 0.01
N = 10
L = 1000
kmin = 0.01
kmax = kmin + dk*N

Bk, Nk = bispectrum(gk, dk, N, L, kmin, kmax)
Bk_exact, Nk_exact = exact_bispectrum(gk, dk, N, L, kmin, kmax)

for i in eachindex(Bk)
    if isnan(Bk[i]) == false
        @test isapprox(Bk[i], Bk_exact[i])
    end
end