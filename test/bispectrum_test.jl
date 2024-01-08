using Test
using FFTW

include("exact_bispectrum.jl")
include("../src/bispectrum_utilities.jl")
include("../src/bispectrum.jl")

N = 512
gr = rand(Float64, N, N, N)
gk = rfft(gr)

gkcut = cut_kgrid(10, gk)

ind = compute_indeces(100);

Bk, Nk = exact_bispectrum(gkcut, 10)
alt_Bk, alt_Nk = bispectrum(gkcut, 10, ind)

Bk = simmetrize_bispectrum(Bk)

for i in 1:10, j in 1:10, k in 1:10
    if alt_Nk[i,j,k] != 0
        @test isapprox(Bk[i,j,k]./Nk[i,j,k], alt_Bk[i,j,k]./alt_Nk[i,j,k])
    end
end
