include("../src/powerspectrum.jl")
include("../src/grid.jl")
include("../src/bispectrum.jl")

function compute_pk_bk(x, y, z, dk, N, L, outfile)
    gr = grid_r(512, x, y, z)
    gk = grid_k(gr)

    kmax = dk*N
    bk_outfile = string("bk_", outfile)
    bk = bispectrum(gk, dk, N, L, kmax)
    write_bispectrum(bk, dk, N, bk_outfile)

    pk_outfile = string("pk_", outfile)
    pk = power_spectrum(gk, dk, N, L)
    write_powerspectrum(pk, dk, pk_outfile)
end
