using JLD2

include("../src/grid.jl")
include("../src/bispectrum.jl")
include("../src/bispectrum_utilities.jl")

function compute_bk(Ngrid, x, y, z, order, outfile, Nf)
    gr = grid_r(Ngrid, x, y, z, order)
    gk = grid_k(gr)
    cutgk = cut_kgrid(Nf, gk)
    
    bk_outfile = string("bk_", outfile)

    ind = compute_indeces(Nf)
    bk, nk = bispectrum(cutgk, Nf, ind)
    
    @save bk_outfile bk nk
end