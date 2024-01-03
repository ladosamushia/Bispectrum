include("../src/powerspectrum.jl")
include("../src/grid.jl")
include("../src/bispectrum.jl")
include("../src/bispectrum_utilities.jl")

"""
    compute_pk_bk(x, y, z, dk, N, L, outfile)

    Compute Bk and Pk given x, y, z, coordinates.

    # Parameters
    - `x, y, z::{Array, Float}`: Cartesian coordinates.
    - `dk::Float`: Bin size.
    - `Nf::Int`: Number of fundamental modes.
    - `L::Float`: Size of the cube.
    - `outfile::string`: Output text file. Will be proceeded by "bk_", "pk_" respectively.
"""
function compute_pk_bk(x, y, z, dk, Nf, L, outfile, Ncounts)
    gr = grid_r(512, x, y, z, 3)
    gk = grid_k(gr)

    bk_outfile = string("bk_", outfile)

    B0 = zeros(Float64, Nf, Nf, Nf)
    #B2 = zeros(Float64, Nf, Nf, Nf)
    B2 = nothing
    ind = compute_indeces(Nf)
    bk = bispectrum(gk, Nf, B0, B2, ind)
    
    B0ave = compute_bispectrum(B0, Nf, 3)
    
    write_bispectrum(B0ave, dk, 0, 33, bk_outfile)

    pk_outfile = string("pk_", outfile)
    pk = power_spectrum(gk, dk, Nf, L)
    write_powerspectrum(pk, dk, pk_outfile)
end

function compute_pk_only(x, y, z, dk, N, L, outfile)
    gr = grid_r(512, x, y, z, 3)
    gk = grid_k(gr)

    pk_outfile = string("pk_", outfile)
    pk = power_spectrum(gk, dk, N, L)
    write_powerspectrum(pk, dk, pk_outfile)
end

