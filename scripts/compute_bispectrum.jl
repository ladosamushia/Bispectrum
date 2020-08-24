include("../src/powerspectrum.jl")
include("../src/grid.jl")
include("../src/bispectrum.jl")

"""
    compute_pk_bk(x, y, z, dk, N, L, outfile)

    Compute Bk and Pk given x, y, z, coordinates.

    # Parameters
    - `x, y, z::{Array, Float}`: Cartesian coordinates.
    - `dk::Float`: Bin size.
    - `N::Int`: Number of bins.
    - `L::Float`: Size of the cube.
    - `outfile::string`: Output text file. Will be proceeded by "bk_", "pk_" respectively.
"""
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
