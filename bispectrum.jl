using Threads

include("utilities.jl")
include("bispectrum_utilities.jl")

"""
    bispectrum(grid_k, dk, N, L, kmax)

    Compute bispectrum.

    # Arguments
    - `grid_k::array{3,complex}`: Fourier grid.
    - `dk::float`: k-bin spacing.
    - `N::Int`: number of k bins.
    - `L:float`: Size of the original box.
    - `kmax:float`: Largest k value to look at.

    # Output
    - `Bk::Array{float}`: Binned bispectrum.
"""
function bispectrum(grid_k, dk, N, L, kmax)
    Nbins = bispectrum_bins(N)
    Bk = zeros(nthreads, Nbins)
    Nk = zeros(nthreads, Nbins)

    Nx, Ny, Nz = size(grid_k)

    kx, ky, kz = Fourier_frequencies(Nz, L)

    Nmax = floor(Int, kmax / (kx[2] - kx[1]))

    @threads for i in 1:Nmax
        loop_over_k1k2!(Nmax, i, Nk, Bk, grid_k, threadid, dk)
    end

    Bk = Bk ./ Nk * L^6 / Nz^9
end 