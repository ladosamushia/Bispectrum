using Base.Threads

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
    Bk = zeros(nthreads(), Nbins)
    Nk = zeros(nthreads(), Nbins)

    Nx, Ny, Nz = size(grid_k)

    kx, ky, kz = Fourier_frequencies(Nz, L)
    k_fundamental = kx[2] - kx[1]
    Nmax = floor(Int, kmax / k_fundamental)
    
    @threads for i in 1:Nmax
        loop_over_k1k2!(Nmax, i, Nk, Bk, grid_k, threadid(), dk / k_fundamental)
    end

    Bk = Bk ./ Nk * (L / Nz)^6 / Nz^3
    Bk = sum(Bk, dims=1)

end 