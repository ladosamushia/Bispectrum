using Base.Threads
using DelimitedFiles

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

    Bk = sum(Bk, dims=1) ./ sum(Nk, dims=1) * (L / Nz)^6 / Nz^3
    return Bk, sum(Nk, dims=1)
end 

"""
    write_bispectrum(Bk, dk, N, ofile)

    Write bispectrum and k triplet to a file.

    # Parameters
    - `Bk::Array{Float}`: Bispectrum.
    - `dk::Float`: Bin width.
    - `N::Float`: Number of k bins.
    - `ofile::string`: Output file name.
"""
function write_bispectrum(Bk, dk, N, ofile)
    kbin = collect(range(dk/2, length=N, step=dk))
    output = zeros(bispectrum_bins(N), 4)
    for i in 1:N
        for j in 1:i
            for k in 1:j
                B_index = tri_index(i, j, k, 1)
                output[B_index, :] = [kbin[i], kbin[j], kbin[k], Bk[B_index]]
            end
        end
    end
    writedlm(ofile, output)
end
