using Base.Threads
using DelimitedFiles

include("utilities.jl")
include("bispectrum_utilities.jl")

"""
    bispectrum(grid_k, dk, N, L, kmin, kmax)

    Compute bispectrum.

    # Arguments
    - `grid_k::array{3,complex}`: Fourier grid.
    - `dk::float`: k-bin spacing.
    - `N::Int`: number of k bins.
    - `L:float`: Size of the original box.
    - `kmin:float`: Smallest k value to look at.
    - `kmax:float`: Largest k value to look at.

    # Output
    - `Bk::Array{float}`: Binned bispectrum.
"""
function bispectrum(grid_k, dk, N, L, kmin, kmax)
    #Nbins = bispectrum_bins(N)
    Bk = zeros(nthreads(), N, N, N, 3)
    Nk = zeros(nthreads(), N, N, N)
    Nx, Ny, Nz = size(grid_k)

    kx, ky, kz = Fourier_frequencies(Nz, L)
    k_fundamental = kx[2] - kx[1]
    Nmax = floor(Int, kmax / k_fundamental)
    println(Nmax, " ", k_fundamental, " ", kmin, " ", kmax, " ", dk)
    @threads for i in 0:Nmax
        loop_over_k1k2!(Nmax, kmin / k_fundamental, kmax, k_fundamental, i, Nk, Bk, grid_k, threadid(), dk / k_fundamental)
    end
   
    Bk = sum(Bk, dims=1) ./ sum(Nk, dims=1) * (L/Nz)^6 /Nz^3
    return Bk
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
function write_bispectrum(Bk, dk, kmin, N, ofile)
    kbin = collect(range(kmin + dk/2, length=N, step=dk))
    f = open(ofile, "w")
    for i in 1:N, j in i:N, k in j:N
        if k < j + i
            writedlm(f, [kbin[i] kbin[j] kbin[k] Bk[1, k, j, i, 1] Bk[1, k, j, i, 2] Bk[1, k, j, i, 3]])
        end
    end
    close(f)
end
