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
function bispectrum(grid_k, N, L, klen, kmu)
    # k1 and k2 can not be longer than this
    maxk12 = ceil(Int, sqrt(3)*N)
    # k3 can not be longer than this 
    maxk3 = ceil(Int, 2*sqrt(3)*N)
    Bk = zeros(nthreads(), maxk12, maxk12, maxk3, 3)
    Nk = zeros(nthreads(), maxk12, maxk12, maxk3)

    # kxy coordinates change from -N to N including 0 (2N+1 total)
    # kz coordinate chagnes from 0 to N (N+1) total
    # Julia indexing starts with 0 (0 -> 1, N -> N+1)
    kxyrange = 1:2*N+1
    kzrange = 1:N+1
    @threads for k12xyz in iterators.product(kxyrange,kxyrange,kzrange,kxyrange,kxyrange,kzrange)
        k1x, k1y, k1z, k2x, k2y, k2z = k12xyz
        tid = threadid()
        k3x = k1x + k2x
        k3y = k1y + k2y
        k3z = k1z + k2z
        @inbounds i1 = klen[k1x,k1y,k1z]
        @inbounds i2 = klen[k2x,k2y,k2z]
        @inbounds i3 = klen[k3x,k3y,k3z]
        @inbounds mu = kmu[k1x,k1y,k1z]
        Bk = real(grid_k[k1x,k1y,k1z]*grid_k[k2x,k2y,k2z]*conj(grid_k[k3x,k3y,k3z]))
        @inbounds Nk[tid, i1, i2, i3] += 1 
        @inbounds Bk[tid, i1, i2, i3, 1] += Bk
        @inbounds Bk[tid, i1, i2, i3, 2] += Bk*(1 - mu^2)/2
        @inbounds Bk[tid, i1, i2, i3, 3] += Bk*(3 - 30*mu^2 + 35*mu^4)/8
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
