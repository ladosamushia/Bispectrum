using DelimitedFiles
using Base.Threads

include("utilities.jl")

"""
    power_spectrum(grid_k, dk, Nkbins, L)

    Compute power spectrum.

    # Arguments
    - `grid_k::array{3,complex}`: Fourier grid.
    - `dk::float`: k-bin spacing.
    - `Nk::Int`: number of k bins.
    - `L:float`: Size of the original box.

    # Output
    - `Pk::Array{float}`: Binned power spectrum.
"""
function power_spectrum(grid_k, dk, Nkbins, L)
    Pk = zeros(nthreads(), Nkbins, 3)
    Nk = zeros(nthreads(), Nkbins)

    Nx, Ny, Nz = size(grid_k)

    kx, ky, kz = Fourier_frequencies(Nz, L)

    @threads for ix in 1:Nx
        loop_over_kykz!(grid_k, Pk, Nk, Nkbins, Ny, Nz, kx, ky, kz, dk, ix, threadid())    
    end

    Pk = sum(Pk, dims=1) ./ sum(Nk, dims=1) * L^3 / Nz^6
end 

function loop_over_kykz!(grid_k, Pk, Nk, Nkbins, Ny, Nz, kx, ky, kz, dk, ix, tid)
    for iy in 1:Ny, iz in 1:Nz
        k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2) 
        mu = kz[iz] / k
        if k > 0 && k <= dk*Nkbins
            ik = ceil(Int, k/dk)
            d_square = abs2(grid_k[ix,iy,iz])
            Pk[tid, ik, 1] += d_square
            Pk[tid, ik, 2] += d_square*(1 - mu^2)/2
            Pk[tid, ik, 3] += d_square*(3 - 30*mu^2 + 35*mu^4)/8
            Nk[tid, ik] += 1
        end
    end   
end

function write_powerspectrum(pk, dk, outfile)
    f = open(outfile, "a")
    k = dk/2
    for i in 1:size(pk)[2]
        writedlm(f, [k pk[1,i,1] pk[1,i,2] pk[1,i,3]])
        k += dk
    end
    close(f)
end 
