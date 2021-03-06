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
    Pk = zeros(nthreads(), Nkbins)
    Nk = zeros(nthreads(), Nkbins)

    Nx, Ny, Nz = size(grid_k)

    kx, ky, kz = Fourier_frequencies(Nz, L)

    @threads for ix in 1:Nx
        loop_over_kykz!(grid_k, Pk, Nk, Nkbins, Ny, Nz, kx, ky, kz, dk, ix, threadid())    
    end

    Pk = sum(Pk, dims=1) ./ sum(Nk, dims=1) * (L / Nz)^3 / Nz^3
end 

function loop_over_kykz!(grid_k, Pk, Nk, Nkbins, Ny, Nz, kx, ky, kz, dk, ix, tid)
    for iy in 1:Ny, iz in 1:Nz
        k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2) 

        if k > 0 && k <= dk*Nkbins
            ik = ceil(Int, k/dk)
            Pk[tid, ik] += abs2(grid_k[ix,iy,iz])
            Nk[tid, ik] += 1
        end
    end   
end

function write_powerspectrum(pk, dk, outfile)
    f = open(outfile, "a")
    k = dk/2
    for p in pk
        writedlm(f, [k p])
        k += dk
    end
    close(f)
end 
