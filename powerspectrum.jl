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
    Pk = zeros(Nkbins)
    Nk = zeros(Nkbins)

    Nx, Ny, Nz = size(grid_k)

    kx, ky, kz = Fourier_frequencies(Nz, L)

    for ix in 1:Nx, iy in 1:Ny, iz in 1:Nz
        k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2) 

        if k > 0 && k <= dk*Nkbins
            ik = ceil(Int, k/dk)
            Pk[ik] += abs2(grid_k[ix,iy,iz])
            Nk[ik] += 1
        end       
    end

    Pk = Pk ./ Nk * (L / Nz)^3 / Nz^3
end 