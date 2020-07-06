using FFTW

"""
Make a grid of particle distribution

Parameters:
    Ngrid: integer
        number of grid cells
    x: array
    y: array
    z: array

Output:
    grid: 3D array
        number of particles on a grid

"""
function grid_r(Ngrid, x, y, z)

    xmax = maximum(x)
    xmin = minimum(x)
    ymax = maximum(y)
    ymin = minimum(y)
    zmax = maximum(z)
    zmin = minimum(z)

    Lx = xmax - xmin
    Ly = ymax - ymin
    Lz = zmax - zmin

    dx = Lx/Ngrid
    dy = Ly/Ngrid
    dz = Lz/Ngrid

    grid = zeros(Ngrid,Ngrid,Ngrid)

    for i in eachindex(x, y, z)

        index_x = ceil(Int, x[i]/dx)
        index_y = ceil(Int, y[i]/dy)
        index_z = ceil(Int, z[i]/dz)

        # Only fill in immediate neighbours
        for step_x = -2:1, step_y = -2:1, step_z = -2:1

            i_x = index_x + step_x
            i_y = index_y + step_y
            i_z = index_z + step_z
            # Enforce periodicity. Stay inside the box.
            if i_x < 1
                i_x = Ngrid
            end
            if i_x > Ngrid
                i_x = 1
            end
            if i_y < 1
                i_y = Ngrid
            end
            if i_y > Ngrid
                i_y = 1
            end
            if i_z < 1
                i_z = Ngrid
            end
            if i_z > Ngrid
                i_z = 1
            end

            dist_x = x[i] - i_x*dx
            dist_y = y[i] - i_y*dy
            dist_z = z[i] - i_z*dz

            s = sqrt(dist_x^2 + dist_y^2 + dist_z^2)/dx

            if s > 2
                weight = 0
            elseif s < 1
                weight = (4 - 6*s^2 + 3*s^3)/6
            else
                weight = (2 - s)^3/6
            end

            grid[i_x, i_y, i_z] += weight

        end

    end

    return grid

end

"""
Fourier transform real space grid.

Parameters:
    grid_r: 3D array

Output:
    grid_k: 3D complex array
"""
function grid_k(grid_r)

    rfft(grid_r)

end

"""
Compute power spectrum.

Parameters:
    grid_k: 3D complex array
        Fourier grid
    dk: float
        k-bin spacing
    Nk: int
        number of k bins
    L: float
        Size of the original box

Output:
    Pk: Array 
        Power spectrum

"""
function power_spectrum(grid_k, dk, Nkbins, L)

    Pk = zeros(Nkbins)
    Nk = zeros(Nkbins)

    Nx, Ny, Nz = size(grid_k)
    dL = L/Ny
    kx = 2*pi*rfftfreq(Ny, dL) # Ny is not a bug
    ky = 2*pi*fftfreq(Ny, dL)
    kz = 2*pi*fftfreq(Nz, dL)

    for ix in 1:Nx, iy in 1:Ny, iz in 1:Nz

        k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2) 
        if k > 0 && k <= dk*Nkbins
            ik = ceil(Int, k/dk)
            Pk[ik] += abs2(grid_k[ix,iy,iz])
            Nk[ik] += 1
        end

    end
    println(Nk)
    Pk = Pk./Nk*L^3/Ny^6

end 

"""
Pre-compute k1, k2, k3 below k_max and satisfying trinagular condition.
This needs to only be done once and will speed up the bispectrum computations.

Parameters: 
    Ngrid: Integer
        The grid size
    dL: Float
        Cell size in real space
    kmax: Float
        Maximum k to go to

Output:
    k_list: Array{9, ?}
        
Output is index of k1, index of k2, k1, k2, k3
"""
function k_pairs(Ngrid, dL, kmax)

    Nx = Int(Ngrid/2)
    Ny = Int(Ngrid)
    Nz = Int(Ngrid)

    kx = 2*pi*rfftfreq(Ngrid, dL) 
    ky = 2*pi*fftfreq(Ngrid, dL)
    kz = 2*pi*fftfreq(Ngrid, dL)
    println(kx)
    k_list = zeros(9, 1)
    for ix1 in 1:Nx, iy1 in 1:Ny, iz1 in 1:Nz
        k1 = sqrt(kx[ix1]^2 + ky[iy1]^2 + kz[iz1]^2)
        if k1 > kmax
            continue
        end
        for ix2 in 1:Nx, iy2 in 1:Ny, iz2 in 1:Nz
            k2 = sqrt(kx[ix2]^2 + ky[iy2]^2 + kz[iz2]^2)
            if k2 > k1
                continue
            end
            k3x = kx[ix1] - kx[ix2];
            k3y = ky[iy1] - ky[iy2];
            k3z = kz[iz1] - kz[iz2];
            k3 = sqrt(k3x^2 + k3y^2 + k3z^2)
            if k3 > k2 || k1 > k2 + k3
                continue
            end 
            k_list = hcat(k_list, [ix1, iy1, iz1, ix2, iy2, iz2, k1, k2, k3])
        end
    end
    return k_list
end

"""
Compute Bispectru of a Fourier grid.

"""
function bispectrum(grid_k, dk, Nkbins, L)

    Nx, Ny, Nz = size(grid_k)
    dL = L/Ny
    kx = 2*pi*rfftfreq(Ny, dL) # Ny is not a bug
    ky = 2*pi*fftfreq(Ny, dL)
    kz = 2*pi*fftfreq(Nz, dL)

    kmax = dk*Nkbins

    for ix1 in 1:Nx, iy1 in 1:Ny, iz1 in 1:Nz
        k1 = sqrt(kx[ix1]^2 + ky[iy1]^2 + kz[iz1]^2)
        if k1 < kmax
            continue
        end
        for ix2 in 1:Nx, iy2 in 1:Ny, iz2 in 1:Nz
            k2 = sqrt(kx[ix2]^2 + ky[iy2]^2 + kz[iz2]^2)
            if k2 < kmax || k2 > k1
                continue
            end
            ix3 = ix1 + ix2
            iy3 = iy1 + iy2
            iz3 = iz1 + iz2
            grid_k[ix1,iy1,iz1]*grid_k[ix2,iy2,iz2]*grid_k[ix3,iy3,iz3]
        end
    end

end
