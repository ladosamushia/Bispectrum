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

    rfftw(grid_r)

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
    dL: float
        grid cell size in real space

Output:
    Pk: Array 
        Power spectrum

"""

function power_spectrum(grid_k, dk, Nk, dL)

    Pk = zeros(Nk)
    Nk = zeros(Nk)

    Nx, Ny, Nz = size(grid_k)
    kx = 2*pi*rfftfreq(Ny, dL) # Ny is not a bug
    ky = 2*pi*fftfreq(Ny, dL)
    kz = 2*pi*fftfreq(Nz, dL)

    for ix in 1:Nx, iy in 1:Ny, iz in 1:Nz

        k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2) 
        ik = ceil(Int, k/dk)
        if ik <= Nk
            Pk[ik] += grid_k[ix,iy,iz]
	    Nk[ik] += 1
        end
    end

    Pk = Pk./Nk

end 
