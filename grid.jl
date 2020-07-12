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

function make_shell!(grid_k, grid_k1, k_grid, kmin, kmax)

    Nx, Ny, Nz = size(grid_k)

    for i in 1:256, j in 1:512, k in 1:512
        if k_grid[i,j,k] >= kmin && k_grid[i,j,k] < kmax
            grid_k1[i,j,k] = grid_k[i,j,k]
        else
            grid_k1[i,j,k] = 0
        end
    end

end
        

"""""
Compute Bispectru of a Fourier grid.

"""
function bispectrum(grid_k, dk, Nkbins, L)

    kbinedges = range(0, step=dk, length=Nkbins+1)

    Nx, Ny, Nz = size(grid_k)
    dL = L/Ny
    kx = 2*pi*rfftfreq(Ny, dL) # Ny is not a bug
    ky = 2*pi*fftfreq(Ny, dL)
    kz = 2*pi*fftfreq(Nz, dL)

    k_grid = zeros(size(grid_k))
    for i in 1:256, j in 1:512, k in 1:512
        k_grid[i,j,k] = sqrt(kx[i]^2 + ky[j]^2 + kz[k]^2)
    end

    kmax = dk*Nkbins

    grid_k1 = zeros(ComplexF64,size(grid_k))
    grid_k2 = zeros(ComplexF64,size(grid_k))
    grid_k3 = zeros(ComplexF64,size(grid_k))

    make_shell!(grid_k, grid_k1, k_grid, kbinedges[1], kbinedges[1+1])
    grid_r1 = irfft(grid_k1, Ny)
    make_shell!(grid_k, grid_k2, k_grid, kbinedges[2], kbinedges[2+1])
    grid_r2 = irfft(grid_k2, Ny)
    make_shell!(grid_k, grid_k3, k_grid, kbinedges[3], kbinedges[3+1])
    grid_r3 = irfft(grid_k3, Ny)
    counter = 1
    Threads.@threads for i1 in 1:Nkbins
        println(i1)
        for i2 in 1:i1
            for i3 in i1-i2+1:i2
                sum(grid_r1.*grid_r2.*grid_r3)
                counter += 1
            end
        end
    end
    println(counter)
end
