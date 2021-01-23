using FFTW

include("../src/utilities.jl")

"""
    Weight(s, order)

    Assign weight to nearby grid points based on distance.

    # Arguments
    - `s::float`: distance to the grid point in units of grid size.
    - `order::Int`: order of interpolation. can be 1 through 4.
"""
function Weight(s, order)
    if order == 1
        if s < 1.0
            return 1.0
        else
            return 0.0
        end
    end
    if order == 2
        if s < 0.5
            return 1.0 - s
        else
            return 0.0
        end
    end
    if order == 3
        if s < 0.5
            return 0.75 - s^2
        elseif s < 1.5
            return 0.5*(1.5 - s)^2
        else
            return 0.0
        end
    end
    if order == 4
        if s < 1.0
            return 1.0/6.0*(4.0 - 6.0*s^2 + 3.0*s^3)
        elseif s < 2.0
            return 1.0/6.0*(2.0 - s)^3
        else
            return 0.0
        end
    end
    return 0.0
end


"""
    grid_r(Ngrid, x, y, z)

    Make a grid of particle distribution.

    # Arguments
    - `Ngrid::Int`: number of grid cells.
    - `x::float`: x coordinate.
    - `y::float`: y coordinate.
    - `z::float`: z coordinate.
    - `order::Int`: assignment order. can be 1 through 4.

    # Output
    - `grid::array{3}`: number of particles in each grid cell.
"""
function grid_r(Ngrid, x, y, z, order)

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

        for j in -2:2
            ix = index_x + j
            ix = wrap_grid(ix, Ngrid)
            sx = distance_to_grid(x[i], dx, ix, Ngrid)
            Wx = Weight(sx, order)
            for k in -2:2
                iy = index_y + k
                iy = wrap_grid(iy, Ngrid)
                sy = distance_to_grid(y[i], dy, iy, Ngrid)
                Wy = Weight(sy, order)
                for l in -2:2
                    iz = index_z + l
                    iz = wrap_grid(iz, Ngrid)
                    sz = distance_to_grid(z[i], dz, iz, Ngrid)
                    Wz = Weight(sz, order)
                    grid[ix, iy, iz] += Wx*Wy*Wz
                end
            end
        end

    end

    return grid

end

"""
    grid_k(grid_r)

    Fourier transform real space grid.

    # Arguments
    - `grid_r::array{3}`

    # Output
    - `grid_k::array{3,complex}`
"""
function grid_k(grid_r)

    rfft(grid_r)

end