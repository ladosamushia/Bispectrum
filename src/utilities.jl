using FFTW

"""
    Fourier_frequencies(N, L)

    Return Fourier modes.

    # Arguments: 
    - `N::Int`: Size of the grid.
    - `L::Float`: Size of the cube (in Mpc or Mpc/h).

    # Output:
    - `kx:array`
    - `ky:array`
    - `kz:array`
"""

function Fourier_frequencies(N, L)
    kx = 2*pi*rfftfreq(N, N/L)
    ky = 2*pi*fftfreq(N, N/L)
    kz = 2*pi*fftfreq(N, N/L)

    return kx, ky, kz
end

"""
    wrap_grid(i, Ngrid)

    Make sure index does not go outside 1 - Ngrid (periodic cube).

    # Arguments:
    - `i::Int`: index.
    - `Ngrid::Int`: size of the grid.

    # Output:
    - `i::Int`: New index wrapped if necessary. Always between 1 - Ngrid.
"""
function wrap_grid(i, Ngrid)
    if i < 1
        i += Ngrid
    elseif i > Ngrid
        i -= Ngrid
    end
    return i
end

"""
    wrap_L(x, L)

    Make sure x is between 0 and L

    # Arguments:
    - `x::float`: coordinate
    - `L::float`: grid size

    # Output:
    - `xnew::float`: New coordinate brought inside 0 - L interval
"""
function wrap_L(x, L)
    xnew = x
    if x < 0
        xnew = L + x
    elseif x > L
        xnew = x - L
    end
    return xnew
end

"""
    distance_to_grid(x, dL, index)

    Return distance to grid point indexed by index. Wraps around periodic cubes properly.

    # Arguments:
    - `x:float`: coordinate of a point
    - `dL::float`: grid size
    - `index::Int`: index of the grid point (dL*index distance away from the origin)
    - `Ngrid::Int`: grid size
"""
function distance_to_grid(x, dL, index, Ngrid)
    s = abs(x - (dL*index - dL/2))/dL
    if s > Ngrid/2
        s = Ngrid - s
    end
    return s
end
