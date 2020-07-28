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