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
    dL = L/N

    kx = 2*pi*rfftfreq(N, dL)
    ky = 2*pi*fftfreq(N, dL)
    kz = 2*pi*fftfreq(N, dL)

    return kx, ky, kz
end