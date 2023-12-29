using DelimitedFiles

include("utilities.jl")
include("bispectrum_utilities.jl")

"""
    bispectrum(grid_k, N, B0, B2, ind)

    Compute bispectrum.

    # Arguments
    - `grid_k::array{3,complex}`: Fourier grid.
    - `N::Int`: Only go up to k = N*k_fundamental
    - `B0::array{3,Float}': Array to store Bispectrum monopole
    - `B2::array{3,Float}': Array to store Bispectrum quadrupole
    - `ind::array{Int8}': Array to store the integer square roots

    B0 and B2 are not normalized (/Nk). This will happen later. There is no need
to recompute Nk as it is a geometric factor and does not change.
    Some entries in B0 and B2 will be zero. The ones that don't satisfy the
triangular equality.
    The isosceles triangles will be double counted.
    The equilateral triangles will be six-fold counted.
"""
function bispectrum(grid_k, N, B, ind)
    for k1x in -N:N, k1y in -N:N, k1z in 0:N, k2x in -N:N, k2y in -N:N, k2z in 0:N
        k1 = k1x^2 + k1y^2 + k1z^2+1
        if k1 >= N^2 continue end
        k2 = k2x^2 + k2y^2 + k2z^2+1
        if k2 >= k1 continue end
        k3x = k1x + k2x
        k3y = k1y + k2y
        k3z = k1z + k2z
        k3 = k3x^2 + k3y^2 + k3z^2+1
        if k3 >= k2 continue end
            ik1 = ind[k1]
            ik2 = ind[k2]
            ik3 = ind[k3]
        @inbounds B[ik1,ik2,ik3] += real(grid_k[k1x+N+1,k1y+N+1,k1z+1]*grid_k[k2x+N+1,k2y+N+1,k2z+1]*conj(grid_k[k3x+N+1,k3y+N+1,k3z+1]))
    end 
    return 0
end 

"""
    write_bispectrum(Bk, dk, N, ofile)

    Write bispectrum and k triplet to a file.

    # Parameters
    - `Bk::Array{Float}`: Bispectrum.
    - `dk::Float`: Bin width.
    - `N::Float`: Number of k bins.
    - `ofile::string`: Output file name.
"""
function write_bispectrum(Bk, dk, kmin, N, ofile)
    kbin = collect(range(kmin + dk/2, length=N, step=dk))
    f = open(ofile, "w")
    for i in 1:N, j in i:N, k in j:N
        if k < j + i
            writedlm(f, [kbin[i] kbin[j] kbin[k] Bk[1, k, j, i, 1] Bk[1, k, j, i, 2] Bk[1, k, j, i, 3]])
        end
    end
    close(f)
end
