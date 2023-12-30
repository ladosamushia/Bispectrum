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
function bispectrum(grid_k, N, B0, B2, ind)
    for k1x in -N:N, k1y in -N:N, k1z in 0:N, k2x in -N:N, k2y in -N:N, k2z in 0:N
        k1 = k1x^2 + k1y^2 + k1z^2
        if k1 > N^2 || k1 == 0 continue end
        k2 = k2x^2 + k2y^2 + k2z^2
        if k2 > k1 || k2 == 0 continue end
        k3x = k1x + k2x
        k3y = k1y + k2y
        k3z = k1z + k2z
        k3 = k3x^2 + k3y^2 + k3z^2
        if k3 > k2 || k3 == 0 continue end
        ik1 = ind[k1]
        ik2 = ind[k2]
        ik3 = ind[k3]
        mu = k1z/k1
        B_tmp = real(grid_k[k1x+N+1,k1y+N+1,k1z+1]*grid_k[k2x+N+1,k2y+N+1,k2z+1]*conj(grid_k[k3x+N+1,k3y+N+1,k3z+1]))
        @inbounds B0[ik1,ik2,ik3] += B_tmp
        @inbounds B2[ik1,ik2,ik3] += B_tmp*(3*mu^2 - 1)/2*2.5
    end 
    return 0
end 

function compute_bispectrum(B, N)
    Ngrid = size(B)[1]
    for i in 1:Ngrid, j in 1:Ngrid, k in 1:Ngrid
        B[i,j,k] /= N[i,j,k]
    end
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
