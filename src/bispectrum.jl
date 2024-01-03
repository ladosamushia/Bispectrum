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
    Ngrid = size(grid_k)[3]
    for k1x in 0:N, k1y in -N:N, k1z in -N:N, k2x in -N:N, k2y in -N:N, k2z in -N:N
        k1 = sqrt(k1x^2 + k1y^2 + k1z^2)
        if k1 > N || k1 == 0 continue end
        k2 = sqrt(k2x^2 + k2y^2 + k2z^2)
        if k2 > k1 || k2 == 0 continue end
        k3x = - k1x - k2x
        k3y = - k1y - k2y
        k3z = - k1z - k2z
        k3 = sqrt(k3x^2 + k3y^2 + k3z^2)
        if k3 > k2 || k3 == 0 continue end
        k1y < 0 ? ik1y = Ngrid + k1y + 1 : ik1y = k1y + 1
        k1z < 0 ? ik1z = Ngrid + k1z + 1 : ik1z = k1z + 1
        B_tmp = grid_k[k1x+1,ik1y,ik1z]
        if k2x < 0
            -k2y < 0 ? ik2y = Ngrid + -k2y + 1 : ik2y = -k2y + 1
            -k2z < 0 ? ik2z = Ngrid + -k2z + 1 : ik2z = -k2z + 1
            B_tmp *= conj(grid_k[-k2x+1,ik2y,ik2z])
        else
            k2y < 0 ? ik2y = Ngrid + k2y + 1 : ik2y = k2y + 1
            k2z < 0 ? ik2z = Ngrid + k2z + 1 : ik2z = k2z + 1
            B_tmp *= grid_k[k2x+1,ik2y,ik2z]
        end
        if k3x < 0
            -k3y < 0 ? ik3y = Ngrid + -k3y + 1 : ik3y = -k3y + 1
            -k3z < 0 ? ik3z = Ngrid + -k3z + 1 : ik3z = -k3z + 1
            B_tmp *= conj(grid_k[-k3x+1,ik3y,ik3z])
        else
            k3y < 0 ? ik3y = Ngrid + k3y + 1 : ik3y = k3y + 1
            k3z < 0 ? ik3z = Ngrid + k3z + 1 : ik3z = k3z + 1
            B_tmp *= grid_k[k3x+1,ik3y,ik3z]
        end
        ik1 = floor(Int, k1)
        ik2 = floor(Int, k2)
        ik3 = floor(Int, k3)
        B0[ik1,ik2,ik3] += real(B_tmp)
        B2[ik1,ik2,ik3] += 1
    end 
    return 0
end 

function powerspectrum(grid_k, N, P0, B2, ind)
    for k1x in 0:N, k1y in -N:N, k1z in -N:N
        k1 = k1x^2 + k1y^2 + k1z^2
        if k1 > N^2 || k1 == 0 continue end
        ik1 = ind[k1]
        mu = k1z/k1
        P_tmp = grid_k[k1x+1,k1y+N+1,k1z+N+1]*conj(grid_k[k1x+1,k1y+N+1,k1z+N+1])
       # println(B_tmp)
        P0[ik1] += P_tmp
       # @inbounds B2[ik1,ik2,ik3] += B_tmp*(3*mu^2 - 1)/2*2.5
    end 
    return 0
end 

function compute_bispectrum(B, N, Nf, Ncounts)
    Nb = div(N,Nf)
    Bnew = zeros(Nb, Nb, Nb)
    Nnew = zeros(Int, Nb, Nb, Nb)
    for i in 1:N, j in 1:N, k in 1:N
        inew = div(i-1,Nf) + 1; if inew > Nb continue end
        jnew = div(j-1,Nf) + 1; if jnew > Nb continue end
        knew = div(k-1,Nf) + 1; if knew > Nb continue end
        Bnew[inew,jnew,knew] += B[i,j,k]
        Nnew[inew,jnew,knew] += Ncounts[i,j,k]
    end
    for i in 1:Nb, j in 1:Nb, k in 1:Nb
        if Nnew[i,j,k] != 0
            Bnew[i,j,k] /= Nnew[i,j,k]
        end
    end
    return Bnew
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
            writedlm(f, [kbin[i] kbin[j] kbin[k] Bk[i, j, k]])
        end
    end
    close(f)
end
