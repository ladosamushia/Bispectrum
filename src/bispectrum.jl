using DelimitedFiles

include("utilities.jl")
include("bispectrum_utilities.jl")

function bispectrum(gk, N, ind)
    Bk = zeros(N, N, N)
    Nk = zeros(N, N, N)

    for ix1 in 0:N, iy1 in -N:N, iz1 in -N:N
        if ix1 == 0 && (iy1 > 0 || (iy1 == 0 && iz1 > 0)) continue end
        l1 = ix1^2 + iy1^2 + iz1^2
        if l1 > N^2 || l1 == 0 continue end
        ind1 = ind[l1]
        for ix2 in -N:N, iy2 in -N:N, iz2 in -N:N
            l2 = ix2^2 + iy2^2 + iz2^2
            if l2 > l1 || l2 == 0 continue end
            ind2 = ind[l2]
            
            ix3 = - ix1 - ix2
            iy3 = - iy1 - iy2
            iz3 = - iz1 - iz2
            l3 = ix3^2 + iy3^2 + iz3^2
            if l3 > l2 || l3 == 0 continue end
            

            ind3 = ind[l3]
            
            s = 1
            if l1 == l2 == l3
                s = 6
            elseif l1 == l2 || l2 == l3
                s = 2
            end
            
            Bk[ind1,ind2,ind3] += real(gk[ix1+N+1,iy1+N+1,iz1+N+1]*gk[ix2+N+1,iy2+N+1,iz2+N+1]*gk[ix3+N+1,iy3+N+1,iz3+N+1])/s
            Nk[ind1,ind2,ind3] += 1/s               
        end           
    end
    return Bk, Nk
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
        inew = div(i-1,Nf) + 1;
        jnew = div(j-1,Nf) + 1;
        knew = div(k-1,Nf) + 1;
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
