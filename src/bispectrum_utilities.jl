function cut_kgrid(N, grid_k)
    Nx = size(grid_k)[1]
    Nyz = size(grid_k)[3]
    kx = Int.(rfftfreq(Nyz, Nyz))
    ky = Int.(fftfreq(Nyz, Nyz))
    kz = Int.(fftfreq(Nyz, Nyz))
    
    cut_grid_k = zeros(Complex{Float64},2*N+1,2*N+1,2*N+1)
    for i in 1:Nx, j in 1:Nyz, k in 1:Nyz
        if kx[i] <= N && abs(ky[j]) <= N && abs(kz[k]) <= N
            cut_grid_k[kx[i]+N+1,ky[j]+N+1,kz[k]+N+1] = grid_k[i,j,k]
            cut_grid_k[-kx[i]+N+1,-ky[j]+N+1,-kz[k]+N+1] = conj(grid_k[i,j,k])
        end
    end
    
    return cut_grid_k
end

function simmetrize_bispectrum(B)
    N = size(B)[3]
    Bsym = zeros(N, N, N)
    for i in 1:N, j in 1:N, k in 1:N
        Bsym[i,j,k] = (B[i,j,k] + B[i,k,j] + B[j,i,k] + B[j,k,i] + B[k,i,j] + B[k,j,i])/6
    end
    return Bsym
end
                
function compute_indeces(N)
    Nsq = N^2
    ind = zeros(Int16, Nsq+1)
    for i in 0:Nsq
        ind[i+1] = isqrt(i)
    end
    return ind
end

function reduce_bispectrum(B, N, Nf, Ncounts)
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

function write_bispectrum(Bk, dk, kmin, N, ofile)
    kbin = collect(range(kmin + dk/2, length=N, step=dk))
    f = open(ofile, "w")
    for i in 1:N, j in 1:N, k in 1:N
        if B[i,j,k] != 0
            writedlm(f, [kbin[i] kbin[j] kbin[k] Bk[i, j, k]])
        end
    end
    close(f)
end

function compute_Nk(N)
    Nk = zeros(N+1, N+1, N+1)
    for ix1 in 0:N, iy1 in -N:N, iz1 in -N:N
        if ix1 == 0 && (iy1 > 0 || (iy1 == 0 && iz1 > 0)) continue end
        l1sq = ix1^2 + iy1^2 + iz1^2
        if l1sq > N^2 || l1sq == 0 continue end
        for ix2 in -l1:l1, iy2 in -l1:l1, iz2 in -l1:l1
            l2sq = ix2^2 + iy2^2 + iz2^2
            if l2sq > l1sq || l2sq == 0 continue end
            ix3 = - ix1 - ix2
            iy3 = - iy1 - iy2
            iz3 = - iz1 - iz2
            l3sq = ix3^2 + iy3^2 + iz3^2
            if l3sq > l2sq || l3sq == 0 continue end
            l1 = ind[l1sq + 1]
            l2 = ind[l2sq + 1]
            l3 = ind[l3sq + 1]
                            
            s = 1
            if l1sq == l2sq == l3sq
                s = 6
            elseif l1sq == l2sq || l2sq == l3sq
                s = 2
            end
                            
            Nk[l1,l2,l3] += 1/s     
        end    
    end           
    return Nk
end