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
    ind = zeros(Int8, Nsq+1)
    for i in 1:Nsq+1
        ind[i] = isqrt(i)
    end
    return ind
end