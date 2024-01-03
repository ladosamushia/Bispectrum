function cut_kgrid(N, grid_k)
    Ngrid = size(grid_k)[3]
    cut_grid_k = zeros(Complex{Float64},N+1,2*N+1,2*N+1)
    for i in 1:N+1, j in 1:2*N+1, k in 1:2*N+1
        j <= N ? jold = Ngrid - N + j : jold = j - N
        k <= N ? kold = Ngrid - N + k : kold = k - N
        cut_grid_k[i,j,k] = grid_k[i,jold,kold]
    end
    return cut_grid_k
end

function compute_indeces(N)
    Nsq = N^2
    ind = zeros(Int8, Nsq+1)
    for i in 1:Nsq+1
        ind[i] = isqrt(i)
    end
    return ind
end

function compute_Nk(N, ind)
    Nk = zeros(Int, N, N, N)
    cc = 0
    for k1x in 0:N, k1y in -N:N, k1z in -N:N, k2x in -N:N, k2y in -N:N, k2z in -N:N
        k1 = k1x^2 + k1y^2 + k1z^2
        if k1 > N^2 || k1 == 0 continue end
        k2 = k2x^2 + k2y^2 + k2z^2
        if k2 > k1 || k2 == 0 continue end
        k3x = - k1x - k2x
        k3y = - k1y - k2y
        k3z = - k1z - k2z
        k3 = k3x^2 + k3y^2 + k3z^2
        if k3 > k2 || k3 == 0 continue end
        ik1 = ind[k1]
        ik2 = ind[k2]
        ik3 = ind[k3]
        Nk[ik1,ik2,ik3] += 1
        cc += 1
    end
    println(cc)
    return Nk
end

function compute_Nk_alt(N, ind)
    Nk = zeros(Int, N, N, N)
    cc = 0
    for k1x in 0:N, k1y in -N:N, k1z in -N:N, k2x in -N:N, k2y in -N:N, k2z in -N:N
        k1 = isqrt(k1x^2 + k1y^2 + k1z^2)
        if k1 > N || k1 == 0 continue end
        k2 = isqrt(k2x^2 + k2y^2 + k2z^2)
        if k2 > k1 || k2 == 0 continue end
        k3x = - k1x - k2x
        k3y = - k1y - k2y
        k3z = - k1z - k2z   
        k3 = isqrt(k3x^2 + k3y^2 + k3z^2) 
        if k3 > k2 || k3 == 0 continue end
        Nk[k1,k2,k3] += 1
        cc += 1
    end
    println(cc)
    return Nk
end