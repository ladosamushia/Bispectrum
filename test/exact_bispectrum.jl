include("../src/bispectrum_utilities.jl")

function exact_bispectrum(gk, N)
    Bk = zeros(N, N, N)
    Nk = zeros(Int, N, N, N)

    for ix1 in -N:N, iy1 in -N:N, iz1 in -N:N
        l1 = sqrt(ix1^2 + iy1^2 + iz1^2)
        if l1 > N || l1 == 0 continue end
        for ix2 in -N:N, iy2 in -N:N, iz2 in -N:N
            l2 = sqrt(ix2^2 + iy2^2 + iz2^2)
            if l2 > N || l2 == 0 continue end
            ix3 = - ix1 - ix2
            iy3 = - iy1 - iy2
            iz3 = - iz1 - iz2
            l3 = sqrt(ix3^2 + iy3^2 + iz3^2) 
            if l3 > N || l3 == 0 continue end
            ind1 = floor(Int, l1)
            ind2 = floor(Int, l2)
            ind3 = floor(Int, l3)

            Bk[ind1,ind2,ind3] += real(gk[ix1+N+1,iy1+N+1,iz1+N+1]*gk[ix2+N+1,iy2+N+1,iz2+N+1]*gk[ix3+N+1,iy3+N+1,iz3+N+1])
            Nk[ind1, ind2, ind3] += 1          
        end           
    end
    return Bk, Nk
end