include("src/bispectrum_utilities.jl")

function exact_bispectrum(grid_k, dk, N, L, kmax)
    Nbins = bispectrum_bins(N)
    Bk = zeros(Nbins)
    Nk = zeros(Nbins)
    Nx, Ny, Nz = size(grid_k)

    kx, ky, kz = Fourier_frequencies(Nz, L)
    k_fundamental = kx[2] - kx[1]
    Nmax = floor(Int, kmax / k_fundamental)

    for ix1 in 1:N, iy1 in -N:N, iz1 in -N:N
        l1 = sqrt(ix1^2 + iy1^2 + iz1^2)
        if l1 > N continue end
        if l1 == 0 continue end
        for ix2 in -N:N, iy2 in -N:N, iz2 in -N:N
            ix3 = - ix1 - ix2
            iy3 = - iy1 - iy2
            iz3 = - iz1 - iz2
            l2 = sqrt(ix2^2 + iy2^2 + iz2^2)
            if l2 == 0 continue end
            if l2 > l1 continue end
            l3 = sqrt(ix3^2 + iy3^2 + iz3^2) 
            if l3 == 0 continue end
            if l3 > l2 continue end
            
            i123 = tri_index(l1, l2, l3, dk / k_fundamental)
            
            i1n, i2n, i3n, j1n, j2n, j3n, k1n, k2n, k3n = get_indeces(ix1, ix2, iz3, ix2, iy2, iz2, ix3, iy3, iz3, Ngrid)
                                    
            Bk_tmp = grid_k[i1n, j1n, k1n]
            if ix2 < 0 Bk_tmp *= conj(grid_k[i2n, j2n, k2n]) else Bk_tmp *= grid_k[i2n, j2n, k2n] end
            if ix3 < 0 Bk_tmp *= conj(grid_k[i3n, j3n, k3n]) else Bk_tmp *= grid_k[i3n, j3n, k3n] end
            
            Nk[i123] += 1 
            Bk[i123] += real(Bk_tmp)           
        end           
    end

    Bk = Bk ./ Nk * (L / Nz)^6 / Nz^3
end