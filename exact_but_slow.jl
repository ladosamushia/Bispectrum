using FFTW
using Base.Threads
using Base.Iterators

include("src/bispectrum_utilities.jl")
include("src/utilities.jl")

Ngrid = 512
L = 1000
dk = 0.01
kmax = 0.1

function bispectrum_exact(grid_k, dk, L, kmax)
    Ngrid = size(grid_k)[3]

    kx, ky, kz = Fourier_frequencies(Ngrid, L)
    kx = copy(ky)
    Nmax = floor(Int, kmax / (kx[2] - kx[1]))

    N = ceil(Int, kmax / dk)
    Nbin = bispectrum_bins(N)
    Bk = zeros(nthreads(), Nbin)
    Nk = zeros(nthreads(), Nbin)

    for i1 in 1:Nmax
        @threads for i in 1:Nmax
            loop_over_k1k2_exact(kx, ky, kz, Nmax, Ngrid, i, threadid(), grid_k, Nk, Bk)
        end
    end
    return sum(Bk, dims=1) ./ sum(Nk, dims=1) * (L / Nz)^6 / Nz^3
end

function loop_over_k1k2_exact(kx, ky, kz, Nmax, Ngrid, i1, tid, grid_k, Nk, Bk)
    kx1 = kx[i1]
    for j1 in flatten((1:Nmax, Ngrid - Nmax + 1:Ngrid))
        ky1 = ky[j1]
        for k1 in flatten((1:Nmax, Ngrid - Nmax + 1:Ngrid))
            kz1 = kz[k1]
            l1 = sqrt(kx1^2 + ky1^2 + kz1^2)
            if l1 > kmax || l1 == 0 continue end
            for i2 in flatten((1:Nmax, Ngrid - Nmax + 1:Ngrid))
                kx2 = kx[i2]
                kx3 = - kx1 - kx2
                for j2 in flatten((1:Nmax, Ngrid - Nmax + 1:Ngrid))
                    ky2 = ky[j2]
                    ky3 = - ky1 - ky2
                    for k2 in flatten((1:Nmax, Ngrid - Nmax + 1:Ngrid))
                        kz2 = kz[k2]
                        kz3 = - kz1 - kz2
                        l2 = sqrt(kx2^2 + ky2^2 + kz2^2)
                        if l2 > l1 || l2 == 0 continue end
                        l3 = sqrt(kx3^2 + ky3^2 + kz3^2)
                        if l3 > l2 || l3 == 0 continue end
                        i123 = tri_index(l1, l2, l3, dk)
                        Nk[tid, i123] += 1

                        i2n, i3, j3, k3 = k3_indeces(i1, i2, j1, j2, k1, k2, Ngrid)
                        
                        Bk_tmp = grid_k[i1, j1, k1]
                        if kx2 < 0 Bk_tmp *= conj(grid_k[i2n, j2, k2]) else Bk_tmp *= grid_k[i2n, j2, k2] end
                        if kx3 < 0 Bk_tmp *= conj(grid_k[i3, j3, k3]) else Bk_tmp *= grid_k[i3, j3, k3] end
                        Bk[tid, i123] += real(Bk_tmp)
                    end
                end
            end
        end
    end
end

function k3_indeces(i1, i2, j1, j2, k1, k2, Ngrid)
    i2n = unwrap_index(i2, Ngrid)
    j2n = unwrap_index(j2, Ngrid)
    k2n = unwrap_index(k2, Ngrid)
    j1n = unwrap_index(j1, Ngrid)
    k1n = unwrap_index(k1, Ngrid)

    i3 = - i1 - i2n
    j3 = - j1n - j2n
    k3 = - k1n - k2n

    i3n = abs(i3) + 1
    j3n = wrap_index(j3, Ngrid)
    k3n = wrap_index(k3, Ngrid)

    i2n = abs(i2n) + 1

    return i2n, i3n, j3n, k3n
end
