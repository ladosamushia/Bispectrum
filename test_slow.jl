using Base.Iterators
using Base.Threads
using FFTW

function get_k3(i1, i2)
    i1n = i1 > 257 ? i1 - 512 - 1 : i1 - 1
    i2n = i2 > 257 ? i2 - 512 - 1 : i2 - 1
    i3n = i1n + i2n
    if i3n < 0
        i3n = i3n + 512 + 1
    else
        i3n = i3n + 1
    end
end

function the_loop(x1, Nmax, k, kmax, dk, Bk, Nk, gk, tid)
    for y1 in flatten((1:Nmax, 512 - Nmax + 1:512)), z1 in flatten((1:Nmax, 512 - Nmax + 1:512))
        k1 = sqrt(k[x1]^2 + k[y1]^2 + k[z1]^2)
        if k1 > kmax || k1 == 0 continue end
        for x2 in 1:Nmax, y2 in flatten((1:Nmax, 512 - Nmax + 1:512)), z2 in flatten((1:Nmax, 512 - Nmax + 1:512))
            k2 = sqrt(k[x2]^2 + k[y2]^2 + k[z2]^2)
            if k2 > k1 || k2 == 0 continue end
            k3 = sqrt((k[x1] + k[x2])^2 + (k[y1] + k[y2])^2 + (k[z1] + k[z2])^2)
            if k3 > k1 || k3 == 0 continue end
            x3 = x1 + x2 - 1
            y3 = get_k3(y1, y2)
            z3 = get_k3(z1, z2)
            i1 = ceil(Int, k1 / dk)
            i2 = ceil(Int, k2 / dk)
            i3 = ceil(Int, k3 / dk)
            Bk[tid, i1, i2, i3] += real(gk[x1, y1, z1] * gk[x2, y2, z2] * conj(gk[x3, y3, z3]))
            Nk[tid, i1, i2, i3] += 1
        end
    end
end

function exact_but_slow(gk, kmax, Nbin)
    Bk = zeros(nthreads(), Nbin, Nbin, Nbin)
    dk = kmax / Nbin
    Nk = zeros(nthreads(), Nbin, Nbin, Nbin)
    k = fftfreq(512, 512/1000)*2*π
    kf = k[2] - k[1]
    Nmax = ceil(Int, kmax / kf)
    @threads for x1 in 1:Nmax 
        the_loop(x1, Nmax, k, kmax, dk, Bk, Nk, gk, threadid()) 
    end
    return Bk, Nk
end
