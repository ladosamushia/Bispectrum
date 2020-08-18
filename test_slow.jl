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

function the_loop(x1, k, kmax, dk, kf, Bk, Nk, gk, tid)
    y1max = floor(Int, sqrt(kmax^2 - k[x1]^2) / kf)
    for y1 in flatten((1:y1max, 512 - y1max + 2:512))
        z1max = floor(Int, sqrt(kmax^2 - k[x1]^2 - k[y1]^2) / kf) + 1
        for z1 in flatten((1:z1max, 512-z1max+2:512))
            k1 = sqrt(k[x1]^2 + k[y1]^2 + k[z1]^2)
            if k1 > kmax || k1 == 0 continue end
            x2max = floor(Int, k1 / kf)
            for x2 in flatten((1:x2max, 512 - x2max + 2:512))
                y2max = floor(Int, sqrt(k1^2 - k[x2]^2 + 1) / kf)
                for y2 in flatten((1:y2max, 512 - y2max + 2:512))
                    z2max = floor(Int, sqrt(k1^2 - k[x2]^2 - k[y2]^2 + 1) / kf)
                    for z2 in flatten((1:z2max, 512 - z2max + 1:512))
                        k2 = sqrt(k[x2]^2 + k[y2]^2 + k[z2]^2)
                        if k2 > k1 || k2 == 0 continue end
                        k3 = sqrt((k[x1] + k[x2])^2 + (k[y1] + k[y2])^2 + (k[z1] + k[z2])^2)
                        if k3 > k1 || k3 == 0 continue end
                        ix2 = x2 > 257 ? 512 - x2 + 2 : x2
                        x3 = get_k3(x1, x2)
                        ix3 = x3 > 257 ? 512 - x3 + 2 : x3
                        y3 = get_k3(y1, y2)
                        z3 = get_k3(z1, z2)
                        i1 = ceil(Int, k1 / dk)
                        i2 = ceil(Int, k2 / dk)
                        i3 = ceil(Int, k3 / dk)
                        Bk_tmp = gk[x1, y1, z1]
                        if x2 > 257 Bk_tmp *= conj(gk[ix2, y2, z2]) else Bk_tmp *= gk[ix2, y2, z2] end
                        if x3 > 257 Bk_tmp *= conj(gk[ix3, y3, z3]) else Bk_tmp *= gk[ix3, y3, z3] end
                        Bk[tid, i1, i2, i3] += real(Bk_tmp)
                        Nk[tid, i1, i2, i3] += 1
                    end
                end
            end
        end
    end
end

function exact_but_slow(gk, kmax, Nbin)
    Bk = zeros(nthreads(), Nbin, Nbin, Nbin)
    dk = kmax / Nbin
    Nk = zeros(nthreads(), Nbin, Nbin, Nbin)
    k = fftfreq(512, 512 / 1000) * 2 * π
    kf = k[2] - k[1]
    x1max = floor(Int, kmax / kf)
    @threads for x1 = 1:x1max
        the_loop(x1, k, kmax, dk, kf, Bk, Nk, gk, threadid())
    end
    return sum(Bk, dims = 1) ./ sum(Nk, dims = 1)
end
