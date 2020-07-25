"""
    bispectrum_bins(N)

    Number of bispectrum bins with N bins in each direction.
"""
function bispectrum_bins(N)

    # k1 >= k2 >= k3 and N bins in k1
    N*(N + 1)*(N + 2)/6

end

function set_k2_min_max(k2min, k2max, l1, i, j, k, i2, j2)
    k2ul = (-l1^2/2 - i*i2 - j*j2)/k
    if k == 0
        if (-l1^2/2 - i*i2 - j*j2) + 1e-10 < 0
            return 1, 0
        end
    elseif k > 0
        if k2ul + 1e-10 > k2max
            # 
        elseif k2ul + 1e-10 < k2min
            return 1, 0
        else
            k2max = ceil(Int, k2ul)
        end
    else
        if k2ul - 1e-10 < k2min
            #
        elseif k2ul - 1e-10 > k2max
            return 1, 0
        else
            k2min = floor(Int, k2ul)
        end
    end
    return k2min, k2max
end

function tri_index(l1, l2, l3, dk)

    k1 = floor(Int, l1/dk)
    k2 = floor(Int, l2/dk)
    k3 = floor(Int, l3/dk)

    ceil(Int, k1*(k1^2 - 1)/6 + k2*(k2 - 1)/2 + k3)

end

function loop_over_k1k2!(N, i, Nk, Bk, grid_k, tid, dk)

    jmin = floor(Int, sqrt(N^2-i^2+1e-10))
    
    for j in -jmin:jmin
    
        kmin = floor(Int, sqrt(N^2-i^2-j^2+1e-10))
    
        for k in -kmin:kmin
    
            l1 = sqrt(i^2+j^2+k^2)
            if l1 == 0 continue end
            i2min = -floor(Int, l1+1e-10)
            i2max = -i2min
    
            for i2 in i2min:i2max
    
                i3 = -i-i2
                j2min = -floor(Int,sqrt(l1^2-i2^2+1e-10))
                j2max = -j2min
    
                for j2 in j2min:j2max
    
                    j3 = -j-j2
                    k2min = -floor(Int, sqrt(l1^2-i2^2-j2^2+1e-10))
                    k2max = -k2min
                    
                    k2min, k2max = set_k2minmax(k2min, k2max, l1, i, j, k, i2, j2)
                   
                    for k2 in k2min:k2max

                        k3 = -k-k2
                        l2 = sqrt(i2^2+j2^2+k2^2)
                        if l2 == 0 continue end
                        l3 = sqrt(i3^2+j3^2+k3^2) 
                        if l3 == 0 continue end
                        if l2 < l3 continue end

                        i123 = tri_index(l1, l2, l3, dk)

                        Nk[tid, i123] += 1 
                        Bk[tid, i123] += 1

                     end

                end

            end

        end

    end

end