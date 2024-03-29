"""
    bispectrum_bins(N)

    Number of bispectrum bins with N bins in each direction.
"""
function bispectrum_bins(N)

    # k1 >= k2 >= k3 and N bins in k1
    Int(N*(N + 1)*(N + 2)/6)

end

"""
    set_k2_min_max(k2min, k2max, l1, i, j, k, i2, j2)

    Find minimum and maximum k2_z values that satisfy k3 < k2 < k1

    # Arguments
    - `k2min::float`: minimum value required by k2 < k1
    - `k2max::float`: maximum value required by k2 < k1
    - `l1::float`: |k1|
    - `i,j,k::float`: k1x, k1y, k1z
    - `i2,j2::float`: k2x, k2y

    # Output
    - `k2min,k2max::float`: New values of k2min/k2max, now accounting for k3 < k2 as well
"""
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

"""
    tri_index(l1, l2, l3, dk)

    1D index to arrange B(k1,k2,k3) measurements.

    # Parameters
    - `l1,l2,l3::float`: lengths of k1, k2, k3
    - `dk::float1: bin size
"""
function tri_index(l1, l2, l3, dk, kmin)
    k1 = ceil(Int, (l1 - kmin)/dk)
    k2 = ceil(Int, (l2 - kmin)/dk)
    k3 = ceil(Int, (l3 - kmin)/dk)

    ceil(Int, k1*(k1^2 - 1)/6 + k2*(k2 - 1)/2 + k3)
end

function wrap_index(index, N)
    index <= -1 ? Int(N + index + 1) : Int(index + 1)
end

function unwrap_index(index, N)
    index > div(N, 2) ? Int(index - N - 1) : Int(index - 1)
end

""" 
    wrap_index(i1, i2, i3, j1, j2, j3, k1, k2, k3, N)

    Go from k-indeces for k1, k2, k3 to their position in the fftfreq order.
    Invert the vectors if kx < 0.

    # Parameters
    - `i1, i2, ... ::Int`: Indeces
    - `N::Int`: Grid size

    # Output
    - `i1n, i2n, ...::Int`: Same indes but in fftfreq order.
"""
function get_indeces(i1, i2, i3, j1, j2, j3, k1, k2, k3, N)
    i1n = i1 + 1
    i2n = abs(i2) + 1
    i3n = abs(i3) + 1
    if i1 < 0 j1 = - j1; k1 = - k1 end
    if i2 < 0 j2 = - j2; k2 = - k2 end
    if i3 < 0 j3 = - j3; k3 = - k3 end
    j1n = wrap_index(j1, N)
    j2n = wrap_index(j2, N)
    j3n = wrap_index(j3, N)
    k1n = wrap_index(k1, N)
    k2n = wrap_index(k2, N)
    k3n = wrap_index(k3, N)

    return i1n, i2n, i3n, j1n, j2n, j3n, k1n, k2n, k3n
end

"""
    loop_over_k1k2!(Nmax, i, Nk, Bk, grid_k, tid, dk)

    Loop over k1y, k1x, k1z, k2x, k2y, k2z values while satisfying k1>k2>k3 condition.

    # Parameters
    - `Nmax::Int`: Maximum index needed (set by kmax)
    - `i::Int`: Index of k1x
    - `Nk::array{Int,Nkbin+1}`: Number of triangles going to that bin
    - `Bk::array{Int,Nkbin+1}`: Binned bispectrum
    - `grid_k::array{Complex,3}`: delta_k
    - `tid::Int`: threadid
    - `dk::Float`: Bin width in units of fundamental frequency

    Modifies Bk, Nk, arrays.
"""
function loop_over_k1k2!(Nmax, k0, kmax, k_fundamental, i, Nk, Bk, grid_k, tid, dk)
    Ngrid = size(grid_k)[3]
    
    jmin = floor(Int, sqrt(kmax^2 / k_fundamental^2 - i^2 + 1e-6))
    
    for j in -jmin:jmin
        kmin = floor(Int, sqrt(kmax^2 / k_fundamental^2 - i^2 - j^2 + 1e-6))
    
        for k in -kmin:kmin
            l1 = sqrt(i^2 + j^2 + k^2)
            if l1 <= k0 continue end
            i2min = -floor(Int, l1)
            i2max = -i2min
    
            for i2 in i2min:i2max
                i3 = -i - i2
                j2min = - Nmax
                j2min = -floor(Int, sqrt(l1^2 - i2^2 + eps(Float32)))
                j2max = -j2min

                for j2 in j2min:j2max
                    j3 = -j - j2
                    k2min = -Nmax
                    k2min = -floor(Int, sqrt(l1^2 - i2^2 - j2^2 + eps(Float32)))
                    k2max = -k2min
                    
                    k2min, k2max = set_k2_min_max(k2min, k2max, l1, i, j, k, i2, j2)
                   
                    for k2 in k2min:k2max
                        k3 = -k-k2
                        l2 = sqrt(i2^2 + j2^2 + k2^2)
                        if l2 <= k0 continue end
                        l3 = sqrt(i3^2 + j3^2 + k3^2) 
                        if l3 > l2 || l3 <= k0 continue end

                        ind1 = ceil(Int, (l1 - k0)/dk)
                        ind2 = ceil(Int, (l2 - k0)/dk)
                        ind3 = ceil(Int, (l3 - k0)/dk)

                        i1n, i2n, i3n, j1n, j2n, j3n, k1n, k2n, k3n = get_indeces(i, i2, i3, j, j2, j3, k, k2, k3, Ngrid)
                        
                        Bk_tmp = grid_k[i1n, j1n, k1n]
                        if i2 < 0 Bk_tmp *= conj(grid_k[i2n, j2n, k2n]) else Bk_tmp *= grid_k[i2n, j2n, k2n] end
                        if i3 < 0 Bk_tmp *= conj(grid_k[i3n, j3n, k3n]) else Bk_tmp *= grid_k[i3n, j3n, k3n] end

                        mu = k/l1
                        Nk[tid, ind1, ind2, ind3] += 1 
                        Bk[tid, ind1, ind2, ind3, 1] += real(Bk_tmp)
                        Bk[tid, ind1, ind2, ind3, 2] += real(Bk_tmp)*(1 - mu^2)/2
                        Bk[tid, ind1, ind2, ind3, 3] += real(Bk_tmp)*(3 - 30*mu^2 + 35*mu^4)/8

                     end
                end
            end
        end
    end
end 
