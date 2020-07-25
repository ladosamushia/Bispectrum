function wrap_index(i, N)
    return i < 0 ? Int(N + i + 1) : Int(i + 1)
end

function count(A, Nmax, dk, Nk)
    Nbins = some_counts(Nk)
    cc = zeros(Nbins)
    nn = zeros(Nbins)
    counter = 0
    for i1 in 0:Nmax
        # k1 < kmax
        jmax1 = floor(Int, sqrt(Nmax^2 - i1^2))
        for j1 in -jmax1:jmax1
            # k1 < kmax
            kmax1 = floor(Int, sqrt(Nmax^2 - i1^2 - j1^2))
            for k1 in -kmax1:kmax1
                k1length = sqrt(i1^2 + j1^2 + k1^2)
                k1_i = ceil(Int, k1length/dk)
                # k2 < k1, k3 < k2. 
                if i1 != 0
                    i2max = min(ceil(Int,-k1length^2-j1*k1length-k1*k1length)/2/i1, floor(Int,k1length))
                else
                    i2max = floor(Int,k1length)
                end
                i2min = -floor(Int,k1length)
                for i2 in i2min:i2max
                    i3 = i1 + i2
                    if j1 != 0
                        j2max = min(ceil(Int,-k1length^2-i1*i2-k1*k1length/2/j1), floor(Int,sqrt(k1length^2 - i2^2)))
                    else
                        j2max = floor(Int,sqrt(k1length^2 - i2^2))
                    end
                    j2min = -floor(Int,sqrt(k1length^2 - i2^2))
                    for j2 in j2min:j2max
                        if j2^2 > k1length^2 - i2^2 continue end
                        j3 = j1 + j2
                        if k1 != 0
                            k2max = min(ceil(Int,-k1length^2-i1*i2-j1*j2)/2/k1, floor(Int,sqrt(k1length^2 - i2^2 - j2^2)))
                        else
                            k2max = floor(Int,sqrt(k1length^2 - i2^2 - j2^2))
                        end
                        k2min = -floor(Int,sqrt(k1length^2 - i2^2 - j2^2))
                        for k2 in k2min:k2max
                            k2length = sqrt(i2^2 + j2^2 + k2^2)
                            k2_i = ceil(Int, k2length/dk)
                            k3 = k1 + k2
                            k3length = sqrt(i3^2 + j3^2 + k3^2)
                            k3_i = ceil(Int, k3length/dk)
                            k_123 = tri_index(k1_i, k2_i, k3_i)
                            iy1 = wrap_index(j1, 512)
                            iz1 = wrap_index(k1, 512)
                            iy2 = wrap_index(j2, 512)
                            iz2 = wrap_index(k2, 512)
                            iy3 = wrap_index(j3, 512)
                            iz3 = wrap_index(k3, 512)
                            if k_123 == 0
                                continue
                            end
                            #cc[k_123] += real(A[i1+1,iy1,iz1]*A[i2+1,iy2,iz2]*conj(A[i3+1,iy3,iz3]))
                            nn[k_123] += 1
                            counter += 1
                        end
                    end
                end
            end
        end
    end
    println(counter)
    return cc, nn
end

function some_test(N)
    n = 0
    for i in -N:N
        for j in -N:N
            for k in -N:N
#            if i == 0 && j == 0 && k == 0 continue end
                for i2 in -N:N
                    i3 = -i-i2
                    for j2 in -N:N
                        j3 = -j-j2
                        for k2 in -N:N
#                            if i2 == 0 && j2 == 0 && k2 == 0 continue end
                            k3 = -k-k2
#                            if i3 == 0 && j3 == 0 && k3 == 0 continue end
                            l1 = sqrt(i^2+j^2+k^2)
                            l2 = sqrt(i2^2+j2^2+k2^2)
                            if l1 < l2 continue end
                            l3 = sqrt(i3^2+j3^2+k3^2)
                            if l2 < l3 continue end
                            if l1 <= N && l2 <= N && l3 <= N
                                n += 1
                                #println(i," ",j," ",k," ",i2," ",j2," ",k2," ")
                            end
                        end
                    end
                end
            end
        end
    end
    println(n)
end

function some_function(N,i,n)
        jmin = floor(Int, sqrt(N^2-i^2+1e-10))
        for j in -jmin:jmin
            kmin = floor(Int, sqrt(N^2-i^2-j^2+1e-10))
            for k in -kmin:kmin
                l1 = sqrt(i^2+j^2+k^2)
                if l1 == 0 continue end
                i2min = -floor(Int, l1+1e-10)
                i2max = -i2min
               """ 
                if i == 0
                    #
                elseif i > 0
                    i2ul = (-l1^2/2 + abs(j)*l1 + abs(k)*l1)/i
                    if i2ul + 1e-10 > i2max
                       # 
                    elseif i2ul + 1e-10 < i2min
                        continue
                    else
                        i2max = ceil(Int, i2ul)
                    end
                else
                    i2ul = (-l1^2/2 + abs(j)*l1 + abs(k)*l1)/i
                    if i2ul - 1e-10 < i2min
                        #
                    elseif i2ul - 1e-10 > i2max
                        continue
                    else
                        i2min = floor(Int, i2ul)
                    end
                end               
               """ 
                for i2 in i2min:i2max
                    i3 = -i-i2
                    j2min = -floor(Int,sqrt(l1^2-i2^2+1e-10))
                    j2max = -j2min
                   """ 
                    if j == 0
                        #
                    elseif j > 0
                        j2ul = (-l1^2/2 - i*i2 + abs(k)*sqrt(l1^2-i2^2+1e-10))/j
                        if j2ul + 1e-10 > j2max
                           # 
                        elseif j2ul + 1e-10 < j2min
                            continue
                        else
                            j2max = ceil(Int, j2ul)
                        end
                    else
                        j2ul = (-l1^2/2 - i*i2 + abs(k)*sqrt(l1^2-i2^2+1e-10))/j
                        if j2ul - 1e-10 < j2min
                            #
                        elseif j2ul - 1e-10 > j2max
                            continue
                        else
                            j2min = floor(Int, j2ul)
                        end
                    end
                   """ 
                    for j2 in j2min:j2max
                        j3 = -j-j2
                        k2min = -floor(Int, sqrt(l1^2-i2^2-j2^2+1e-10))
                        k2max = -k2min
                        
                        k2ul = (-l1^2/2 - i*i2 - j*j2)/k
                        if k == 0
                            if (-l1^2/2 - i*i2 - j*j2) + 1e-10 < 0
                                continue
                            end
                        elseif k > 0
                            if k2ul + 1e-10 > k2max
                               # 
                            elseif k2ul + 1e-10 < k2min
                                continue
                            else
                                k2max = ceil(Int, k2ul)
                            end
                        else
                            if k2ul - 1e-10 < k2min
                                #
                            elseif k2ul - 1e-10 > k2max
                                continue
                            else
                                k2min = floor(Int, k2ul)
                            end
                        end
                       
                        for k2 in k2min:k2max
                            k3 = -k-k2
                            l2 = sqrt(i2^2+j2^2+k2^2)
                            l3 = sqrt(i3^2+j3^2+k3^2)
                            if l2 < l3 continue end
                            if l2 == 0 continue end
                            if l3 == 0 continue end
                            n[i+1] += 1 
                            #println(i," ",j," ",k," ",i2," ",j2," ",k2," ")
                        end
                    end
                end
            end
        end
end

function some_test_2(N)
    n = zeros(N+1)
    Threads.@threads for i in 0:N
        some_function(N,i,n)
    end
    println(sum(n))
end

function all_pairs(N)
    n = 0
    for i in 0:N
        jmin = floor(Int, sqrt(N^2-i^2))
        for j in -jmin:jmin
            kmin = floor(Int, sqrt(N^2-i^2-j^2))
            for k in -kmin:kmin
                if i == 0 && j == 0 && k == 0 continue end
                l1 = sqrt(i^2+j^2+k^2)
                i2min = -floor(Int, l1)
                i2max = -i2min
                for i2 in i2min:i2max
                    i3 = -i-i2
                    j2min = -floor(Int,sqrt(l1^2-i2^2))
                    j2max = -j2min
                    for j2 in j2min:j2max
                        j3 = -j-j2
                        k2min = -floor(Int, sqrt(l1^2-i2^2-j2^2))
                        k2max = -k2min
                        for k2 in k2min:k2max
                            k3 = -k-k2
                            l2 = sqrt(i2^2+j2^2+k2^2)
                            l3 = sqrt(i3^2+j3^2+k3^2)
                            #if l2 < l3 continue end
                            if i2 == 0 && j2 == 0 && k2 == 0 continue end
                            if i3 == 0 && j3 == 0 && k3 == 0 continue end
                            n += 1
                            println(i," ",j," ",k," ",i2," ",j2," ",k2," ")
                        end
                    end
                end
            end
        end
    end
    println(n)
end

function some_counts(N)

    n = 0
    for i in 1:N
        for j in 1:i
            for k in 1:j
                n += 1
            end
        end
    end
    return n
end

function tri_index(k1, k2, k3)
    Int(k1*(k1^2 - 1)/6 + k2*(k2-1)/2 + k3)
end
