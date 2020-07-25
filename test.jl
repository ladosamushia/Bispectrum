using FFTW

L = 1000
N = 512

kx = rfftfreq(N, N/L)*2*pi
ky = fftfreq(N, N/L)*2*pi
kz = fftfreq(N, N/L)*2*pi

gr = rand(N,N,N)
gk = rfft(gr)
"""
function tri(a, M)
    Mx = Int(M/2+1)
    for i1 in 1:Mx, j1 in 1:M, k1 in 1:M
        for i2 in 1:Mx, j2 in 1:M, k2 in 1:M
            a[i1,j1,k1]*a[i2,j2,k2]
        end
    end
end

k = zeros(257,512,512)
for i in 1:257, j in 1:512, l in 1:512
    k[i,j,l] = sqrt(kx[i]^2+ky[j]^2+kz[l]^2)
end
"""
function count()
    counter = 0
    for i1 in 1:49 
        println(i1)
        NN = ceil(Int, sqrt(49^2-i1^2))
        for j1 in vcat(1:NN,512-NN+1:512) 
            MM = minimum([1, ceil(Int, sqrt(49^2-i1^2-j1^2))])
            for l1 in vcat(1:MM,512-MM+1:512)
                for i2 in 1:49 
                    NN2 = ceil(Int, sqrt(49^2-i2^2))
                    for j2 in vcat(1:NN2,512-NN2+1:512) 
                        MM2 = minimum([1, ceil(Int, sqrt(49^2-i2^2-j2^2))])
                        for l2 in vcat(1:MM,512-MM+1:512)
                            counter += 1
                        end
                    end
                end
            end
        end
    end
    println(counter)
end
