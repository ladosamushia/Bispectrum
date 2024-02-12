using DelimitedFiles

include("utilities.jl")
include("bispectrum_utilities.jl")

function bispectrum(gk, N, ind)
    Bk = zeros(N+1, N+1, N+1)
    Nk = zeros(N+1, N+1, N+1)

    for ix1 in 0:N
        iy1max = N^2 - ix1^2 + 1
        l12x = ix1^2
        for iy1 in -ind[iy1max]:ind[iy1max]
            iz1max = iy1max - iy1^2
            l12xy = l12x + iy1^2
            for iz1 in -ind[iz1max]:ind[iz1max]
                if ix1 == 0 && (iy1 > 0 || (iy1 == 0 && iz1 > 0)) continue end
                l12 = l12xy + iz1^2
                if l12 == 0 continue end
                l1 = ind[l12 + 1]
                ix2max = l12 + 1
                ix2min = ind[iy1^2+iz1^2+1]
                for ix2 in -ind[ix2max]:ix2min
                    iy2max = ix2max - ix2^2
                    l22x = ix2^2
                    for iy2 in -ind[iy2max]:ind[iy2max]
                        l22xy = l22x + iy2^2
                        iz2max = iy2max - iy2^2
                        for iz2 in -ind[iz2max]:ind[iz2max]
                            if ix1*ix2 + iy1*iy2 + iz1*iz2 >= 0 continue end
                            l22 = l22xy + iz2^2
                            if l22 < l12/4 || l22 == 0 continue end

                            ix3 = - ix1 - ix2
                            iy3 = - iy1 - iy2
                            iz3 = - iz1 - iz2
                            l32 = ix3^2 + iy3^2 + iz3^2
                            if l32 > l22 || l32 == 0 continue end
            

                            l2 = ind[l22 + 1]
                            l3 = ind[l32 + 1]
                            
                            s = 1
                            if l12 == l22 == l32
                                s = 6
                            elseif l12 == l22 || l22 == l32
                                s = 2
                            end
                            Bk[l1,l2,l3] += real(gk[ix1+N+1,iy1+N+1,iz1+N+1]*gk[ix2+N+1,iy2+N+1,iz2+N+1]*gk[ix3+N+1,iy3+N+1,iz3+N+1])/s
                            Nk[l1,l2,l3] += 1/s     
                        end
                    end
                end
            end
        end    
    end           
    return Bk, Nk
end

function alt_bispectrum(gk, N, ind)
    Bk = zeros(N+1, N+1, N+1)
    Nk = zeros(N+1, N+1, N+1)

    for ix1 in 0:N, 
        iy1m = N^2 - ix1^2
        for iy1 in -ind[iy1m+1]:ind[iy1m+1]
            iz1m = iy1m - iy1^2
            for iz1 in -ind[iz1m+1]:ind[iz1m+1]
                if ix1 == 0 && (iy1 > 0 || (iy1 == 0 && iz1 > 0)) continue end
                l12 = ix1^2 + iz1^2 + iy1^2
                if l12 == 0 continue end
                l1 = ind[l12 + 1]
                l1yz = sqrt(iy1^2+iz1^2+1)
                ix2max = floor(Int, l1yz*0.86602540378-ix1*0.5)
                if ix1 > l1/2
                    ix2min = -l1
                else
                    ix2min = ceil(Int, -ix1*0.5-l1yz*0.86602540378)
                end
                for ix2 in ix2min:ix2max
                    Ryz2 = l12 - ix2^2
                    Ryz = ind[Ryz2 + 1]
                    Tyz = l12/2 + ix1*ix2
                    iRyz1 = iz1^2 + iy1^2
                    D2 = -iz1^2*(Tyz^2 - Ryz2*iRyz1)
                    if D2 < 0 || (iy1^2 + iz1^2) == 0
                        y2min = -Ryz
                        y2max = Ryz
                    else
                        if Tyz + iy1*(-Ryz) <= 0
                            y2min = -Ryz
                        else  
                            y2min = (-Tyz*iy1 - sqrt(D2))/iRyz1
                        end
                        if Tyz + iy1*(Ryz) <= 0
                            y2max = Ryz
                        else
                            y2max = (-Tyz*iy1 + sqrt(D2))/iRyz1
                        end
                    end
                    for iy2 in ceil(Int, y2min):floor(Int, y2max)
                        Tz = Tyz + iy1*iy2
                        Rz2 = l12 - ix2^2 - iy2^2
                        Rz = ind[Rz2+1]
                        if Tz + iz1*(-Rz) <= 0 || iz1 == 0
                            z2min = -Rz
                        else
                            z2min = -Tz/iz1
                        end
                        if Tz + iz1*(Rz) <= 0 || iz1 == 0
                            z2max = Rz
                        else
                            z2max = -Tz/iz1
                        end
                        for iz2 in ceil(Int, z2min):floor(Int, z2max)
                            l22 = ix2^2 + iy2^2 + iz2^2
                            if l22 == 0 continue end
                            ix3 = - ix1 - ix2
                            iy3 = - iy1 - iy2
                            iz3 = - iz1 - iz2
                            l32 = ix3^2 + iy3^2 + iz3^2
                            if l32 == 0 continue end
                            l2 = ind[l22 + 1]
                            l3 = ind[l32 + 1]
                            s = 1
                            if l12 == l22 == l32
                                s = 6
                            elseif l12 == l22 || l22 == l32
                                s = 2
                            end

                            Bk[l1,l2,l3] += real(gk[ix1+N+1,iy1+N+1,iz1+N+1]*gk[ix2+N+1,iy2+N+1,iz2+N+1]*gk[ix3+N+1,iy3+N+1,iz3+N+1])/s
                            Nk[l1,l2,l3] += 1/s     
                        end
                    end
                end
            end
        end    
    end           
    return Bk, Nk
end
