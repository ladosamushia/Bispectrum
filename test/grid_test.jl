using Test
using Random

include("../src/grid.jl")

x = [0, 4.5, 10]
y = [0, 4.5, 10]
z = [0, 4.5, 10]
gr = grid_r(10, x, y, z, 4)

for i in -1:1, j in -1:1, k in -1:1
    sx = abs(i)
    if sx > 2
        w = 0
    elseif sx < 1
        w = (4 - 6*sx^2 + 3*sx^3)/6
    else
        w = (2 - sx)^3/6
    end
    sy = abs(j)
    if sy > 2
        w *= 0
    elseif sy < 1
        w *= (4 - 6*sy^2 + 3*sy^3)/6
    else
        w *= (2 - sy)^3/6
    end   
    sz = abs(k)
    if sz > 2
        w *= 0
    elseif sz < 1
        w *= (4 - 6*sz^2 + 3*sz^3)/6
    else
        w *= (2 - sz)^3/6
    end                      
    @test isapprox(gr[5 + i, 5 + j, 5 + k], w/size(x)[1]*size(gr)[1]^3 - 1)
end

gr = zeros(10,10,10)
gk = grid_k(gr)
@test all(gk .== 0)

gr = ones(10,10,10)
gk = grid_k(gr)
for i in 2:6
    for j in 1:10
        for k in 1:10
            @test gk[i, j, k] == 0
        end
    end
end
@test isapprox(Real(gk[1, 1, 1]), length(gr); atol=1e-6)

Npart = 100000
x = rand(Npart)*1000
y = rand(Npart)*1000
z = rand(Npart)*1000
gr = grid_r(512, x, y, z, 1)
@test isapprox(sum(gr), 0; atol=1e-6)
gr = grid_r(512, x, y, z, 2)
@test isapprox(sum(gr), 0; atol=1e-6)
gr = grid_r(512, x, y, z, 3)
@test isapprox(sum(gr), 0; atol=1e-6)
gr = grid_r(512, x, y, z, 4)
@test isapprox(sum(gr), 0; atol=1e-6)