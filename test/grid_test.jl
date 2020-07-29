using Test

include("../src/grid.jl")

x = [0, 5, 10]
y = [0, 5, 10]
z = [0, 5, 10]
gr = grid_r(10, x, y, z)

for i in -3:3
    for j in -3:3
        for k in -3:3
            l = sqrt(i^2 + j^2 + k^2)
            if l > 2
                w = 0
            elseif l < 1
                w = (4 - 6*l^2 + 3*l^3)/6
            else
                w = (2 - l)^3/6
            end
            @test isapprox(gr[5 + i, 5 + j, 5 + k], w)
        end
    end
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
@test gk[1, 1, 1] == 1000