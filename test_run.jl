include("src/bispectrum_utilities.jl")
include("src/powerspectrum.jl")
include("src/bispectrum.jl")
include("test.jl")

n, Nk = number_of_triangles(5, 1, 5)
println(n)
println(Nk)
println(sum(Nk))

println(bispectrum_bins(5))
Ngrid = 64

grid_k = zeros(ComplexF32, div(Ngrid, 2) + 1, Ngrid, Ngrid)
grid_k[1 + 1, 4 + 1, 2 + 1] = 1 + 2im
grid_k[3 + 1, -2 + Ngrid + 1, 3 + 1] = 2 - 3im
grid_k[2 + 1, -2 + Ngrid + 1, -5 + Ngrid + 1] = 1 - 1im

L = 1000

bk, nk = bispectrum(grid_k, Ngrid / L * 2 * pi, 5, L, 5 * Ngrid / L * 2 * pi)
println(bk)
println(sum(nk))


Nk = zeros(1, 1000)
Bk = zeros(1, 1000)
loop_over_k1k2!(6, 2, Nk, Bk, grid_k, 1, 1)