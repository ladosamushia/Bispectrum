using Test

include("../src/bispectrum_utilities.jl")

@test bispectrum_bins(2) == 4
@test bispectrum_bins(3) == 10

counter = 1
for i in 1:3, j in 1:i, k in 1:j
    @test tri_index(i, j, k, 1) == counter
    global counter += 1
end

#set_k2_min_max(k2min, k2max, l1, i, j, k, i2, j2)

@test wrap_index(0, 257) == 1
@test wrap_index(5, 257) == 6
@test wrap_index(-3, 257) == 255
@test wrap_index(-1, 257) == 257

@test get_indeces(3, -2, 0, 7, 257, -257, -15, 15, -1, 257) == (3, 3, 1, 8, 258, 1, 257 - 14, 16, 257)