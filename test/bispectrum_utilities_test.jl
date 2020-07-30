using Test

include("../src/bispectrum_utilities.jl")

@test bispectrum_bins(2) == 4
@test bispectrum_bins(3) == 10

#set_k2_min_max(k2min, k2max, l1, i, j, k, i2, j2)
