using Test
using FFTW

include("../src/bispectrum_utilities.jl")

function compute_indeces_test(N)
    Nsq = N^2
    ind = zeros(Int8, Nsq+1)
    for i in 1:Nsq+1
        ind[i] = isqrt(i)
    end
    return ind
end

ind1 = compute_indeces(100)
ind2 = compute_indeces_test(100)

for i in 1:100
    @test ind1[i] == ind2[i]
end

gr = rand(512, 512, 512)
gkr = rfft(gr)
gkc = fft(gr)

gkcut = cut_kgrid(10, gkr)