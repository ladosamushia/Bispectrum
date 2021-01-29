using PyCall
using DelimitedFiles
using JLD
using Random

include("read_quijote_snapshots.jl")
include("compute_bispectrum.jl")

function quijote_pk_bk(input_file, output_file, dk, N, L)
    println(input_file)
    println("reading snapshot")
    xyz = py"read_quijote_snapshot"(input_file)
    println("gridding")
    compute_bk_pk(view(xyz,:,1), view(xyz,:,2), view(xyz,:,3), dk, N, L, output_file)
end

function quijote_pk_bk_subsample(input_file, output_file, dk, N, L)
    println(input_file)
    println("reading snapshot")
    xyz = py"read_quijote_snapshot"(input_file)
    # subsample. I want to achieve n ~ 10e-4
    i_subsample = randsubseq(MersenneTwister(), 1:size(xyz)[1], 100000.0/size(xyz)[1])
    println("gridding")
    compute_bk_pk(view(xyz,i_subsample,1), view(xyz,i_subsample,2), view(xyz,i_subsample,3), dk, N, L, output_file)
end
