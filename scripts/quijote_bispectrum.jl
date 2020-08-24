using PyCall
using DelimitedFiles
using JLD

include("read_quijote_snapshots.jl")
include("compute_bispectrum.jl")

function quijote_pk_bk(input_file, output_file, dk, N, L)
    println(input_file)
    println("reading snapshot")
    xyz = py"read_quijote_snapshot"(input_file)
    println("gridding")
    compute_bk_pk(view(xyz,:,1), view(xyz,:,2), view(xyz,:,3), dk, N, L, output_file)
end
