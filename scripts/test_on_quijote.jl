using PyCall
using DelimitedFiles
using JLD

include("read_quijote_snapshots.jl")
include("../src/powerspectrum.jl")
include("../src/grid.jl")
include("../src/bispectrum.jl")

function quijote_pk_bk(filename, i)
    println(filename)
    bk_filename = string("/projects/QUIJOTE/Lado/Bk_",string(i),".jld")
    pk_filename = string("/projects/QUIJOTE/Lado/Pk_",string(i),".jld")
    if isfile(bk_filename)
        println("file already exists")
        return
    end
    println("reading snapshot")
    xyz = py"read_quijote_snapshot"(filename)
    println("gridding")
    gr = grid_r(512, view(xyz,:,1), view(xyz,:,2), view(xyz,:,3))
    println("Fourier transforming")
    gk = grid_k(gr)
    println("Computing Bk")
    bk = bispectrum(gk, 0.01, 50, 1000, 0.5)
    save(bk_filename, "bk", bk)
    println("Computing Pk")
    pk = power_spectrum(gk, 0.01, 50, 1000)
    save(pk_filename, "pk", pk)
end

for i in 0:999
    filename = string("/projects/QUIJOTE/Snapshots/fiducial/",string(i),"/snapdir_002/snap_002")
    quijote_pk_bk(filename, i)
end
