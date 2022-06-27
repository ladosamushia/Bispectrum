using HDF5
using("compute_bispectrum.jl")

function read_molino(filename)
    molinofile = h5open(filename)
    molinopos = read(molinofile["pos"])
    molinovel = read(molinofile["vel"])
    x = molinopos[1,:]
    y = molinopos[2,:]
    z = molinopos[3,:] + molinovel[3,:]/2997.92458
    return x, y, z
end

function compute_molino_bk_pk(input_file, output_file, dk, N, L)
    x, y, z = read_molino(input_file)
    compute_pk_bk(x, y, z, dk, N, L, output_file)
end