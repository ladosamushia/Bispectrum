using HDF5

include("./compute_bispectrum.jl")

function read_Abacus(dirname)
    first_file = true
    x = nothing
    y = nothing
    z = nothing
    for filename in readdir(dirname)
        if filename[end-4:end] == ".hdf5" 
            println(filename)
            f = h5open(string(dirname, filename))
            pos = read(f["Data"]["pos"])
            xnew = pos[1,:]
            ynew = pos[2,:]
            znew = pos[3,:]
            if first_file
                first_file = false
                x = xnew
                y = ynew
                z = znew
            else
                x = [x; xnew]
                y = [y; ynew]
                z = [z; znew]
            end
        end
    end
    return x, y, z
end

function compute_Abacus_bk(dirname, output_file, dk, N, L)
    x, y, z = read_Abacus(dirname)
    compute_pk_bk(x, y, z, dk, N, L, output_file)
end
