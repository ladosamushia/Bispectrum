using FITSIO

include("./compute_bispectrum.jl")

function read_Abacus(dirname)
    first_file = true
    x = nothing
    y = nothing
    z = nothing
    for filename in readdir(dirname)
        if (filename[end-4:end] == ".fits")
            f = FITS(string(dirname, filename))
            xnew = read(f[2], "x")
            ynew = read(f[2], "y")
            znew = read(f[2], "z")
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

function compute_Abacus_bk(dirname, Ngrid, order, outfile, Nf)
    x, y, z = read_Abacus(dirname)
    compute_bk(Ngrid, x, y, z, order, outfile, Nf)
end
