<<<<<<< HEAD
using FITSIO

include("./compute_bispectrum.jl")

# Only read files that start with fnamebeg

function read_Abacus(dirname, fnamebeg)
    first_file = true
    x = nothing
    y = nothing
    z = nothing
    for filename in readdir(dirname)
        if filename[end-4:end] == ".fits" & filename[1:length(fnamebeg)] ==
fnamebeg
            println(filename)
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

function compute_Abacus_bk(dirname, fnamebeg, output_file, dk, N, L)
    x, y, z = read_Abacus(dirname, fnamebeg)
    compute_pk_bk(x, y, z, dk, N, L, output_file)
end
=======
using FITSIO

include("./compute_bispectrum.jl")

# Only read files that start with fnamebeg

function read_Abacus(dirname, fnamebeg)
    first_file = true
    x = nothing
    y = nothing
    z = nothing
    for filename in readdir(dirname)
        if filename[end-4:end] == ".fits" && filename[1:length(fnamebeg)] ==
fnamebeg
            println(filename)
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

function compute_Abacus_bk(dirname, fnamebeg, output_file, dk, N, L)
    x, y, z = read_Abacus(dirname, fnamebeg)
    compute_pk_bk(x, y, z, dk, N, L, output_file)
end
>>>>>>> multipoles
