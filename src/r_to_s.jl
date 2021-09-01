include("utilities.jl")

"""
   r_to_s(z, vz, a, H, L)

   Change the real-space position to the redshift-space.

    #Arguments
    - `z::Array{float}`: real-space coordinate
    - `vz::Array{float}`: redshift-space coordinate
    - `a::float`: the scale factor
    - `H::float`: Hubble parameter at that redshift in units of h
    - `L::float`: Size of the simulation
"""

function r_to_s(z, vz, a, H, L)
   znew = copy(z) 
   for i in 1:length(z)
       zz = z[i] + vz[i]/a/H
       znew[i] = wrap_L(zz, L)
    end
    return znew
end 
