include("../src/r_to_s.jl")

z = [0.3, 56, 998.1, 502]
v = [-3.0, 3, 5.0, 4]
a = 0.4
H = 0.5
L = 1000

znew = r_to_s(z, v, a, H, L)

@test isapprox(znew[1],985.3)
@test isapprox(znew[2],71.0)
@test isapprox(znew[3],23.1)
@test isapprox(znew[4],522.0)
