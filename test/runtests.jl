using Test, SafeTestsets

begin
    @safetestset "Utilities" begin include("utilities_test.jl") end
    @safetestset "Grid" begin include("grid_test.jl") end
    @safetestset "Bispectrum utilities" begin include("bispectrum_utilities_test.jl") end
    @safetestset "Power Spectrum" begin include("powerspectrum_test.jl") end
    @safetestset "Bispectrum" begin include("bispectrum_test.jl") end
    @safetestset "r_to_s" begin include("r_to_s_test.jl") end
end
