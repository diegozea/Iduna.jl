using Iduna
using Test
using Aqua

@testset "Iduna.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Iduna)
    end
    # Write your tests here.
end
