using Iduna
using Test
using TestItemRunner

@testset "Iduna.jl" begin
    include("aqua.jl")
end

@run_package_tests
