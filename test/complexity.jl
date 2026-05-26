using CodeComplexity

@testset "Code quality (CodeComplexity.jl)" begin
    @test isempty(check_measure(CyclomaticComplexity(), Iduna))
    @test isempty(check_measure(CognitiveComplexity(), Iduna))
end
