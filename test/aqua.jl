using Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(Iduna)
end
