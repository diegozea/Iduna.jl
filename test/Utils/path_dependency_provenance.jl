using Iduna, JSON, Test

include("path_dependency_fixture.jl")

@testset "Path dependency scientific assets" begin
    mktempdir() do root
        fixture = path_dependency_fixture(root)
        before = run_path_dependency_process(fixture, "before")
        edit_path_dependency_matrix(fixture)
        after = run_path_dependency_process(fixture, "after")
        println("Scientific asset before: ", JSON.json(before))
        println("Scientific asset after: ", JSON.json(after))
        @test before["package_path"] == after["package_path"] == fixture.package
        @test before["blosum62_aa"] == 4
        @test after["blosum62_aa"] == -99
        @test before["source_sha256"] != after["source_sha256"]
        @test before["implementation_hash"] != after["implementation_hash"]
        @test before["action"] == "run"
        @test after["action"] == "rebuild"
        @test before["fresh"] != after["fresh"]
        @test after["cached"] == after["fresh"]
        @test after["cached"]["status"] == "warn"
        @test !after["cached"]["identical"]
        @test after["cached"]["insertions"] == after["cached"]["deletions"] == 4
    end
end
