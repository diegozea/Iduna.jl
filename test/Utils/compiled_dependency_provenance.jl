@testset "Compiled path dependency scientific assets" begin
    mktempdir() do root
        fixture = path_dependency_fixture(root)
        run_path_dependency_process(fixture, "warmup"; compiled = true, warmup = true)
        edit_path_dependency_matrix(fixture)
        before = run_path_dependency_process(fixture, "before"; compiled = true)
        same = run_path_dependency_process(fixture, "same"; compiled = true)
        state_path = joinpath(root, "shared", "validation", "pid_10.00", "stage_state.json")
        state = JSON.parsefile(state_path)
        delete!(state["identity"]["implementation"], "protein_alignment_sha256")
        state["identity"]["implementation"]["contract"] = 1
        state["identity_hash"] = Iduna.Utils._identity_hash(state["identity"])
        Iduna.Utils.write_json(state_path, state)
        obsolete = run_path_dependency_process(fixture, "obsolete"; compiled = true)
        # Retain only this fixture's old cache; shared installations are never changed.
        compiled = joinpath(
            fixture.depot, "compiled", "v$(VERSION.major).$(VERSION.minor)",
            "BioAlignments")
        @test isdir(compiled)
        mv(compiled, joinpath(root, "retained-compiled-BioAlignments"))
        after = run_path_dependency_process(fixture, "after"; compiled = true)
        repeated = run_path_dependency_process(fixture, "repeated"; compiled = true)
        println("Compiled scientific asset before: ", JSON.json(before))
        println("Compiled scientific asset after: ", JSON.json(after))
        @test before["package_path"] == after["package_path"] == fixture.package
        @test before["blosum62_aa"] == before["model_aa"] == 4
        @test after["blosum62_aa"] == after["model_aa"] == -99
        @test before["source_sha256"] == after["source_sha256"]
        @test before["model_sha256"] != after["model_sha256"]
        @test before["implementation_hash"] != after["implementation_hash"]
        @test same["implementation_hash"] == before["implementation_hash"]
        @test repeated["implementation_hash"] == after["implementation_hash"]
        @test before["action"] == "run"
        @test same["action"] == repeated["action"] == "reuse"
        @test obsolete["action"] == "rebuild"
        @test obsolete["implementation_hash"] == before["implementation_hash"]
        @test after["action"] == "rebuild"
        @test before["fresh"] != after["fresh"]
        for result in (before, same, obsolete, after, repeated)
            @test result["cached"] == result["fresh"]
        end
        @test after["cached"]["status"] == "warn"
        @test !after["cached"]["identical"]
        @test after["cached"]["insertions"] == after["cached"]["deletions"] == 4
    end
end
