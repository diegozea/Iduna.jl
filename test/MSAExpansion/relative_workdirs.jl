@testset "expansion relative workdirs and relocation" begin
    cd(tmp) do
        fresh = Iduna.MSAExpansion.expand_msa(target, success_seed, "relative_expansion";
            mmseqs_db = success_db, mmseqs_threads = 1)
        @test fresh.status === :ok
        @test fresh.workdir == joinpath(tmp, "relative_expansion")
        @test isfile(fresh.match_stockholm)
        for requested in (fresh.workdir, "relative_expansion")
            cached = Iduna.MSAExpansion.expand_msa(target, success_seed, requested;
                mmseqs_db = success_db, mmseqs_threads = 1)
            @test cached.status === :skipped
            @test cached.match_stockholm == fresh.match_stockholm
        end
        cd(fresh.workdir) do
            cached = Iduna.MSAExpansion.expand_msa(target, success_seed, ".";
                mmseqs_db = success_db, mmseqs_threads = 1)
            @test cached.status === :skipped
            @test cached.match_stockholm == fresh.match_stockholm
        end
        mv(fresh.workdir, "relocated_expansion")
        relocated = Iduna.MSAExpansion.expand_msa(
            target, success_seed, "relocated_expansion";
            mmseqs_db = success_db, mmseqs_threads = 1)
        @test relocated.status === :skipped
        @test isfile(relocated.match_stockholm)
        @test startswith(relocated.match_stockholm, abspath("relocated_expansion"))
        state = Iduna.Utils._read_stage_state(relocated.run_dir)
        @test all(!isabspath(path) for path in keys(state["output_digests"]))
    end
end
