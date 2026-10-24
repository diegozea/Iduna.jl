@testset "load_result preserves legacy custom-PID artifacts" begin
    archive = joinpath(@__DIR__, "fixtures", "pre59-custom-pid")
    contents(root) = Dict(relpath(joinpath(dir, file), root) => read(joinpath(dir, file))
    for (dir, _, files) in walkdir(root) for file in files)
    original = contents(archive)
    @test length(original) == 11
    @test !isfile(joinpath(archive, "thoraxe_msa", "candidate_summary.csv"))

    for layout in (:original, :moved, :coexisting_precise)
        @testset "$layout" begin
            mktempdir() do tmp
                workdir = joinpath(tmp, "original")
                cp(archive, workdir)
                if layout === :moved
                    moved = joinpath(tmp, "moved")
                    mv(workdir, moved)
                    @test !ispath(workdir)
                    workdir = moved
                elseif layout === :coexisting_precise
                    for dir in ("candidates", "runs")
                        root = joinpath(workdir, "thoraxe_msa", dir)
                        cp(joinpath(root, "pid_30.00"), joinpath(root, "pid_30.001"))
                    end
                end
                before = contents(workdir)
                loaded = Iduna.load_result(workdir)
                candidate = joinpath("thoraxe_msa", "candidates", "pid_30.00")
                sequence = joinpath(candidate, "sequences", "candidate_sequences_full.fasta")
                species = joinpath(candidate, "species", "candidate_species_full.txt")
                thoraxe = joinpath("thoraxe_msa", "runs", "pid_30.00", "full", "thoraxe")
                seed = only(loaded.thoraxe_msa.seeds)
                @test loaded.workdir == workdir
                @test seed.pid == 30.001
                @test seed.stockholm_path == joinpath(candidate, "candidate_msa_full.sto")
                @test loaded.thoraxe_msa.sequence_fastas == [sequence]
                @test loaded.thoraxe_msa.species_files == [species]
                @test loaded.thoraxe_msa.thoraxe_dirs == [thoraxe]
                @test all(isfile(joinpath(workdir, path))
                for path in loaded.thoraxe_msa.sequence_fastas)
                @test all(isfile(joinpath(workdir, path))
                for path in loaded.thoraxe_msa.species_files)
                @test all(isdir(joinpath(workdir, path))
                for path in loaded.thoraxe_msa.thoraxe_dirs)
                @test read(joinpath(workdir, sequence), String) == ">ENSG1\nACDE\n"
                @test read(joinpath(workdir, species), String) == "homo_sapiens\n"
                @test isfile(joinpath(workdir, thoraxe, "path_table.csv"))
                @test Iduna.ResultsValidation.stringsequence(Iduna.load_seed_msa(loaded), 1) ==
                      "ACDE"
                @test contents(workdir) == before
            end
        end
    end
    @test contents(archive) == original
end
