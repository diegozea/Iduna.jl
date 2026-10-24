using Iduna, JSON, Test

@testset "Path dependency scientific assets" begin
    mktempdir() do root
        utils = Iduna.Utils
        pkg = utils.Pkg
        package = joinpath(root, "BioAlignments")
        cp(pkgdir(utils.BioAlignments), package)
        environment = mkpath(joinpath(root, "environment"))
        project = joinpath(environment, "Project.toml")
        manifest_path = joinpath(environment, "Manifest.toml")
        cp(Base.active_project(), project)
        manifest = pkg.TOML.parsefile(pkg.Types.Context().env.manifest_file)
        selected = pkg.dependencies()
        for entries in values(manifest["deps"]), entry in entries

            if haskey(entry, "path")
                entry["path"] = selected[Base.UUID(entry["uuid"])].source
            end
        end
        dependency = only(manifest["deps"]["BioAlignments"])
        delete!(dependency, "git-tree-sha1")
        dependency["path"] = package
        open(manifest_path, "w") do io
            pkg.TOML.print(io, manifest; sorted = true)
        end
        write(joinpath(root, "seed.sto"), "# STOCKHOLM 1.0\nseq1 AAAA\nother AAAT\n//\n")
        write(joinpath(root, "uniprot.fasta"), ">P20963\nAAAA\n")
        stable_files = [project, manifest_path, joinpath(package, "Project.toml")]
        for (directory, _, names) in walkdir(joinpath(package, "src")), name in names

            endswith(name, ".jl") && push!(stable_files, joinpath(directory, name))
        end
        original_hashes = Dict(path => utils._file_sha256(path) for path in stable_files)
        matrix = joinpath(package, "src", "data", "submat", "BLOSUM62")
        original = read(matrix, String)
        @test occursin("A  4 -1", original)
        probe = joinpath(@__DIR__, "path_dependency_process.jl")
        for phase in ("before", "after")
            if phase == "after"
                chmod(matrix, 0o644)
                write(matrix, replace(original, "A  4 -1" => "A -99 -1"; count = 1))
            end
            # Reload the matrix in a fresh process without reusing compiled dependency data.
            command = `$(Base.julia_cmd()) --startup-file=no --compiled-modules=no --compile=min -O0 --code-coverage=none --project=$environment $probe $root $phase`
            log = joinpath(root, phase * ".log")
            process = open(log, "w") do io
                run(pipeline(ignorestatus(command); stdout = io, stderr = io))
            end
            success(process) || print(read(log, String))
            @test success(process)
            @test Dict(path => utils._file_sha256(path) for path in stable_files) ==
                  original_hashes
        end
        before, after = [JSON.parsefile(joinpath(root, phase * ".json"))
                         for phase in ("before", "after")]
        println("Scientific asset before: ", JSON.json(before))
        println("Scientific asset after: ", JSON.json(after))
        @test before["package_path"] == after["package_path"] == package
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
