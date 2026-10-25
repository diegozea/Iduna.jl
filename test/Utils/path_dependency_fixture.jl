function path_dependency_fixture(root)
    utils = Iduna.Utils
    pkg = utils.Pkg
    package = joinpath(root, "BioAlignments")
    cp(pkgdir(utils.BioAlignments), package)
    environment = mkpath(joinpath(root, "environment"))
    project = joinpath(environment, "Project.toml")
    project_data = pkg.TOML.parsefile(Base.active_project())
    project_data["deps"]["BioAlignments"] = string(Base.PkgId(utils.BioAlignments).uuid)
    open(project, "w") do io
        pkg.TOML.print(io, project_data; sorted = true)
    end
    manifest_path = joinpath(environment, "Manifest.toml")
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
    stable_hashes = Dict(path => utils._file_sha256(path) for path in stable_files)
    matrix = joinpath(package, "src", "data", "submat", "BLOSUM62")
    original = read(matrix, String)
    @test occursin("A  4 -1", original)
    depot = mkpath(joinpath(root, "depot"))
    return (; root, environment, package, matrix, original, stable_hashes, depot)
end

function run_path_dependency_process(fixture, phase; compiled = false, warmup = false)
    probe = joinpath(@__DIR__, "path_dependency_process.jl")
    mode = compiled ? `--compiled-modules=yes --pkgimages=existing` :
           `--compiled-modules=no --compile=min -O0`
    args = warmup ?
           `-e 'using BioAlignments; @assert BioAlignments.BLOSUM62[Char(65), Char(65)] == 4'` :
           `$probe $(fixture.root) $phase`
    command = `$(Base.julia_cmd()) --startup-file=no $mode --code-coverage=none --project=$(fixture.environment) $args`
    command = addenv(command,
        "JULIA_DEPOT_PATH" =>
            join([fixture.depot, DEPOT_PATH...], Sys.iswindows() ? ';' : ':'))
    log = joinpath(fixture.root, phase * ".log")
    process = open(log, "w") do io
        run(pipeline(ignorestatus(command); stdout = io, stderr = io))
    end
    success(process) || print(read(log, String))
    @test success(process)
    @test Dict(path => Iduna.Utils._file_sha256(path)
    for path in keys(fixture.stable_hashes)) ==
          fixture.stable_hashes
    return warmup ? nothing : JSON.parsefile(joinpath(fixture.root, phase * ".json"))
end

function edit_path_dependency_matrix(fixture)
    chmod(fixture.matrix, 0o644)
    write(fixture.matrix, replace(fixture.original, "A  4 -1" => "A -99 -1"; count = 1))
end
