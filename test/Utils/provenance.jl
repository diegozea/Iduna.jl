module ImplementationProvenanceTests

using Test, JSON, Iduna
using Iduna.Utils: SHA, Pkg
using Iduna.Utils: _source_fingerprint, _implementation_identity, _executable_identity,
                   _mmseqs_snapshot_files, _read_mmseqs_snapshot, _mmseqs_snapshot_identity
@testset "Implementation provenance" begin
    mktempdir() do root
        @test _source_fingerprint(root) === nothing
        @test_logs (:warn, r"source identity unavailable") @test _implementation_identity(;
            source_sha256 = nothing) === nothing
        write(joinpath(root, "Project.toml"), "name = \"Tiny\"\n")
        source = joinpath(mkpath(joinpath(root, "src")), "Tiny.jl")
        write(source, "tiny() = 1\n")
        before = _source_fingerprint(root)
        write(joinpath(root, "README.md"), "Unrelated documentation")
        @test _source_fingerprint(root) == before
        write(source, "tiny() = 2\n")
        @test _source_fingerprint(root) != before
        extension = joinpath(mkpath(joinpath(root, "ext")), "TinyExt.jl")
        write(extension, "extension() = 1\n")
        with_extension = _source_fingerprint(root)
        write(extension, "extension() = 2\n")
        @test _source_fingerprint(root) != with_extension
        for directory in ("src", "ext")
            asset = joinpath(mkpath(joinpath(root, directory, "data")), "scores")
            before_asset = _source_fingerprint(root)
            write(asset, "A 4\n")
            with_asset = _source_fingerprint(root)
            @test with_asset != before_asset
            write(asset, "A -99\n")
            @test _source_fingerprint(root) != with_asset
            rm(asset)
            @test _source_fingerprint(root) == before_asset
        end
        before_coverage = _source_fingerprint(root)
        write(source * ".cov", "generated Julia coverage\n")
        write(source * ".123.cov", "generated Julia coverage\n")
        @test _source_fingerprint(root) == before_coverage
        write(joinpath(root, "src", "data", "scores.cov"), "scientific data\n")
        @test _source_fingerprint(root) != before_coverage
        id = _implementation_identity((JSON, SHA); source_sha256 = before)
        @test id.source_sha256 == before
        @test id.julia == string(VERSION)
        @test only(filter(dep -> dep.name == "JSON", id.dependencies)).version ==
              string(Base.pkgversion(JSON))
        @test _implementation_identity((JSON, SHA); source_sha256 = "changed") != id
        executable = joinpath(root, "tool")
        @test _executable_identity(executable) === nothing
        write(executable, "abcd")
        original = _executable_identity(executable)
        @test original.sha256 == bytes2hex(SHA.sha256("abcd"))
        write(executable, "wxyz")
        @test _executable_identity(executable) != original
        link = joinpath(root, "tool-link")
        symlink(executable, link)
        @test _executable_identity(link) == _executable_identity(executable)
        chmod(executable, 0)
        try
            @test _executable_identity(executable) === nothing
        finally
            chmod(executable, 0o600)
        end
    end
end

@testset "Immutable MMseqs snapshots" begin
    mktempdir() do tmp
        db = joinpath(tmp, "db")
        manifest = db * ".iduna.json"
        for suffix in ("", "_seq", "_aln", "_h", "_seq_h"),
            tail in ("", ".dbtype", ".index")

            write(db * suffix * tail, "abcd")
        end
        function snapshot!()
            files = _mmseqs_snapshot_files(db, manifest)
            record = Dict("schema_version" => 1, "immutable" => true,
                "snapshot_id" => "tiny-test-build",
                "components" =>
                    Dict(file => bytes2hex(open(SHA.sha256, joinpath(tmp, file)))
                    for file in files))
            write(manifest, JSON.json(record))
            return record
        end
        @test_logs (:warn, r"no immutable") @test _mmseqs_snapshot_identity(db, nothing) ===
                                                  nothing
        snapshot!()
        first_identity = _read_mmseqs_snapshot(db, manifest)
        @test _read_mmseqs_snapshot(db, manifest) == first_identity
        for suffix in ("_seq", "_aln")
            write(db * suffix, "wxyz")
            snapshot!()
            next_identity = _read_mmseqs_snapshot(db, manifest)
            @test next_identity != first_identity
            first_identity = next_identity
        end
        for malformed in ("[]", "null", "{", "{}")
            write(manifest, malformed)
            @test_logs (:warn, r"invalid snapshot") @test _mmseqs_snapshot_identity(db, manifest) ===
                                                          nothing
        end
        for update in (x -> x["immutable"] = false,
            x -> x["schema_version"] = 2,
            x -> x["schema_version"] = true,
            x -> x["schema_version"] = 1.0,
            x -> x["snapshot_id"] = " ",
            x -> delete!(x["components"], "db_seq"),
            x -> x["components"]["db_seq"] = "not-a-hash")
            record = snapshot!()
            update(record)
            write(manifest, JSON.json(record))
            @test_throws ErrorException _read_mmseqs_snapshot(db, manifest)
        end
        snapshot!()
        mv(db, db * ".0")
        snapshot!()
        @test _read_mmseqs_snapshot(db, manifest).snapshot_id == "tiny-test-build"
        mv(db * ".0", db * ".1")
        snapshot!()
        @test_throws ErrorException _read_mmseqs_snapshot(db, manifest)
        mv(db * ".1", db * ".0")
        rm(db * "_h.index")
        snapshot!()
        @test_throws ErrorException _read_mmseqs_snapshot(db, manifest)
    end
    @test Iduna.Utils._require_mmseqs_components(
        ["db_β" * suffix * tail for suffix in ("", "_seq", "_aln", "_h", "_seq_h")
         for tail in (".0", ".dbtype", ".index")], "db_β") === nothing
end

module IdentityProbeTool
const artifact_dir = "/tmp/test-artifact-tree"
const command = Ref{Any}(`missing`)
probe() = command[] isa Exception ? throw(command[]) : command[]
end

@testset "Installed tool identities" begin
    Utils = Iduna.Utils
    source = Utils._loaded_source_identity()
    @test source == Utils._source_fingerprint(pkgdir(Iduna))
    @test Utils._loaded_source_identity() == source
    base = Utils._implementation_identity()
    @test issubset(
        Set([
            "BioAlignments", "BioSequences", "MIToS", "StatsBase", "Statistics", "Random",
            "CSV", "DataFrames", "CodecZlib", "JSON", "HTTP"]),
        Set(dep.name for dep in base.dependencies))
    @test Utils._implementation_identity((JSON, JSON)) == base
    @test any(dep -> dep.name == "FastaIO", base.dependencies)
    @test any(dep -> dep.name == "AliasTables", base.dependencies)
    @test_logs (:warn, r"Julia dependency identity unavailable") @test Utils._implementation_identity((IdentityProbeTool,)) ===
                                                                       nothing
    withenv("OMP_NUM_THREADS" => "3") do
        @test Utils._implementation_identity().thread_environment["OMP_NUM_THREADS"] == "3"
    end
    for (mod, names) in ((Iduna.EPLI.HHsuite_jll, (:hhmake, :hhalign)),
        (Iduna.MSAExpansion.HMMER_jll, (:hmmbuild, :hmmalign, :esl_alimask)),
        (Iduna.MSAExpansion.MMseqs2_jll, (:mmseqs,)))
        actual = Utils._jll_tool_identity(mod, names)
        @test actual !== nothing
        @test all(isfile(value.path) for value in values(actual.executables))
    end
    @test Utils._jll_tool_identity(IdentityProbeTool, (:probe,)) === nothing
    IdentityProbeTool.command[] = ErrorException("probe failure")
    @test_logs (:warn, r"External tool identity unavailable") @test Utils._jll_tool_identity(
        IdentityProbeTool, (:probe,)) === nothing
    IdentityProbeTool.command[] = InterruptException()
    @test_throws InterruptException Utils._jll_tool_identity(IdentityProbeTool, (:probe,))
    mktempdir() do dir
        tool = joinpath(dir, "tool")
        write(tool, "abcd")
        IdentityProbeTool.command[] = Cmd([tool])
        @test Utils._jll_tool_identity(IdentityProbeTool, (:probe,)) === nothing
        write(tool, "wxyz")
        @test Utils._jll_tool_identity(IdentityProbeTool, (:probe,)) === nothing
    end
end

const selected = Pkg.dependencies()
const common = (
    Iduna.Utils.BioAlignments,
    Iduna.Utils.BioSequences,
    Iduna.Utils.MIToS,
    Iduna.Utils.StatsBase,
    Iduna.Utils.Statistics,
    Iduna.Utils.Random,
    Iduna.Utils.CSV,
    Iduna.Utils.DataFrames,
    Iduna.Utils.CodecZlib,
    Iduna.Utils.JSON,
    Iduna.Utils.HTTP
)

@testset "Selected and loaded dependency identities" begin
    actual = Iduna.Utils._julia_dependency_identity(common)
    @test actual !== nothing
    @test actual.julia == string(VERSION)
    byname = Dict(record.name => record for record in actual.packages)
    @test byname["FastaIO"].version ==
          string(selected[Base.UUID("a0c94c4b-ebed-5953-b5fc-82fe598ac79f")].version)
    @test haskey(byname, "AliasTables")
    @test length(actual.packages) > length(common)
    @test issorted([record.uuid for record in actual.packages])
    @test length(unique(record.uuid for record in actual.packages)) ==
          length(actual.packages)
    @test Iduna.Utils._julia_dependency_identity(reverse(common)) == actual
    @test byname["Random"].julia == string(VERSION)
    @test byname["Statistics"].tree_hash !== nothing
    local_package = Iduna.Utils._julia_dependency_identity((Iduna,))
    @test only(filter(record -> record.name == "Iduna", local_package.packages)).source_sha256 !==
          nothing
    extra = Iduna.Utils._julia_dependency_identity((common..., Iduna.ThorAxeMSA.ThorAxe))
    @test any(record -> record.name == "ThorAxe", extra.packages)
    println(
        "Saved environment: ",
        length(actual.packages),
        " common dependencies; ",
        length(extra.packages),
        " with ThorAxe."
    )

    # A tiny independent DAG: a real loaded root -> an unloaded registry leaf -> a path leaf.
    root = Iduna.Utils.BioAlignments
    uuid = Base.PkgId(root).uuid
    leaf = Base.UUID("ffffffff-0000-0000-0000-000000000001")
    pathleaf = Base.UUID("ffffffff-0000-0000-0000-000000000002")
    info = selected[uuid]
    info_tuple = (;
        (field => getfield(info, field) for field in fieldnames(typeof(info)))...)
    mktempdir() do tmp
        write(joinpath(tmp, "Project.toml"), "name = \"LocalLeaf\"\nversion = \"1.0.0\"\n")
        mkpath(joinpath(tmp, "src"))
        source = joinpath(tmp, "src", "LocalLeaf.jl")
        write(source, "const value = 1\n")
        graph = Dict{Base.UUID, Any}(
            uuid => merge(info_tuple, (; dependencies = Dict("RegistryLeaf" => leaf))),
            leaf => (;
                name = "RegistryLeaf",
                version = v"1.0.0",
                tree_hash = "tree-a",
                is_tracking_path = false,
                source = tmp,
                dependencies = Dict("LocalLeaf" => pathleaf)
            ),
            pathleaf => (;
                name = "LocalLeaf",
                version = v"1.0.0",
                tree_hash = nothing,
                is_tracking_path = true,
                source = tmp,
                dependencies = Dict{String, Base.UUID}()
            )
        )
        initial = Iduna.Utils._julia_dependency_identity((root,); graph)
        @test initial !== nothing
        @test Set(record.name for record in initial.packages) ==
              Set(["BioAlignments", "RegistryLeaf", "LocalLeaf"])
        registry_record = only(filter(record -> record.name == "RegistryLeaf", initial.packages))
        @test registry_record.version == "1.0.0"
        @test registry_record.tree_hash == "tree-a"
        changed = copy(graph)
        changed[leaf] = merge(graph[leaf], (; version = v"1.1.0"))
        @test Iduna.Utils._julia_dependency_identity((root,); graph = changed) != initial
        changed[leaf] = merge(graph[leaf], (; tree_hash = "tree-b"))
        @test Iduna.Utils._julia_dependency_identity((root,); graph = changed) != initial
        write(source, "const value = 2\n")
        @test Iduna.Utils._julia_dependency_identity((root,); graph) != initial
        write(source, "const value = 1\n")
        @test Iduna.Utils._julia_dependency_identity((root,); graph) == initial
        matrix = joinpath(mkpath(joinpath(tmp, "src", "data")), "matrix")
        write(matrix, "A 4\n")
        with_matrix = Iduna.Utils._julia_dependency_identity((root,); graph)
        @test with_matrix != initial
        write(matrix, "A -99\n")
        @test Iduna.Utils._julia_dependency_identity((root,); graph) != with_matrix
        rm(matrix)
        @test Iduna.Utils._julia_dependency_identity((root,); graph) == initial
        changed[leaf] = merge(graph[leaf], (; tree_hash = nothing))
        @test Iduna.Utils._julia_dependency_identity((root,); graph = changed) === nothing
        changed[leaf] = merge(graph[leaf], (; version = nothing))
        @test Iduna.Utils._julia_dependency_identity((root,); graph = changed) === nothing
        @test Iduna.Utils._julia_dependency_identity((root,); graph = Dict()) === nothing
        for patch in ((version = v"99.0.0",), (source = tmp,), (name = "WrongRoot",))
            changed = copy(graph)
            changed[uuid] = merge(graph[uuid], patch)
            @test Iduna.Utils._julia_dependency_identity((root,); graph = changed) ===
                  nothing
        end
        missing_source = copy(graph)
        missing_source[pathleaf] = merge(graph[pathleaf], (;
            source = joinpath(tmp, "absent")))
        @test Iduna.Utils._julia_dependency_identity((root,); graph = missing_source) ===
              nothing
        broken = copy(graph)
        delete!(broken, pathleaf)
        @test Iduna.Utils._julia_dependency_identity((root,); graph = broken) === nothing
        # Cycles and diamonds terminate and each package is recorded once.
        cyclic = copy(graph)
        cyclic[pathleaf] = merge(
            graph[pathleaf],
            (; dependencies = Dict("Root" => uuid, "RegistryLeaf" => leaf))
        )
        @test Iduna.Utils._julia_dependency_identity((root,); graph = cyclic) == initial
    end
    # Detect in-process environment updates even when a direct root is unchanged.
    stale = copy(selected)
    stats_id = Base.PkgId(Iduna.Utils.StatsBase).uuid
    alias_id = selected[stats_id].dependencies["AliasTables"]
    alias_info = selected[alias_id]
    alias_fields = (;
        (
        field => getfield(alias_info, field)
    for
    field in fieldnames(typeof(alias_info))
    )...
    )
    mixed = Dict{Base.UUID, Any}(stale)
    mixed[alias_id] = merge(alias_fields, (; version = v"99.0.0"))
    @test Iduna.Utils._julia_dependency_identity((Iduna.Utils.StatsBase,); graph = mixed) ===
          nothing
    mismatched_stdlib = Dict{Base.UUID, Any}(selected)
    random_uuid = Base.PkgId(Iduna.Utils.Random).uuid
    random_info = selected[random_uuid]
    @test Iduna.Utils._is_bundled_stdlib(random_uuid, random_info)
    @test !Iduna.Utils._is_bundled_stdlib(Base.PkgId(Iduna.Utils.Statistics).uuid,
        selected[Base.PkgId(Iduna.Utils.Statistics).uuid])
    @test !Iduna.Utils._is_bundled_stdlib(random_uuid,
        (; is_tracking_path = true, tree_hash = nothing))
    random_fields = (;
        (field => getfield(random_info, field)
    for field in fieldnames(typeof(random_info)))...)
    mismatched_stdlib[random_uuid] = merge(random_fields, (;
        source = "/nonexistent/wrong/Random"))
    @test Iduna.Utils._julia_dependency_identity((Iduna.Utils.Random,); graph = mismatched_stdlib) ===
          nothing
    # Stdlib records need no source hash; loaded-module consistency is checked separately.
    std_record = (;
        name = random_info.name,
        version = random_info.version,
        tree_hash = nothing,
        is_tracking_path = false,
        source = "/nonexistent/stdlib/source",
        dependencies = Dict{String, Base.UUID}()
    )
    std = Iduna.Utils._dependency_record(random_uuid, std_record)
    @test std !== nothing
    @test std.julia == string(VERSION)
    @test std.source_sha256 === nothing
    @test std.tree_hash === nothing
end

include("source_asset_symlinks.jl")
include("protein_alignment_provenance.jl")
include("path_dependency_provenance.jl")
include("compiled_dependency_provenance.jl")

end
