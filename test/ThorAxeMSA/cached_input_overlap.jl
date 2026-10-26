function cached_input_tree_snapshot(root)
    snapshot = Dict{String, Any}()
    for (dir, dirs, files) in walkdir(root), name in vcat(dirs, files)

        path = joinpath(dir, name)
        snapshot[relpath(path, root)] = islink(path) ? (:link, readlink(path)) :
                                        isdir(path) ? (:directory, nothing) :
                                        (:file, read(path))
    end
    return snapshot
end

function check_cached_input_overlap_rejection(root, source, workdir, entrypoint, overwrite;
        recorded_state = false)
    thoraxe = Iduna.ThorAxeMSA
    dest = joinpath(workdir, "thoraxe_input")
    target = Iduna.ResolvedTarget(;
        input_id = "ENST000001", input_kind = :ensembl_transcript,
        ensembl_gene_id = "ENSG000001", transcript_id = "ENST000001")
    metadata = thoraxe._expected_transcript_query_metadata(target;
        specieslist = "homo_sapiens", orthology = "1:1", source_kind = "cached_input")
    expected = merge(metadata,
        (; source_path = abspath(source),
            source_fingerprint = thoraxe._bundle_fingerprint(source)))
    if isdir(dest)
        # This metadata-only cache was reused before strict manifest checks.
        thoraxe._write_transcript_query_metadata!(dest, expected)
    end
    stage_dir = thoraxe._thoraxe_input_stage_dir(workdir)
    if recorded_state
        thoraxe._write_thoraxe_input_state(workdir, dest, :done,
            thoraxe._thoraxe_input_identity(dest, expected); action = :run)
    else
        @test Iduna.Utils._read_stage_state(stage_dir) === nothing
    end
    if isdir(workdir)
        # Stop a buggy public call before any biological tool can run.
        write(joinpath(workdir, "thoraxe_msa"), "Intentional MSA stage blocker.\n")
        write(joinpath(workdir, "sibling.txt"), "Preserve unrelated work.\n")
    end
    before = cached_input_tree_snapshot(root)
    err = try
        if entrypoint === :public
            thoraxe.build_thoraxe_msa(target, workdir;
                cached_thoraxe_input_dir = source, overwrite,
                pid_thresholds = [10.0], pid_sample_count = 0, pid_sample_seed = 7,
                specieslist = "homo_sapiens", specieslist_filter = false,
                biomart_datasets_filter = false)
        else
            thoraxe._ensure_transcript_query(target, workdir;
                cached_input_dir = source, specieslist = "homo_sapiens", overwrite)
        end
        nothing
    catch exception
        exception
    end
    @test err isa ArgumentError
    @test occursin("overlap", sprint(showerror, err))
    @test occursin(source, sprint(showerror, err))
    @test occursin(dest, sprint(showerror, err))
    unchanged = cached_input_tree_snapshot(root) == before
    @test unchanged
    @test thoraxe._has_valid_ensembl_bundle(source)
    @test isfile(joinpath(source, "user-notes.txt")) &&
          read(joinpath(source, "user-notes.txt"), String) == "Preserve supplied source.\n"
end

@testset "cached input rejects a source inside its destination" begin
    for entrypoint in (:producer, :public), overwrite in (false, true),
        recorded_state in (false, true),
        layout in (:direct, :source_parent_alias, :workdir_parent_alias,
            :outward_source_alias, :outward_source_parent_alias)
        @testset "$entrypoint / overwrite=$overwrite / recorded=$recorded_state / $layout" begin
            mktempdir() do tmp
                physical = mkpath(joinpath(tmp, "physical"))
                alias = joinpath(tmp, "alias")
                symlink(physical, alias; dir_target = true)
                workdir = layout === :workdir_parent_alias ? alias : physical
                dest = write_test_ensembl_bundle(joinpath(workdir, "thoraxe_input"))
                source = layout === :source_parent_alias ?
                         joinpath(alias, "thoraxe_input", "original") :
                         joinpath(physical, "thoraxe_input", "original")
                if layout in (:outward_source_alias, :outward_source_parent_alias)
                    external = mkpath(joinpath(tmp, "external"))
                    symlink(external, source; dir_target = true)
                    if layout === :outward_source_parent_alias
                        # Reach the nested access path through another parent alias.
                        source = joinpath(alias, "thoraxe_input", "original", "bundle")
                    end
                end
                write_test_ensembl_bundle(source)
                write(joinpath(source, "user-notes.txt"), "Preserve supplied source.\n")
                @test !samefile(source, dest)
                check_cached_input_overlap_rejection(
                    tmp, source, workdir, entrypoint, overwrite; recorded_state)
            end
        end
    end
end

@testset "cached input rejects a destination inside its source" begin
    for entrypoint in (:producer, :public), overwrite in (false, true),
        layout in (:existing, :missing_destination, :missing_workdir_parent_alias,
            :destination_symlink)
        @testset "$entrypoint / overwrite=$overwrite / $layout" begin
            mktempdir() do tmp
                source = write_test_ensembl_bundle(joinpath(tmp, "source"))
                write(joinpath(source, "user-notes.txt"), "Preserve supplied source.\n")
                workdir = joinpath(source, "work")
                if layout === :missing_workdir_parent_alias
                    alias = joinpath(tmp, "alias")
                    symlink(source, alias; dir_target = true)
                    workdir = joinpath(alias, "missing", "work")
                else
                    mkpath(workdir)
                    if layout === :existing
                        write_test_ensembl_bundle(joinpath(workdir, "thoraxe_input"))
                    elseif layout === :destination_symlink
                        external = write_test_ensembl_bundle(joinpath(tmp, "external"))
                        symlink(external, joinpath(workdir, "thoraxe_input"); dir_target = true)
                    end
                end
                check_cached_input_overlap_rejection(tmp, source, workdir, entrypoint,
                    overwrite; recorded_state = layout in (:existing, :destination_symlink))
            end
        end
    end
end

@testset "cached input allows disjoint directories" begin
    for layout in (:prefix_sibling, :missing_workdir)
        mktempdir() do tmp
            workdir = layout === :prefix_sibling ? tmp : joinpath(tmp, "new", "work")
            source = write_test_ensembl_bundle(joinpath(tmp, "thoraxe_input_archive"))
            write(joinpath(source, "user-notes.txt"), "Preserve supplied source.\n")
            before = cached_input_tree_snapshot(source)
            target = Iduna.ResolvedTarget(;
                input_id = "ENST000001", input_kind = :ensembl_transcript,
                ensembl_gene_id = "ENSG000001", transcript_id = "ENST000001")
            dest = joinpath(workdir, "thoraxe_input")
            @test Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
                cached_input_dir = source, specieslist = "homo_sapiens") == dest
            @test cached_input_tree_snapshot(source) == before
            @test Iduna.ThorAxeMSA._has_valid_ensembl_bundle(dest)
            @test read(joinpath(dest, "user-notes.txt"), String) ==
                  "Preserve supplied source.\n"
        end
    end
end
