@testset "cached input rebuild preserves the supplied bundle" begin
    for entrypoint in (:producer, :public),
        layout in (:source_parent_alias, :workdir_parent_alias, :identical, :distinct)

        @testset "$entrypoint / $layout" begin
            mktempdir() do tmp
                thoraxe = Iduna.ThorAxeMSA
                physical = mkpath(joinpath(tmp, "physical"))
                alias = joinpath(tmp, "alias")
                if layout in (:source_parent_alias, :workdir_parent_alias)
                    symlink(physical, alias; dir_target = true)
                end
                workdir = layout === :workdir_parent_alias ? alias : physical
                dest = joinpath(workdir, "thoraxe_input")
                source = layout === :source_parent_alias ?
                         joinpath(alias, "thoraxe_input") :
                         layout === :distinct ? joinpath(tmp, "source") :
                         joinpath(physical, "thoraxe_input")
                write_test_ensembl_bundle(source)
                write(joinpath(source, "user-notes.txt"), "Keep this scientific archive note.\n")
                target = Iduna.ResolvedTarget(;
                    input_id = "ENST000001",
                    input_kind = :ensembl_transcript,
                    ensembl_gene_id = "ENSG000001",
                    transcript_id = "ENST000001")
                metadata = thoraxe._expected_transcript_query_metadata(target;
                    specieslist = "homo_sapiens", orthology = "1:1",
                    source_kind = "cached_input")
                expected = merge(metadata,
                    (;
                        source_path = abspath(source),
                        source_fingerprint = thoraxe._bundle_fingerprint(source)))
                thoraxe._write_transcript_query_metadata!(source, expected)
                if layout === :distinct
                    cp(source, dest)
                    write(joinpath(dest, "stale-destination.txt"), "Replace this old copy.\n")
                end
                @test samefile(source, dest) == (layout !== :distinct)
                @test (abspath(source) == abspath(dest)) == (layout === :identical)
                source_metadata = read(thoraxe._transcript_query_metadata_path(source))
                hashes = Dict(relpath(joinpath(dir, file), source) =>
                                  Iduna.Utils._file_sha256(joinpath(dir, file))
                for (dir, _, files) in walkdir(source)
                for file in files
                if file != thoraxe._TRANSCRIPT_QUERY_METADATA_FILE)
                preserved(root) = all(isfile(joinpath(root, path)) &&
                                          Iduna.Utils._file_sha256(joinpath(root, path)) ==
                                          digest
                for (path, digest) in hashes)
                sibling = joinpath(workdir, "sibling.txt")
                write(sibling, "Unrelated workdir file.\n")
                stage_dir = thoraxe._thoraxe_input_stage_dir(workdir)
                @test Iduna.Utils._read_stage_state(stage_dir) === nothing
                identity = thoraxe._thoraxe_input_identity(dest, expected)
                @test !thoraxe._classify_thoraxe_input_stage(
                    workdir, dest, identity; overwrite = false).reusable

                if entrypoint === :public
                    # Stop after the real input stage, before any biological tool runs.
                    blocker = joinpath(workdir, "thoraxe_msa")
                    write(blocker, "Intentional MSA stage blocker.\n")
                    err = try
                        thoraxe.build_thoraxe_msa(target, workdir;
                            cached_thoraxe_input_dir = source, overwrite = false,
                            pid_thresholds = [10.0], pid_sample_count = 0,
                            pid_sample_seed = 7, specieslist = "homo_sapiens",
                            specieslist_filter = false, biomart_datasets_filter = false)
                        nothing
                    catch exception
                        exception
                    end
                    @test err isa Base.IOError
                    @test err isa Base.IOError && err.code == Base.UV_EEXIST
                    @test occursin(blocker, sprint(showerror, err))
                    @test read(blocker, String) == "Intentional MSA stage blocker.\n"
                else
                    result = try
                        thoraxe._ensure_transcript_query(target, workdir;
                            cached_input_dir = source, specieslist = "homo_sapiens",
                            overwrite = false)
                    catch exception
                        exception
                    end
                    @test result == dest
                end
                @test preserved(source)
                @test preserved(dest)
                @test read(sibling, String) == "Unrelated workdir file.\n"
                @test !isfile(joinpath(dest, "stale-destination.txt"))
                state = Iduna.Utils._read_stage_state(stage_dir)
                @test state["status"] == "done"
                @test state["action"] == "rebuild"
                state["status"] == "done" || return
                identity = thoraxe._thoraxe_input_identity(dest, expected)
                @test thoraxe._classify_thoraxe_input_stage(
                    workdir, dest, identity; overwrite = false).reusable
                @test thoraxe._ensure_transcript_query(target, workdir;
                    cached_input_dir = source, specieslist = "homo_sapiens",
                    overwrite = true) == dest
                @test Iduna.Utils._read_stage_state(stage_dir)["action"] == "rebuild"
                @test preserved(source)
                @test preserved(dest)
                written = thoraxe.JSON.parse(read(
                    thoraxe._transcript_query_metadata_path(dest), String))
                @test written["source_fingerprint"] == expected.source_fingerprint
                @test written["input_fingerprint"] == thoraxe._bundle_fingerprint(dest)
                if layout === :distinct
                    @test read(thoraxe._transcript_query_metadata_path(source)) ==
                          source_metadata
                end

                @test thoraxe._ensure_transcript_query(target, workdir;
                    cached_input_dir = source, specieslist = "homo_sapiens") == dest
                @test Iduna.Utils._read_stage_state(stage_dir)["action"] == "reuse"
                # A bad digest still forces rebuilding after physical-identity handling.
                if layout === :distinct
                    write(joinpath(dest, "Ensembl", "sequences.fasta"), "corrupt\n")
                else
                    write(thoraxe._transcript_query_metadata_path(dest), "{}")
                end
                @test !thoraxe._classify_thoraxe_input_stage(
                    workdir, dest, identity; overwrite = false).reusable
                @test thoraxe._ensure_transcript_query(target, workdir;
                    cached_input_dir = source, specieslist = "homo_sapiens") == dest
                @test Iduna.Utils._read_stage_state(stage_dir)["action"] == "rebuild"
                @test preserved(source)
                @test preserved(dest)
                @test thoraxe._classify_thoraxe_input_stage(
                    workdir, dest, identity; overwrite = false).reusable
            end
        end
    end
end
