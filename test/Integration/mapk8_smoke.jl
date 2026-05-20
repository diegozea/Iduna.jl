@testset "MAPK8 live integration smoke test" begin
    function with_exponential_retry(f; attempts::Integer = 3)
        retry(f;
            delays = Base.ExponentialBackOff(;
                n = attempts - 1,
                first_delay = 10.0,
                max_delay = 60.0,
                factor = 2.0,
                jitter = 0.2))()
    end

    function write_mapk8_smoke_fasta(path::AbstractString)
        mkpath(dirname(path))
        open(path, "w") do io
            println(io, ">mapk8_human_smoke")
            println(io,
                "MSRSKRDNNFYSVEIGDSTFTVLKRYQNLKPIGSGAQGIVCAAYDAILERNVAIKKLSRPFQNQTHAKRAYRELVLMKCVNHKNIIGLLNVFTPQKSLEEFQDVYIVMELMDANLCQVIQ")
            println(io, ">mapk8_mouse_smoke")
            println(io,
                "MSRSKRDNNFYSVEIGDSTFTVLKRYQNLKPIGSGAQGIVCAAYDAILERNVAIKKLSRPFQNQTHAKRAYRELVLMKCVNHKNIIGLLNVFTPQKSLEEFQDVYIVMELMDANLCQVV")
        end
        return path
    end

    function build_mmseqs_smoke_db(tmp::AbstractString)
        mmseqs = Iduna.MSAExpansion.MMseqs2_jll.mmseqs()
        fasta = write_mapk8_smoke_fasta(joinpath(tmp, "mapk8_smoke.fasta"))
        db = joinpath(tmp, "mapk8_smoke_db")
        seq_db = string(db, "_seq")
        aln_db = string(db, "_aln")
        mmseqs_tmp = joinpath(tmp, "mmseqs_tmp")
        logs_dir = joinpath(tmp, "logs")
        mkpath(mmseqs_tmp)
        mkpath(logs_dir)

        Iduna.Utils.run_logged(`$(mmseqs) createdb $fasta $db`;
            stdout_path = joinpath(logs_dir, "createdb_stdout.log"),
            stderr_path = joinpath(logs_dir, "createdb_stderr.log"))
        Iduna.Utils.run_logged(`$(mmseqs) createdb $fasta $seq_db`;
            stdout_path = joinpath(logs_dir, "createdb_seq_stdout.log"),
            stderr_path = joinpath(logs_dir, "createdb_seq_stderr.log"))
        Iduna.Utils.run_logged(
            `$(mmseqs) search $db $seq_db $aln_db $mmseqs_tmp -a --threads 1`;
            stdout_path = joinpath(logs_dir, "search_stdout.log"),
            stderr_path = joinpath(logs_dir, "search_stderr.log"))
        return db
    end

    function checked_mapk8_smoke_result(result)
        artifact(path) = Iduna.Utils._resolve_artifact_path(path, result.workdir)
        required_files = (
            joinpath(result.workdir, "target.json"),
            joinpath(result.workdir, "result.json"),
            artifact(result.thoraxe_msa.pid_summary),
            artifact(result.thoraxe_msa.seeds[1].stockholm_path),
            artifact(result.expansions[1].match_stockholm),
            artifact(result.expansions[1].full_stockholm),
            artifact(result.expansions[1].a3m_path),
            artifact(result.expansions[1].hits_fasta),
            artifact(result.validations[1].stats_path)
        )
        missing = filter(!isfile, required_files)
        isempty(missing) || error("MAPK8 smoke test missing files: $(join(missing, ", "))")

        seed_msa = load_seed_msa(result)
        expanded_msa = load_expanded_msa(result)
        return (;
            status = result.status,
            uniprot_id = result.target.uniprot_id,
            transcript_core = Iduna.Utils.strip_ensembl_version(result.target.transcript_id),
            gene_core = Iduna.Utils.strip_ensembl_version(result.target.ensembl_gene_id),
            seed_nseq = Iduna.ResultsValidation.nsequences(seed_msa),
            expanded_nseq = Iduna.ResultsValidation.nsequences(expanded_msa)
        )
    end

    checked = with_exponential_retry() do
        mktempdir() do tmp
            mmseqs_db = build_mmseqs_smoke_db(tmp)
            result = iduna("ENST00000374179";
                uniprot_id = "P45983",
                specieslist = "homo_sapiens,mus_musculus",
                specieslist_filter = false,
                biomart_datasets_filter = false,
                workdir = joinpath(tmp, "iduna_mapk8"),
                mmseqs_db,
                overwrite = true,
                pid_thresholds = [10.0],
                threads = 1,
                transcript_query_retries = 3,
                transcript_query_timeout_seconds = 240,
                transcript_query_timeout_max_seconds = 480,
                thoraxe_timeout_seconds = 480)
            return checked_mapk8_smoke_result(result)
        end
    end

    @test checked.status in (:ok, :warn)
    @test checked.uniprot_id == "P45983"
    @test checked.transcript_core == "ENST00000374179"
    @test checked.gene_core == "ENSG00000107643"
    @test checked.seed_nseq >= 1
    @test checked.expanded_nseq >= 1
end
