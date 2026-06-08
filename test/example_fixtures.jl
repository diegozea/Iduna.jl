module ExampleFixtures

using Iduna

export annotated_stockholm,
       artifact_path,
       iduna_args_after_app,
       iduna_args_after_module,
       iduna_stage_kwargs

function artifact_path(path::AbstractString, workdir::AbstractString)
    return Iduna.Utils._resolve_artifact_path(path, workdir)
end

function _resolved_target(input_id::AbstractString, workdir::AbstractString)
    sequence_path = joinpath(workdir, "sequences", "uniprot", "$(input_id).fasta")
    mkpath(dirname(sequence_path))
    write(sequence_path, ">$(input_id)\nAC\n")
    return Iduna.ResolvedTarget(;
        input_id,
        input_kind = Iduna.Utils.id_kind(input_id),
        uniprot_id = Iduna.Utils.is_uniprot_id(input_id) ? input_id : nothing,
        ensembl_gene_id = "ENSG00000120948.20",
        transcript_id = Iduna.Utils.is_ensembl_transcript_id(input_id) ?
                        input_id : "ENST00000240185.8",
        species = "homo_sapiens",
        uniprot_sequence_path = sequence_path)
end

function _thoraxe_msa(_target, workdir::AbstractString; kwargs...)
    pid = 10.0
    pid_dir = Iduna.Utils.format_pid_dir(pid)
    candidate_dir = joinpath(workdir, "thoraxe_msa", "candidates", pid_dir)
    sequence_dir = joinpath(candidate_dir, "sequences")
    species_dir = joinpath(candidate_dir, "species")
    mkpath(sequence_dir)
    mkpath(species_dir)

    seed_stockholm = joinpath(candidate_dir, "candidate_msa_full.sto")
    seed_fasta = joinpath(candidate_dir, "candidate_msa_full.fasta")
    sequence_fasta = joinpath(sequence_dir, "candidate_sequences_full.fasta")
    species_file = joinpath(species_dir, "candidate_species_full.txt")
    blocks_path = joinpath(candidate_dir, "candidate_msa_full_s_exon_blocks.tsv")
    summary_path = joinpath(workdir, "thoraxe_msa", "candidate_summary.csv")

    write(seed_stockholm, "# STOCKHOLM 1.0\nseed AC\n//\n")
    write(seed_fasta, ">seed\nAC\n")
    write(sequence_fasta, ">seed\nAC\n")
    write(species_file, "homo_sapiens\n")
    write(blocks_path,
        "alignment\tpid\tcode\ts_exon_id\tstart_col\tend_col\tn_columns\n")
    write(summary_path,
        "pid,selected,median_identity,mean_identity,stockholm_path,fasta_path\n" *
        "$(pid),true,1.0,1.0,$(seed_stockholm),$(seed_fasta)\n")

    seed = Iduna.SeedSelection(;
        pid,
        median_identity = 1.0,
        mean_identity = 1.0,
        stockholm_path = seed_stockholm,
        fasta_path = seed_fasta,
        s_exon_blocks_tsv = blocks_path,
        summary_path)
    return Iduna.ThorAxeMSAResult(;
        input_dir = joinpath(workdir, "thoraxe_input"),
        thoraxe_dirs = [joinpath(workdir, "thoraxe_msa", "runs", pid_dir)],
        msa_dir = joinpath(workdir, "thoraxe_msa"),
        baseline_fastas = [seed_fasta],
        baseline_stockholms = [seed_stockholm],
        sequence_fastas = [sequence_fasta],
        species_files = [species_file],
        pid_summary = summary_path,
        seeds = [seed],
        logs_dir = joinpath(workdir, "logs", "thoraxe"),
        pid_sample_count = 1,
        pid_sample_fraction = 1.0,
        pid_sample_seed = UInt64(7),
        sampling_strategy = :common)
end

function _expansion(target, seed, workdir::AbstractString; kwargs...)
    run_dir = joinpath(workdir, "expansion", target.ensembl_gene_id,
        target.transcript_id, Iduna.Utils.format_pid_dir(seed.pid))
    expanded_dir = joinpath(run_dir, "expanded_msa")
    mkpath(expanded_dir)

    hits_fasta = joinpath(expanded_dir, "$(target.transcript_id)_hits_raw.fasta")
    full_stockholm = joinpath(expanded_dir, "$(target.transcript_id)_full.sto")
    match_stockholm = joinpath(expanded_dir, "$(target.transcript_id)_matchonly.sto")
    a3m_path = joinpath(expanded_dir, "$(target.transcript_id).a3m")
    blocks_path = joinpath(expanded_dir, "$(target.transcript_id)_s_exon_blocks.tsv")

    write(hits_fasta, ">seed\nAC\n")
    write(full_stockholm, "# STOCKHOLM 1.0\nseed AC\n//\n")
    write(match_stockholm, "# STOCKHOLM 1.0\nseed AC\n//\n")
    write(a3m_path, ">seed\nAC\n")
    write(blocks_path,
        "alignment\tpid\tcode\ts_exon_id\tstart_col\tend_col\tn_columns\n")

    return Iduna.ExpansionResult(;
        run_dir,
        seed_stockholm = seed.stockholm_path,
        seed_fasta = seed.fasta_path,
        hits_fasta,
        full_stockholm,
        match_stockholm,
        a3m_path,
        s_exon_blocks_tsv = blocks_path,
        db_dir = joinpath(run_dir, "dbs"),
        hmm_dir = joinpath(run_dir, "hmm"),
        logs_dir = joinpath(run_dir, "logs"),
        n_hits = 1,
        n_new_hits = 0)
end

function _validation(_target, seed, expansion, workdir::AbstractString)
    stats_path = joinpath(workdir, "validation", Iduna.Utils.format_pid_dir(seed.pid),
        "stats.csv")
    mkpath(dirname(stats_path))
    write(stats_path, "metric,value\n")
    return Iduna.ValidationResult(;
        stats_path,
        seed_nseq = 1,
        seed_ncol = 2,
        expanded_nseq = expansion === nothing ? nothing : 1,
        expanded_ncol = expansion === nothing ? nothing : 2)
end

function _forbidden_expansion(args...; kwargs...)
    error("This example should not run the MSA expansion stage.")
end

function iduna_stage_kwargs(; fail_on_expansion::Bool = false)
    return (;
        _resolve_target = (input_id; workdir, kwargs...) -> _resolved_target(
            input_id, workdir),
        _build_thoraxe_msa = _thoraxe_msa,
        _expand_msa = fail_on_expansion ? _forbidden_expansion : _expansion,
        _validate_results = _validation)
end

function iduna_args_after_module(command::AbstractVector{<:AbstractString})
    module_index = findfirst(==("Iduna"), command)
    module_index === nothing && error("Command does not contain the Iduna module.")
    return String.(command[(module_index + 1):end])
end

function iduna_args_after_app(command::AbstractVector{<:AbstractString})
    separator = findfirst(==("--"), command)
    if separator !== nothing
        return String.(command[(separator + 1):end])
    end
    app_index = findfirst(==("iduna"), command)
    app_index === nothing && error("Command does not contain the iduna app.")
    return String.(command[(app_index + 1):end])
end

function annotated_stockholm(path::AbstractString)
    write(path,
        "# STOCKHOLM 1.0\n" *
        "#=GF SExonCodeMap \"a\"=>\"1_0\",\"b\"=>\"12_2\"\n" *
        "seed ACDE\n" *
        "#=GC SExonCode aabb\n" *
        "//\n")
    return path
end

end
