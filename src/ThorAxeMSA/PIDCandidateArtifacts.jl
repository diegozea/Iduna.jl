# PID Candidate Artifacts
# -----------------------

function _candidate_pid_dir(workdir::AbstractString, pid::Real)
    joinpath(_thoraxe_candidates_dir(workdir), format_pid_dir(pid))
end

function _candidate_sample_label(sample_idx::Integer)
    sample_idx == 0 && return "full"
    return "species_subset_$(lpad(string(sample_idx), 3, '0'))"
end

function _pid_sample_paths(workdir::AbstractString, pid::Real, sample_idx::Integer)
    candidate_dir = _candidate_pid_dir(workdir, pid)
    label = _candidate_sample_label(sample_idx)
    msa_name = sample_idx == 0 ? "candidate_msa_full" : "candidate_msa_$(label)"
    sequence_name = sample_idx == 0 ? "candidate_sequences_full" :
                    "candidate_sequences_$(label)"
    species_name = sample_idx == 0 ? "candidate_species_full" :
                   "candidate_$(label)"
    return (;
        fasta_path = joinpath(candidate_dir, "$(msa_name).fasta"),
        stockholm_path = joinpath(candidate_dir, "$(msa_name).sto"),
        s_exon_blocks_tsv = joinpath(candidate_dir, "$(msa_name)_s_exon_blocks.tsv"),
        sequence_fasta = joinpath(candidate_dir, "sequences", "$(sequence_name).fasta"),
        species_file = joinpath(candidate_dir, "species", "$(species_name).txt")
    )
end

function _shared_sample_species_file(workdir::AbstractString, sample_idx::Integer)
    return joinpath(_thoraxe_sample_species_dir(workdir),
        "candidate_$(_candidate_sample_label(sample_idx)).txt")
end

function _pid_scores_path(workdir::AbstractString, pid::Real)
    joinpath(_candidate_pid_dir(workdir, pid), "scores.csv")
end

function _pid_sample_run_root(workdir::AbstractString, pid::Real, sample_idx::Integer)
    joinpath(_thoraxe_pid_runs_dir(workdir), format_pid_dir(pid),
        _candidate_sample_label(sample_idx))
end

function _pid_sample_thoraxe_dir(workdir::AbstractString, pid::Real, sample_idx::Integer)
    joinpath(_pid_sample_run_root(workdir, pid, sample_idx), "thoraxe")
end

function _write_species_file(path::AbstractString,
        species::AbstractVector{<:AbstractString};
        overwrite::Bool = false)
    if !overwrite && isfile(path)
        return path
    end
    mkpath(dirname(path))
    open(path, "w") do io
        for item in species
            println(io, String(item))
        end
    end
    return path
end

function _symlink_species_file(link_path::AbstractString,
        target_path::AbstractString;
        overwrite::Bool = false)
    if !overwrite && islink(link_path)
        current = readlink(link_path)
        current_path = isabspath(current) ? current : joinpath(dirname(link_path), current)
        abspath(current_path) == abspath(target_path) && return link_path
    end
    if !overwrite && isfile(link_path)
        return link_path
    end
    if !overwrite && ispath(link_path)
        return link_path
    end
    mkpath(dirname(link_path))
    ispath(link_path) && rm(link_path; force = true)
    # A relative link keeps the work directory movable.
    target = relpath(target_path, dirname(link_path))
    symlink(target, link_path)
    return link_path
end

function _write_candidate_sample_inputs(paths,
        msa::AbstractMultipleSequenceAlignment,
        species::AbstractVector{<:AbstractString},
        indices::AbstractVector{<:Integer};
        overwrite::Bool = false)
    if !overwrite && isfile(paths.sequence_fasta) && isfile(paths.species_file)
        return paths.sequence_fasta, paths.species_file
    end
    length(species) == nsequences(msa) ||
        error("Species list length does not match the MSA sequence count.")
    seq_names = collect(sequencenames(msa))
    names = String.(seq_names)
    mkpath(dirname(paths.sequence_fasta))
    write_fasta(paths.sequence_fasta,
        [(names[i], replace(stringsequence(msa, seq_names[i]), '-' => "", '.' => ""))
         for i in indices])
    _write_species_file(paths.species_file, String[species[i] for i in indices];
        overwrite = true)
    return paths.sequence_fasta, paths.species_file
end

function _write_seed_alignment_outputs(paths,
        msa::AbstractMultipleSequenceAlignment,
        pid::Real)
    write_file(paths.fasta_path, msa, FASTA)
    write_file(paths.stockholm_path, msa, Stockholm)
    write_s_exon_blocks_tsv(paths.s_exon_blocks_tsv, msa;
        alignment = "seed",
        pid = Float64(pid))
    return paths
end

function _ensure_seed_blocks_tsv(paths, pid::Real)
    isfile(paths.s_exon_blocks_tsv) && return paths.s_exon_blocks_tsv
    isfile(paths.stockholm_path) || return paths.s_exon_blocks_tsv
    msa = read_file(paths.stockholm_path, Stockholm; keepinserts = true)
    has_s_exon_annotations(msa) || return paths.s_exon_blocks_tsv
    write_s_exon_blocks_tsv(paths.s_exon_blocks_tsv, msa;
        alignment = "seed",
        pid = Float64(pid))
    return paths.s_exon_blocks_tsv
end

function _seed_has_s_exon_annotations(stockholm_path::AbstractString)
    isfile(stockholm_path) || return false
    msa = read_file(stockholm_path, Stockholm; keepinserts = true)
    return has_s_exon_annotations(msa)
end

function _reference_index(msa::AbstractMultipleSequenceAlignment,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    names = String.(sequencenames(msa))
    for id in (gene_id, transcript_id)
        name = resolve_sequence_name(msa, id)
        name === nothing && continue
        idx = findfirst(==(String(name)), names)
        idx === nothing || return idx
    end
    error("Could not find a reference sequence for $(gene_id) / $(transcript_id).")
end

function _species_sample(species::AbstractVector{<:AbstractString},
        reference_species::AbstractString,
        fraction::Real,
        rng::MersenneTwister)
    normalized_species = _unique_nonempty_species(String.(species))
    reference = String(reference_species)
    reference in normalized_species ||
        error("Reference species $(reference) is not present in the sampling universe.")
    selectable = [item for item in normalized_species if item != reference]
    isempty(selectable) && return [reference]
    n_keep = clamp(round(Int, Float64(fraction) * length(selectable)), 1, length(selectable))
    return vcat(reference, sample(rng, selectable, n_keep; replace = false))
end

function _reference_species(candidate, gene_id::AbstractString, transcript_id::AbstractString)
    reference_idx = _reference_index(candidate.msa, gene_id, transcript_id)
    reference_idx <= length(candidate.species) ||
        error("Reference sequence index is outside the candidate species list.")
    return String(candidate.species[reference_idx])
end

function _candidate_species_set(candidate)
    return Set(_unique_nonempty_species(String.(candidate.species)))
end

function _common_sampling_universe(records::AbstractVector, gene_id::AbstractString,
        transcript_id::AbstractString)
    isempty(records) && return (species = String[], reference_species = nothing)
    first_record = first(records)
    ordered = _unique_nonempty_species(String.(first_record.candidate.species))
    common = _candidate_species_set(first_record.candidate)
    for record in Iterators.drop(records, 1)
        intersect!(common, _candidate_species_set(record.candidate))
    end
    # Shared samples compare PID candidates using only species present in every candidate.
    species = [item for item in ordered if item in common]
    reference_species = _reference_species(first_record.candidate, gene_id, transcript_id)
    reference_species in species ||
        error("The eligible PID candidates share no reference species for common sampling.")
    return (; species, reference_species)
end

function _input_sampling_universe(target::ResolvedTarget,
        records::AbstractVector,
        effective_specieslist::Union{Nothing, AbstractString})
    species = _parse_specieslist(effective_specieslist)
    species === nothing &&
        error("sampling_strategy=:input requires a non-empty effective specieslist; pass specieslist or use sampling_strategy=:common.")
    if !isempty(records)
        reference_species = _reference_species(first(records).candidate,
            target.ensembl_gene_id, target.transcript_id)
    else
        query_species = _normalize_species_name(target.species)
        reference_species = query_species === nothing ? first(species) : query_species
    end
    reference_species in species ||
        error("sampling_strategy=:input requires the reference species $(reference_species) to be present in the effective specieslist.")
    return (; species, reference_species)
end

function _shared_sampling_universe(strategy::Symbol,
        target::ResolvedTarget,
        records::AbstractVector,
        effective_specieslist::Union{Nothing, AbstractString})
    if strategy === :common
        return _common_sampling_universe(records, target.ensembl_gene_id,
            target.transcript_id)
    elseif strategy === :input
        return _input_sampling_universe(target, records, effective_specieslist)
    end
    error("No shared sampling universe is defined for sampling_strategy=$(strategy).")
end

function _write_shared_species_samples(workdir::AbstractString,
        species::AbstractVector{<:AbstractString},
        reference_species::AbstractString;
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        overwrite::Bool = false)
    for sample_idx in 1:sample_count
        rng = _sample_rng(sample_seed, sample_idx)
        sampled_species = _species_sample(species, reference_species, sample_fraction, rng)
        _write_species_file(_shared_sample_species_file(workdir, sample_idx),
            sampled_species; overwrite)
    end
    return nothing
end

function _link_pid_candidate_samples(workdir::AbstractString,
        pids::AbstractVector{<:Real};
        sample_count::Integer,
        overwrite::Bool = false)
    for pid in pids
        for sample_idx in 1:sample_count
            paths = _pid_sample_paths(workdir, pid, sample_idx)
            _symlink_species_file(paths.species_file,
                _shared_sample_species_file(workdir, sample_idx); overwrite)
        end
    end
    return nothing
end

function _ensure_pid_candidate_samples(workdir::AbstractString,
        pid::Real,
        msa::AbstractMultipleSequenceAlignment,
        species::AbstractVector{<:AbstractString};
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        overwrite::Bool = false,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    reference_idx = _reference_index(msa, gene_id, transcript_id)
    @debug "Preparing ThorAxe PID species sample inputs." pid sample_count sample_fraction
    for sample_idx in 1:sample_count
        paths = _pid_sample_paths(workdir, pid, sample_idx)
        rng = _sample_rng(sample_seed, sample_idx)
        indices = _sample_indices(nsequences(msa), reference_idx, sample_fraction, rng)
        _write_candidate_sample_inputs(paths, msa, species, indices; overwrite)
    end
    return nothing
end

function _run_kept_thoraxe_pid_msa!(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        paths,
        thoraxe_dir::AbstractString,
        path_table::AbstractString,
        pid::Real,
        specieslist::Union{Nothing, AbstractString},
        sample_idx::Integer,
        runner,
        overwrite::Bool,
        thoraxe_fn::Function = ThorAxe.thoraxe)
    if overwrite || !isfile(path_table) || !_has_phylosofs_outputs(thoraxe_dir)
        run_root = _pid_sample_run_root(workdir, pid, sample_idx)
        isdir(run_root) && safe_rm(run_root, workdir)
        mkpath(dirname(run_root))
        thoraxe_fn(
            input_dir, run_root; identity = Float64(pid), specieslist = specieslist,
            phylosofs = true, runner = runner)
    end
    pid_msa,
    _ = assemble_transcript_msa(thoraxe_dir,
        target.ensembl_gene_id, target.transcript_id)
    _write_seed_alignment_outputs(paths, pid_msa, pid)
    return paths.fasta_path, paths.stockholm_path, thoraxe_dir
end

function _log_full_pid_msa(message::AbstractString,
        target::ResolvedTarget,
        pid::Real,
        sample_idx::Integer,
        extra = NamedTuple();
        level::Symbol = :info)
    sample_idx == 0 || return nothing
    details = (;
        gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        pid,
        sample = "full",
        extra...)
    if level === :debug
        @debug message details...
    else
        @info message details...
    end
    return nothing
end

function _run_thoraxe_pid_msa(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid::Real,
        specieslist::Union{Nothing, AbstractString},
        sample_idx::Integer;
        overwrite::Bool = false,
        keep_thoraxe_dir::Bool = false,
        thoraxe_fn::Function = ThorAxe.thoraxe)
    paths = _pid_sample_paths(workdir, pid, sample_idx)
    thoraxe_dir = _pid_sample_thoraxe_dir(workdir, pid, sample_idx)
    path_table = joinpath(thoraxe_dir, "path_table.csv")
    if !overwrite && isfile(paths.fasta_path) && isfile(paths.stockholm_path) &&
       _seed_has_s_exon_annotations(paths.stockholm_path) &&
       (!keep_thoraxe_dir || (isfile(path_table) && _has_phylosofs_outputs(thoraxe_dir)))
        _log_full_pid_msa("Reusing ThorAxe PID MSA.", target, pid, sample_idx,
            (; fasta_path = paths.fasta_path, stockholm_path = paths.stockholm_path))
        _ensure_seed_blocks_tsv(paths, pid)
        return paths.fasta_path, paths.stockholm_path, thoraxe_dir
    end

    mkpath(_candidate_pid_dir(workdir, pid))
    pid_label = format_pid(pid)
    sample_label = "pid$(pid_label)_sample$(sample_idx)"
    stdout_log = joinpath(_thoraxe_logs_dir(workdir), "$(sample_label)_stdout.log")
    stderr_log = joinpath(_thoraxe_logs_dir(workdir), "$(sample_label)_stderr.log")
    runner = _thoraxe_runner(stdout_log, stderr_log)
    _log_full_pid_msa("Running ThorAxe PID MSA command.", target, pid, sample_idx,
        (; _specieslist_log_details(specieslist)..., stdout_log, stderr_log);
        level = :debug)

    if keep_thoraxe_dir
        return _run_kept_thoraxe_pid_msa!(target, input_dir, workdir, paths,
            thoraxe_dir, path_table, pid, specieslist, sample_idx, runner, overwrite,
            thoraxe_fn)
    end

    tmp_root = joinpath(_thoraxe_msa_dir(workdir), "tmp")
    mkpath(tmp_root)
    mktempdir(tmp_root; prefix = "$(sample_label)_") do tmp
        thoraxe_fn(
            input_dir, tmp; identity = Float64(pid), specieslist = specieslist,
            phylosofs = true, runner = runner)
        pid_msa,
        _ = assemble_transcript_msa(joinpath(tmp, "thoraxe"),
            target.ensembl_gene_id, target.transcript_id)
        _write_seed_alignment_outputs(paths, pid_msa, pid)
    end
    return paths.fasta_path, paths.stockholm_path, thoraxe_dir
end

function _generate_pid_candidate(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid::Real,
        specieslist::Union{Nothing, AbstractString};
        overwrite::Bool = false,
        thoraxe_fn::Function = ThorAxe.thoraxe)
    @debug "Generating ThorAxe PID candidate." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id pid
    paths = _pid_sample_paths(workdir, pid, 0)
    fasta_path, sto_path,
    thoraxe_dir = _run_thoraxe_pid_msa(
        target, input_dir, workdir, pid, specieslist, 0;
        overwrite, keep_thoraxe_dir = true, thoraxe_fn)
    msa,
    species = assemble_transcript_msa(thoraxe_dir,
        target.ensembl_gene_id, target.transcript_id)
    _write_seed_alignment_outputs(paths, msa, pid)
    _write_candidate_sample_inputs(
        paths, msa, species, collect(1:nsequences(msa)); overwrite)
    return (; fasta_path, stockholm_path = sto_path, thoraxe_dir,
        s_exon_blocks_tsv = paths.s_exon_blocks_tsv,
        sequence_fasta = paths.sequence_fasta, species_file = paths.species_file,
        msa, species)
end
