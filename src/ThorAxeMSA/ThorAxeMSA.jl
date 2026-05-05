module ThorAxeMSA

using CSV
using DataFrames
using HHsuite_jll
using MIToS.MSA
using Statistics: mean, median
using ThorAxe

using ..Utils: DEFAULT_PID_THRESHOLDS, ResolvedTarget, SeedSelection, ThorAxeMSAResult,
               fasta_sequence, format_pid, protein_alignment_stats, resolve_sequence_name,
               safe_rm, write_fasta

export assemble_transcript_msa,
       build_thoraxe_msa,
       compute_identity_against_reference,
       select_best_seed

struct _CommandTimeoutError <: Exception
    command::String
    seconds::Float64
    stdout_log::String
    stderr_log::String
end

function Base.showerror(io::IO, err::_CommandTimeoutError)
    print(io,
        "ThorAxe command timed out after $(err.seconds) seconds: $(err.command). ",
        "See $(err.stdout_log) and $(err.stderr_log).")
end

const _REQUIRED_ENSEMBL_FILES = (
    "sequences.fasta",
    "exonstable.tsv",
    "tree.nh",
    "ensembl_version.csv",
    "tsl.csv"
)

_normalize_species_name(species::Nothing) = nothing
function _normalize_species_name(species::AbstractString)
    lowercase(replace(strip(String(species)), ' ' => '_'))
end

_thoraxe_input_dir(workdir::AbstractString) = joinpath(workdir, "thoraxe_input")
_thoraxe_output_dir(workdir::AbstractString) = joinpath(workdir, "thoraxe")
_thoraxe_msa_dir(workdir::AbstractString) = joinpath(workdir, "thoraxe_msa")
_thoraxe_seed_dir(workdir::AbstractString) = joinpath(_thoraxe_msa_dir(workdir), "seeds")
_thoraxe_logs_dir(workdir::AbstractString) = joinpath(workdir, "logs", "thoraxe")

function _has_valid_ensembl_bundle(bundle_root::AbstractString)::Bool
    ensembl_dir = joinpath(bundle_root, "Ensembl")
    isdir(ensembl_dir) || return false
    for file in _REQUIRED_ENSEMBL_FILES
        path = joinpath(ensembl_dir, file)
        isfile(path) && filesize(path) > 0 || return false
    end
    return true
end

function _open_logs(f::Function, stdout_path::AbstractString, stderr_path::AbstractString)
    mkpath(dirname(stdout_path))
    mkpath(dirname(stderr_path))
    open(stdout_path, "w") do out_io
        open(stderr_path, "w") do err_io
            f(out_io, err_io)
        end
    end
end

function _normalize_timeout(timeout_seconds::Union{Nothing, Real})
    timeout_seconds === nothing && return nothing
    value = Float64(timeout_seconds)
    value > 0 || error("ThorAxe timeout values must be positive.")
    return value
end

function _wait_logged_command(process, timeout_seconds::Union{Nothing, Float64})
    if timeout_seconds === nothing
        wait(process)
        return false
    end

    deadline = time() + timeout_seconds
    while process_running(process)
        if time() >= deadline
            kill(process)
            wait(process)
            return true
        end
        sleep(0.25)
    end
    wait(process)
    return false
end

function _run_logged_command(command::Cmd,
        stdout_log::AbstractString,
        stderr_log::AbstractString;
        timeout_seconds::Union{Nothing, Real} = nothing)
    timeout = _normalize_timeout(timeout_seconds)
    _open_logs(stdout_log, stderr_log) do out_io, err_io
        process = run(pipeline(command; stdout = out_io, stderr = err_io), wait = false)
        timed_out = _wait_logged_command(process, timeout)
        timed_out && throw(_CommandTimeoutError(
            string(command), timeout::Float64, String(stdout_log), String(stderr_log)))
        success(process) ||
            error("ThorAxe command failed: $(command). See $(stderr_log).")
    end
    return nothing
end

function _thoraxe_runner(stdout_log::AbstractString,
        stderr_log::AbstractString;
        timeout_seconds::Union{Nothing, Real} = nothing)
    # ThorAxe.jl builds the exact CLI command. The runner only controls logging,
    # timeouts, and error messages around that command.
    return command -> _run_logged_command(
        command, stdout_log, stderr_log; timeout_seconds = timeout_seconds)
end

function _next_timeout(timeout_seconds::Union{Nothing, Real},
        timeout_max_seconds::Union{Nothing, Real})
    timeout_seconds === nothing && return nothing
    current = _normalize_timeout(timeout_seconds)::Float64
    timeout_max_seconds === nothing && return current
    return min(2 * current, _normalize_timeout(timeout_max_seconds)::Float64)
end

function _retry_wait_seconds(attempt::Integer)
    return min(30.0, 2.0^(attempt - 1))
end

function _normalized_specieslist(specieslist::Union{Nothing, AbstractString})
    specieslist === nothing && return nothing
    stripped = strip(String(specieslist))
    return isempty(stripped) ? nothing : stripped
end

function _run_transcript_query_once(gene_core::AbstractString,
        workdir::AbstractString,
        species::Union{Nothing, AbstractString},
        specieslist::Union{Nothing, AbstractString},
        stdout_log::AbstractString,
        stderr_log::AbstractString;
        timeout_seconds::Union{Nothing, Real})
    runner = _thoraxe_runner(stdout_log, stderr_log; timeout_seconds = timeout_seconds)
    cd(workdir) do
        if species === nothing
            ThorAxe.transcript_query(gene_core; specieslist = specieslist, runner = runner)
        else
            ThorAxe.transcript_query(
                gene_core; species = species, specieslist = specieslist, runner = runner)
        end
    end
    return nothing
end

function _ensure_cached_thoraxe_input(source_dir::AbstractString,
        workdir::AbstractString;
        overwrite::Bool = false)
    dest = _thoraxe_input_dir(workdir)
    if !overwrite && _has_valid_ensembl_bundle(dest)
        return dest
    end

    source = abspath(String(source_dir))
    _has_valid_ensembl_bundle(source) ||
        error("Cached ThorAxe input at $(source) is not a valid transcript_query bundle.")
    if abspath(dest) != source
        isdir(dest) && safe_rm(dest, workdir)
        mkpath(dirname(dest))
        cp(source, dest; force = true)
    end
    _has_valid_ensembl_bundle(dest) ||
        error("Copied ThorAxe input bundle at $(dest) is incomplete.")
    return dest
end

function _ensure_transcript_query(target::ResolvedTarget, workdir::AbstractString;
        specieslist::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        cached_input_dir::Union{Nothing, AbstractString} = nothing,
        timeout_seconds::Union{Nothing, Real} = 180,
        timeout_max_seconds::Union{Nothing, Real} = 240,
        max_retries::Integer = 2,
        allow_specieslist_timeout_fallback::Bool = true)
    if cached_input_dir !== nothing
        return _ensure_cached_thoraxe_input(cached_input_dir, workdir; overwrite)
    end

    input_dir = _thoraxe_input_dir(workdir)
    if !overwrite && _has_valid_ensembl_bundle(input_dir)
        return input_dir
    end

    isdir(input_dir) && safe_rm(input_dir, workdir)
    gene_core = first(split(target.ensembl_gene_id, "."; limit = 2))
    tmp_gene_dir = joinpath(workdir, gene_core)
    isdir(tmp_gene_dir) && safe_rm(tmp_gene_dir, workdir)

    logs = _thoraxe_logs_dir(workdir)
    stdout_log = joinpath(logs, "transcript_query_stdout.log")
    stderr_log = joinpath(logs, "transcript_query_stderr.log")
    species = _normalize_species_name(target.species)

    attempts = max(Int(max_retries), 1)
    active_specieslist = _normalized_specieslist(specieslist)
    current_timeout = timeout_seconds
    for attempt in 1:attempts
        isdir(tmp_gene_dir) && safe_rm(tmp_gene_dir, workdir)
        try
            _run_transcript_query_once(gene_core, workdir, species, active_specieslist,
                stdout_log, stderr_log; timeout_seconds = current_timeout)
            if _has_valid_ensembl_bundle(tmp_gene_dir)
                break
            elseif attempt < attempts
                @warn "transcript_query produced an invalid bundle; retrying." gene=gene_core attempt
                active_specieslist = nothing
                sleep(_retry_wait_seconds(attempt))
                continue
            end
        catch err
            if err isa _CommandTimeoutError && _has_valid_ensembl_bundle(tmp_gene_dir)
                @warn "transcript_query timed out but produced a usable Ensembl bundle." gene=gene_core attempt
                break
            end
            if err isa _CommandTimeoutError && attempt < attempts
                if active_specieslist !== nothing && allow_specieslist_timeout_fallback
                    @warn "transcript_query timed out with a species list; retrying without it." gene=gene_core attempt
                    active_specieslist = nothing
                    current_timeout = timeout_max_seconds === nothing ?
                                      timeout_seconds : timeout_max_seconds
                else
                    current_timeout = _next_timeout(current_timeout, timeout_max_seconds)
                end
                sleep(_retry_wait_seconds(attempt))
                continue
            end
            rethrow(err)
        end
    end

    _has_valid_ensembl_bundle(tmp_gene_dir) ||
        error("transcript_query did not create a valid Ensembl bundle at $(tmp_gene_dir). See $(stderr_log).")
    mv(tmp_gene_dir, input_dir; force = true)
    return input_dir
end

function _ensure_baseline_thoraxe(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString;
        specieslist::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        timeout_seconds::Union{Nothing, Real} = nothing)
    thoraxe_dir = _thoraxe_output_dir(workdir)
    if !overwrite && isfile(joinpath(thoraxe_dir, "path_table.csv")) &&
       isfile(joinpath(thoraxe_dir, "s_exon_table.csv"))
        return thoraxe_dir
    end

    isdir(thoraxe_dir) && safe_rm(thoraxe_dir, workdir)
    stdout_log = joinpath(_thoraxe_logs_dir(workdir), "baseline_stdout.log")
    stderr_log = joinpath(_thoraxe_logs_dir(workdir), "baseline_stderr.log")
    runner = _thoraxe_runner(stdout_log, stderr_log; timeout_seconds = timeout_seconds)
    ThorAxe.thoraxe(input_dir, workdir; specieslist = specieslist, runner = runner)

    isfile(joinpath(thoraxe_dir, "path_table.csv")) ||
        error("ThorAxe baseline did not write path_table.csv in $(thoraxe_dir). See $(stderr_log).")
    return thoraxe_dir
end

function _join_msas_consistently(
        lhs::MSAType, rhs::MSAType) where {MSAType <:
                                           AbstractMultipleSequenceAlignment}
    common = intersect(sequencenames(lhs), sequencenames(rhs))
    isempty(common) && error("Cannot join s-exon MSAs that share no sequence names.")
    return join_msas(lhs, rhs)
end

function assemble_transcript_msa(thoraxe_dir::AbstractString,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    path_table_path = joinpath(thoraxe_dir, "path_table.csv")
    isfile(path_table_path) ||
        error("Missing ThorAxe path_table.csv at $(path_table_path).")
    path_table = DataFrame(CSV.File(path_table_path))

    transcript_core = first(split(transcript_id, "."; limit = 2))
    matches = findall(eachrow(path_table)) do row
        ids = split(String(row.TranscriptIDCluster), "/")
        transcript_core in ids || String(transcript_id) in ids
    end
    isempty(matches) &&
        error("Transcript $(transcript_id) was not found in $(path_table_path).")
    length(matches) == 1 ||
        error("Transcript $(transcript_id) matched more than one ThorAxe path.")

    transcript_path = String(path_table.Path[only(matches)])
    exon_tokens = split(replace(transcript_path, "start/" => "", "/stop" => ""), "/")
    exon_files = String[]
    for exon in exon_tokens
        startswith(exon, "0_") && continue
        push!(exon_files, joinpath(thoraxe_dir, "msa", "msa_s_exon_$(exon).fasta"))
    end
    isempty(exon_files) && error("No s-exon FASTA files were found for $(transcript_id).")
    all(isfile, exon_files) || error("At least one expected s-exon MSA is missing.")

    exon_msas = [read_file(file, FASTA) for file in exon_files]
    transcript_msa = reduce(_join_msas_consistently, exon_msas)

    reference = resolve_sequence_name(transcript_msa, gene_id)
    reference === nothing &&
        error("Could not find $(gene_id) in the reconstructed transcript MSA.")
    setreference!(transcript_msa, reference)

    species = fill("unknown", nsequences(transcript_msa))
    s_exon_path = joinpath(thoraxe_dir, "s_exon_table.csv")
    if isfile(s_exon_path)
        table = DataFrame(CSV.File(s_exon_path))
        lookup = Dict(String(row.GeneID) => String(row.Species) for row in eachrow(table))
        species = [get(lookup, String(name), "unknown")
                   for name in sequencenames(transcript_msa)]
    end
    return transcript_msa, species
end

function _save_baseline_msa(workdir::AbstractString,
        msa::AbstractMultipleSequenceAlignment,
        species::AbstractVector{<:AbstractString};
        overwrite::Bool = false)
    seed_dir = _thoraxe_seed_dir(workdir)
    mkpath(seed_dir)
    fasta_path = joinpath(seed_dir, "msa_0.fasta")
    sto_path = joinpath(seed_dir, "msa_0.sto")
    sequences_path = joinpath(seed_dir, "sequences_0.fasta")
    species_path = joinpath(seed_dir, "species_0.txt")
    if !overwrite && all(isfile, (fasta_path, sto_path, sequences_path, species_path))
        return fasta_path, sto_path, sequences_path, species_path
    end

    write_file(fasta_path, msa, FASTA)
    write_file(sto_path, msa, Stockholm)
    write_fasta(sequences_path,
        [(String(name), replace(stringsequence(msa, name), '-' => ""))
         for name in sequencenames(msa)])
    open(species_path, "w") do io
        for item in species
            println(io, String(item))
        end
    end
    return fasta_path, sto_path, sequences_path, species_path
end

function _read_single_fasta_sequence(path::AbstractString)
    seq = fasta_sequence(read(path, String))
    seq === nothing && error("FASTA file $(path) did not contain a sequence.")
    return seq
end

function _extract_reference_sequence(msa::AbstractMultipleSequenceAlignment,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    for id in (gene_id, transcript_id)
        name = resolve_sequence_name(msa, id)
        if name !== nothing
            sequence = replace(stringsequence(msa, name), '-' => "", '.' => "")
            return name, uppercase(String(sequence))
        end
    end
    error("Could not find a reference sequence for $(gene_id) / $(transcript_id).")
end

function _compare_protein_sequences(query_seq::AbstractString, reference_seq::AbstractString)
    protein_alignment_stats(query_seq, reference_seq)
end

function _validate_transcript_translation(target::ResolvedTarget,
        msa::AbstractMultipleSequenceAlignment)
    target.uniprot_sequence_path === nothing && return String[]
    isfile(target.uniprot_sequence_path) || return String[
        "UniProt sequence file is missing; skipped ThorAxe transcript validation."]

    query_name,
    query_seq = _extract_reference_sequence(msa,
        target.ensembl_gene_id, target.transcript_id)
    reference_seq = _read_single_fasta_sequence(target.uniprot_sequence_path)
    stats = _compare_protein_sequences(query_seq, reference_seq)
    warnings = String[]
    if stats.insertions != 0 || stats.deletions != 0
        error(
            "ThorAxe transcript $(target.transcript_id) has indels versus UniProt $(target.uniprot_id) ",
            "(insertions=$(stats.insertions), deletions=$(stats.deletions), query=$(query_name)).")
    elseif stats.mismatches != 0
        push!(warnings,
            "ThorAxe transcript $(target.transcript_id) has $(stats.mismatches) substitutions versus UniProt $(target.uniprot_id).")
    end
    return warnings
end

function _get_codes(output::AbstractString)
    in_alignment = false
    query_line = true
    cols = 1:0
    indices = Int[]
    positions = Int[]
    codes = Char[]
    for line in split(output, '\n')
        if startswith(line, "Probab=")
            in_alignment = true
            continue
        end
        in_alignment || continue
        isempty(line) && continue
        if query_line
            m = match(r"\S+\s+(\d+)\s+(\S+)\s+(\d+)", line)
            if m !== nothing
                start = parse(Int, m.captures[1])
                seq = m.captures[2]
                stop = parse(Int, m.captures[3])
                cols = findfirst(seq, line)
                indices = Int[]
                for residue in seq
                    if residue == '-'
                        push!(indices, 0)
                    else
                        push!(indices, start)
                        start += 1
                    end
                end
                start - 1 == stop || error("Could not parse HHsuite alignment positions.")
            end
            query_line = false
        elseif startswith(line, "Confidence")
            query_line = true
        elseif occursin(r"^\s", line) && !isempty(indices) && !isempty(cols)
            append!(positions, indices)
            append!(codes, collect(line[cols]))
        end
    end
    return positions, codes
end

function _identity_from_codes(positions::Vector{Int}, codes::Vector{Char})
    seen = Dict{Int, Bool}()
    for (pos, code) in zip(positions, codes)
        pos == 0 && continue
        seen[pos] = get(seen, pos, false) || code == '|'
    end
    isempty(seen) && return 0.0
    return 100 * count(identity, values(seen)) / length(seen)
end

function compute_identity_against_reference(reference_fasta::AbstractString,
        sample_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        label::AbstractString = "seed")
    mktempdir() do tmp
        ref_hhm = joinpath(tmp, "reference.hhm")
        sample_hhm = joinpath(tmp, "sample.hhm")
        out_path = joinpath(tmp, "hhalign.out")
        run(pipeline(
            `$(HHsuite_jll.hhmake()) -add_cons -M 100 -i $reference_fasta -o $ref_hhm`,
            stdout = devnull, stderr = devnull))
        run(pipeline(
            `$(HHsuite_jll.hhmake()) -add_cons -M 100 -i $sample_fasta -o $sample_hhm`,
            stdout = devnull, stderr = devnull))
        run(pipeline(
            `$(HHsuite_jll.hhalign()) -glob -M 100 -i $ref_hhm -t $sample_hhm -o $out_path`,
            stdout = devnull, stderr = devnull))
        output = read(out_path, String)
        if logs_dir !== nothing
            mkpath(logs_dir)
            cp(out_path, joinpath(logs_dir, "$(label)_hhalign.out"); force = true)
        end
        positions, codes = _get_codes(output)
        return _identity_from_codes(positions, codes)
    end
end

function _generate_pid_seed(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid::Real,
        species_file::AbstractString;
        overwrite::Bool = false,
        timeout_seconds::Union{Nothing, Real} = nothing)
    seed_dir = _thoraxe_seed_dir(workdir)
    pid_label = format_pid(pid)
    fasta_path = joinpath(seed_dir, "thoraxe_pid$(pid_label)_msa_0.fasta")
    sto_path = joinpath(seed_dir, "thoraxe_pid$(pid_label)_msa_0.sto")
    if !overwrite && isfile(fasta_path) && isfile(sto_path)
        return fasta_path, sto_path
    end

    tmp_root = joinpath(_thoraxe_msa_dir(workdir), "tmp")
    mkpath(tmp_root)
    mktempdir(tmp_root; prefix = "thoraxe_pid$(pid_label)_") do tmp
        stdout_log = joinpath(_thoraxe_logs_dir(workdir), "pid$(pid_label)_stdout.log")
        stderr_log = joinpath(_thoraxe_logs_dir(workdir), "pid$(pid_label)_stderr.log")
        runner = _thoraxe_runner(stdout_log, stderr_log; timeout_seconds = timeout_seconds)
        ThorAxe.thoraxe(
            input_dir, tmp; identity = Float64(pid), specieslist = species_file,
            runner = runner)
        pid_msa,
        _ = assemble_transcript_msa(joinpath(tmp, "thoraxe"),
            target.ensembl_gene_id, target.transcript_id)
        write_file(fasta_path, pid_msa, FASTA)
        write_file(sto_path, pid_msa, Stockholm)
    end
    return fasta_path, sto_path
end

function _summarize_pid_scores(rows::Vector{NamedTuple}, path::AbstractString)
    df = DataFrame(rows)
    sort!(df, [:pid])
    CSV.write(path, df)
    return path
end

function select_best_seed(summary_path::AbstractString)
    df = DataFrame(CSV.File(summary_path))
    isempty(df) && error("Cannot select a best seed from an empty PID summary.")
    sort!(df, [:median_identity, :mean_identity, :pid]; rev = [true, true, false])
    row = first(eachrow(df))
    return SeedSelection(;
        pid = Float64(row.pid),
        median_identity = Float64(row.median_identity),
        mean_identity = Float64(row.mean_identity),
        stockholm_path = String(row.stockholm_path),
        fasta_path = row.fasta_path === missing ? nothing : String(row.fasta_path),
        summary_path
    )
end

function build_thoraxe_msa(target::ResolvedTarget, workdir::AbstractString;
        pid_thresholds::AbstractVector{<:Real} = DEFAULT_PID_THRESHOLDS,
        specieslist::Union{Nothing, AbstractString} = nothing,
        cached_thoraxe_input_dir::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        transcript_query_timeout_seconds::Union{Nothing, Real} = 180,
        transcript_query_timeout_max_seconds::Union{Nothing, Real} = 240,
        transcript_query_retries::Integer = 2,
        allow_specieslist_timeout_fallback::Bool = true,
        thoraxe_timeout_seconds::Union{Nothing, Real} = nothing)
    input_dir = _ensure_transcript_query(target, workdir;
        specieslist, overwrite, cached_input_dir = cached_thoraxe_input_dir,
        timeout_seconds = transcript_query_timeout_seconds,
        timeout_max_seconds = transcript_query_timeout_max_seconds,
        max_retries = transcript_query_retries,
        allow_specieslist_timeout_fallback)
    thoraxe_dir = _ensure_baseline_thoraxe(
        target, input_dir, workdir; specieslist, overwrite,
        timeout_seconds = thoraxe_timeout_seconds)
    msa,
    species = assemble_transcript_msa(thoraxe_dir, target.ensembl_gene_id, target.transcript_id)
    warnings = _validate_transcript_translation(target, msa)
    baseline_fasta, baseline_sto,
    sequence_fasta, species_file = _save_baseline_msa(workdir, msa, species; overwrite)

    summary_path = joinpath(_thoraxe_msa_dir(workdir), "best_seed.csv")
    score_rows = NamedTuple[]
    for pid in Float64.(pid_thresholds)
        fasta_path,
        sto_path = _generate_pid_seed(
            target, input_dir, workdir, pid, species_file;
            overwrite, timeout_seconds = thoraxe_timeout_seconds)
        identity = compute_identity_against_reference(baseline_fasta, fasta_path;
            logs_dir = joinpath(_thoraxe_logs_dir(workdir), "hhalign"),
            label = "pid$(format_pid(pid))")
        push!(score_rows,
            (;
                gene_id = target.ensembl_gene_id,
                transcript_id = target.transcript_id,
                pid,
                mean_identity = identity,
                median_identity = identity,
                n_samples = 1,
                n_sequences_msa0 = nsequences(msa),
                fasta_path,
                stockholm_path = sto_path
            ))
    end
    _summarize_pid_scores(score_rows, summary_path)
    best = select_best_seed(summary_path)

    return ThorAxeMSAResult(;
        input_dir,
        thoraxe_dir,
        msa_dir = _thoraxe_msa_dir(workdir),
        baseline_fasta,
        baseline_stockholm = baseline_sto,
        sequence_fasta,
        species_file,
        pid_summary = summary_path,
        best_seed = best,
        logs_dir = _thoraxe_logs_dir(workdir),
        warnings,
        status = :ok
    )
end

end
