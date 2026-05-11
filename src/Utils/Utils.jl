module Utils

using BioAlignments: AffineGapScoreModel, BLOSUM62, GlobalAlignment, alignment,
                     count_deletions, count_insertions, count_mismatches, pairalign
using BioSequences: LongAA
using CodecZlib: GzipDecompressor
import HTTP
import JSON3
using MIToS.MSA: AbstractMultipleSequenceAlignment, sequencenames
using Printf: @sprintf

include("Types.jl")

export DEFAULT_PID_THRESHOLDS,
       ExpansionResult,
       IdunaResult,
       ResolvedTarget,
       SeedSelection,
       ThorAxeMSAResult,
       ValidationResult,
       decode_body,
       ensure_mmseqs_db,
       fasta_sequence,
       format_fasta,
       id_kind,
       is_ensembl_transcript_id,
       is_uniprot_id,
       prepare_output_dir,
       protein_alignment_stats,
       resolve_sequence_name,
       result_summary,
       run_logged,
       safe_rm,
       sequence_name_variants,
       strip_ensembl_version,
       write_fasta,
       write_json,
       write_text,
       format_pid

const DEFAULT_PID_THRESHOLDS = Float64[10, 20, 30, 60, 80]
const _PROTEIN_ALIGNMENT_SCORE_MODEL = AffineGapScoreModel(BLOSUM62, gap_open = -10, gap_extend = -1)
const _TRANSIENT_HTTP_STATUSES = Set([429, 500, 502, 503, 504])

struct _RetryableHTTPStatus <: Exception
    response::HTTP.Response
end

# Julia's retry helper retries exceptions, so transient statuses are wrapped.
_is_transient_http_status(status::Integer)::Bool = status in _TRANSIENT_HTTP_STATUSES
_is_retryable_http_exception(_state, err)::Bool = err isa _RetryableHTTPStatus

function _http_get_with_retries(url::AbstractString,
        headers;
        retries::Integer = 4,
        sleep_seconds::Real = 1.5,
        max_delay::Real = 30.0,
        http_get::Function = HTTP.get)
    attempts = max(Int(retries), 1)
    try
        # Return the final HTTP response so callers decide what status means success.
        return retry(;
            delays = Base.ExponentialBackOff(;
                n = attempts - 1,
                first_delay = Float64(sleep_seconds),
                max_delay = Float64(max_delay),
                factor = 2.0,
                jitter = 0.0),
            check = _is_retryable_http_exception) do
            resp = http_get(url; headers = headers, retry = false, status_exception = false)
            _is_transient_http_status(resp.status) && throw(_RetryableHTTPStatus(resp))
            return resp
        end()
    catch err
        err isa _RetryableHTTPStatus && return err.response
        rethrow()
    end
end

strip_ensembl_version(id::AbstractString)::String = String(split(String(id), '.'; limit = 2)[1])

function _strip_numeric_prefix(id::AbstractString)
    parts = split(String(id), '|'; limit = 2)
    return length(parts) == 2 && all(isdigit, parts[1]) ? String(parts[2]) : String(id)
end

function sequence_name_variants(id::AbstractString)
    variants = String[]
    for base in (String(id), _strip_numeric_prefix(id))
        push!(variants, base)
        core = strip_ensembl_version(base)
        core == base || push!(variants, core)
    end
    return unique(variants)
end

function resolve_sequence_name(msa::AbstractMultipleSequenceAlignment,
        ids;
        fallback::Bool = false)
    names = String.(sequencenames(msa))
    lookup = Dict{String, String}()
    # Store both full and version-stripped names for Ensembl IDs.
    for name in names
        for variant in sequence_name_variants(name)
            haskey(lookup, variant) || (lookup[variant] = name)
        end
    end

    iterable_ids = ids isa AbstractString ? (ids,) : ids
    for id in iterable_ids
        for variant in sequence_name_variants(id)
            haskey(lookup, variant) && return lookup[variant]
        end
    end
    return fallback && !isempty(names) ? first(names) : nothing
end

function protein_alignment_stats(query_seq::AbstractString,
        reference_seq::AbstractString;
        include_alignment::Bool = false)
    result = pairalign(GlobalAlignment(), LongAA(uppercase(String(reference_seq))),
        LongAA(uppercase(String(query_seq))), _PROTEIN_ALIGNMENT_SCORE_MODEL)
    aln = alignment(result)
    mismatches = count_mismatches(aln)
    insertions = count_insertions(aln)
    deletions = count_deletions(aln)
    stats = (;
        identical = mismatches == 0 && insertions == 0 && deletions == 0,
        mismatches,
        insertions,
        deletions
    )
    return include_alignment ? merge((; aln), stats) : stats
end

is_ensembl_transcript_id(id::AbstractString)::Bool = occursin(r"^ENS[A-Z]*T[0-9]+(\.[0-9]+)?$", String(id))

# UniProt documents accessions as either [OPQ][0-9][A-Z0-9]{3}[0-9] or
# [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}. The optional suffix accepts
# UniProt isoform IDs, which are written as accession-number plus "-<number>".
is_uniprot_id(id::AbstractString)::Bool = occursin(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9])(-[0-9]+)?$",
    String(id))

function id_kind(id::AbstractString)::Symbol
    is_ensembl_transcript_id(id) && return :ensembl_transcript
    is_uniprot_id(id) && return :uniprot
    error("Input ID $(id) is not recognized as a UniProt accession or Ensembl transcript ID.")
end

format_pid(pid::Real) = @sprintf("%.1f", Float64(pid))

function decode_body(resp::HTTP.Response)::String
    body = resp.body
    if length(body) >= 2 && body[1] == 0x1f && body[2] == 0x8b
        return String(transcode(GzipDecompressor, body))
    end
    return String(body)
end

function prepare_output_dir(input_id::AbstractString;
        workdir::Union{Nothing, AbstractString} = nothing,
        output_dir::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false)
    if workdir !== nothing && output_dir !== nothing &&
       abspath(workdir) != abspath(output_dir)
        error("Use either workdir or output_dir, or pass the same path for both.")
    end
    root = workdir === nothing ? output_dir : workdir
    root = root === nothing ? abspath(String(input_id)) : abspath(String(root))

    if isfile(root)
        error("Output path $(root) exists and is a file.")
    end
    # The workdir itself is never deleted. When overwrite=true, each pipeline
    # stage removes only the package-owned subdirectories it is about to rebuild.
    mkpath(root)
    return root
end

function safe_rm(path::AbstractString, root::AbstractString)
    abs_path = abspath(path)
    abs_root = abspath(root)
    rel = relpath(abs_path, abs_root)
    # Deletions are allowed only inside the active work directory.
    if rel == "." || startswith(rel, "..") || isabspath(rel)
        error("Refusing to remove $(abs_path) because it is outside workdir $(abs_root).")
    end
    rm(abs_path; recursive = true, force = true)
    return nothing
end

function run_logged(cmd::Cmd;
        stdout_path::AbstractString,
        stderr_path::AbstractString,
        workdir::Union{Nothing, AbstractString} = nothing)
    mkpath(dirname(stdout_path))
    mkpath(dirname(stderr_path))
    open(stdout_path, "w") do out_io
        open(stderr_path, "w") do err_io
            if workdir === nothing
                run(pipeline(cmd; stdout = out_io, stderr = err_io))
            else
                cd(workdir) do
                    run(pipeline(cmd; stdout = out_io, stderr = err_io))
                end
            end
        end
    end
    return nothing
end

function fasta_sequence(content::AbstractString)::Union{Nothing, String}
    seq = String[]
    for line in split(content, '\n')
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, '>') && continue
        push!(seq, stripped)
    end
    isempty(seq) && return nothing
    return uppercase(join(seq))
end

function _wrap_sequence(seq::AbstractString; width::Int = 60)
    io = IOBuffer()
    i = firstindex(seq)
    while i <= lastindex(seq)
        j = min(i + width - 1, lastindex(seq))
        println(io, seq[i:j])
        i = j + 1
    end
    return String(take!(io))
end

format_fasta(id::AbstractString,
    seq::AbstractString)::String = string(">", id, "\n", _wrap_sequence(uppercase(String(seq))))

function write_fasta(path::AbstractString, records)
    mkpath(dirname(path))
    open(path, "w") do io
        for (name, seq) in records
            println(io, '>', name)
            println(io, _wrap_sequence(uppercase(String(seq))))
        end
    end
    return path
end

function write_text(path::AbstractString, text::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        print(io, text)
    end
    return path
end

function write_json(path::AbstractString, obj)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.pretty(io, obj)
        println(io)
    end
    return path
end

function ensure_mmseqs_db(db::AbstractString)
    for db_path in (String(db), string(db, "_aln"), string(db, "_seq"))
        dbtype = string(db_path, ".dbtype")
        isfile(dbtype) ||
            error("MMseqs2 database at $(db_path) looks incomplete; missing $(dbtype).")
    end
    return String(db)
end

function result_summary(result::IdunaResult)
    return (;
        input_id = result.input_id,
        workdir = result.workdir,
        status = String(result.status),
        warnings = result.warnings,
        target = (;
            input_kind = String(result.target.input_kind),
            uniprot_id = result.target.uniprot_id,
            ensembl_gene_id = result.target.ensembl_gene_id,
            transcript_id = result.target.transcript_id,
            ensembl_protein_id = result.target.ensembl_protein_id,
            species = result.target.species,
            sequence_validated = result.target.sequence_validated,
            mapping_confirmed = result.target.mapping_confirmed
        ),
        thoraxe_msa = (;
            baseline_fasta = result.thoraxe_msa.baseline_fasta,
            baseline_stockholm = result.thoraxe_msa.baseline_stockholm,
            pid_summary = result.thoraxe_msa.pid_summary,
            best_pid = result.thoraxe_msa.best_seed.pid,
            best_seed_stockholm = result.thoraxe_msa.best_seed.stockholm_path,
            warnings = result.thoraxe_msa.warnings,
            status = String(result.thoraxe_msa.status)
        ),
        expansion = result.expansion === nothing ? nothing :
                    (;
            match_stockholm = result.expansion.match_stockholm,
            full_stockholm = result.expansion.full_stockholm,
            a3m_path = result.expansion.a3m_path,
            hits_fasta = result.expansion.hits_fasta,
            n_hits = result.expansion.n_hits,
            n_new_hits = result.expansion.n_new_hits,
            status = String(result.expansion.status)
        ),
        validation = (;
            stats_path = result.validation.stats_path,
            seed_nseq = result.validation.seed_nseq,
            expanded_nseq = result.validation.expanded_nseq,
            aln_identical = result.validation.aln_identical,
            status = String(result.validation.status)
        )
    )
end

end
