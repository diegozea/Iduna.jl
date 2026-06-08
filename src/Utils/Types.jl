using AutoPrettyPrinting: @def_pprint

"""
    ResolvedTarget

Resolved IDs and sequence checks for one input protein or transcript.

# Fields

- `input_id::String`: ID given by the user.
- `input_kind::Symbol`: Kind of input ID, such as `:uniprot` or
  `:ensembl_transcript`.
- `uniprot_id::Union{Nothing, String} = nothing`: UniProt accession. `nothing`
  means no UniProt accession was supplied or resolved.
- `ensembl_gene_id::String`: Ensembl gene ID used by ThorAxe.
- `transcript_id::String`: Ensembl transcript ID used by ThorAxe.
- `ensembl_protein_id::Union{Nothing, String} = nothing`: Ensembl protein ID.
  `nothing` means no Ensembl protein ID was supplied or resolved.
- `species::Union{Nothing, String} = nothing`: Species name. `nothing` means it
  could not be resolved.
- `uniprot_sequence_path::Union{Nothing, String} = nothing`: Path to the saved
  UniProt protein sequence. `nothing` means no UniProt sequence was saved.
- `ensembl_protein_sequence_path::Union{Nothing, String} = nothing`: Path to the
  saved Ensembl protein sequence. `nothing` means no Ensembl sequence was saved.
- `sequence_validated::Union{Nothing, Bool} = nothing`: Whether the UniProt and
  Ensembl protein sequences matched. `nothing` means the comparison was not run.
- `mapping_confirmed::Union{Nothing, Bool} = nothing`: Whether the
  UniProt-to-Ensembl mapping was confirmed. `nothing` means it was not checked.
- `workdir::Union{Nothing, String} = nothing`: Work directory used for paths
  stored in this result. `nothing` means paths are not tied to a work directory.
- `warnings::Vector{String} = String[]`: Non-fatal problems found while resolving
  the target. An empty vector means no warnings were recorded.
"""
Base.@kwdef struct ResolvedTarget
    input_id::String
    input_kind::Symbol
    uniprot_id::Union{Nothing, String} = nothing
    ensembl_gene_id::String
    transcript_id::String
    ensembl_protein_id::Union{Nothing, String} = nothing
    species::Union{Nothing, String} = nothing
    uniprot_sequence_path::Union{Nothing, String} = nothing
    ensembl_protein_sequence_path::Union{Nothing, String} = nothing
    sequence_validated::Union{Nothing, Bool} = nothing
    mapping_confirmed::Union{Nothing, Bool} = nothing
    workdir::Union{Nothing, String} = nothing
    warnings::Vector{String} = String[]
end

"""
    SeedSelection

A ThorAxe percent identity (PID) seed selected for later steps.

# Fields

- `pid::Float64`: Percent identity threshold used to build this seed.
- `median_identity::Union{Missing, Float64}`: Median identity score used during
  seed selection. `missing` means the score was unavailable.
- `mean_identity::Union{Missing, Float64}`: Mean identity score used during seed
  selection. `missing` means the score was unavailable.
- `stockholm_path::String`: Path to the selected seed MSA in Stockholm format.
- `fasta_path::Union{Nothing, String} = nothing`: Path to the selected seed MSA
  in FASTA format. `nothing` means no FASTA copy is available.
- `s_exon_blocks_tsv::Union{Nothing, String} = nothing`: Path to the table that
  maps MSA columns to s-exons. `nothing` means no table is available.
- `summary_path::String`: Path to the seed-selection summary table.
- `used_fallback_dir::Bool = false`: Whether a fallback ThorAxe output directory
  was used.
- `workdir::Union{Nothing, String} = nothing`: Work directory used for paths
  stored in this result. `nothing` means paths are not tied to a work directory.
"""
Base.@kwdef struct SeedSelection
    pid::Float64
    median_identity::Union{Missing, Float64}
    mean_identity::Union{Missing, Float64}
    stockholm_path::String
    fasta_path::Union{Nothing, String} = nothing
    s_exon_blocks_tsv::Union{Nothing, String} = nothing
    summary_path::String
    used_fallback_dir::Bool = false
    workdir::Union{Nothing, String} = nothing
end

"""
    ThorAxeMSAResult

Paths and metadata from the ThorAxe MSA-building stage.

# Fields

- `input_dir::String`: Directory with the transcript-query input bundle.
- `thoraxe_dirs::Vector{String}`: ThorAxe run directories used to build candidates.
- `msa_dir::String`: Directory with Iduna's ThorAxe MSA outputs.
- `baseline_fastas::Vector{String}`: Full candidate MSAs in FASTA format.
- `baseline_stockholms::Vector{String}`: Full candidate MSAs in Stockholm format.
- `sequence_fastas::Vector{String}`: Protein sequence files used for candidates.
- `species_files::Vector{String}`: Species list files used for candidates.
- `pid_summary::String`: CSV table describing percent identity (PID) candidates
  and selected seeds.
- `seeds::Vector{SeedSelection}`: Seed MSAs selected for validation and expansion.
- `logs_dir::String`: Directory with ThorAxe logs.
- `pid_sample_count::Int = 0`: Number of species samples used to score each PID.
  `0` means sampling was disabled or no sampling count was recorded.
- `pid_sample_fraction::Float64 = 1.0`: Fraction of non-reference species in each
  sample. `1.0` means all available non-reference species were kept.
- `pid_sample_seed::UInt64 = UInt64(0)`: Random seed used for PID sampling. `0`
  means no random seed was recorded.
- `sampling_strategy::Symbol = :independent`: How species samples were shared
  across PID values. Accepted values are `:common`, `:independent`, and `:input`;
  this constructor default mainly supports manually built or legacy results.
- `warnings::Vector{String} = String[]`: Non-fatal problems from the ThorAxe MSA
  stage. An empty vector means no warnings were recorded.
- `status::Symbol = :ok`: Stage status, such as `:ok`, `:warn`, or `:error`.
"""
Base.@kwdef struct ThorAxeMSAResult
    input_dir::String
    thoraxe_dirs::Vector{String}
    msa_dir::String
    baseline_fastas::Vector{String}
    baseline_stockholms::Vector{String}
    sequence_fastas::Vector{String}
    species_files::Vector{String}
    pid_summary::String
    seeds::Vector{SeedSelection}
    logs_dir::String
    pid_sample_count::Int = 0
    pid_sample_fraction::Float64 = 1.0
    pid_sample_seed::UInt64 = UInt64(0)
    sampling_strategy::Symbol = :independent
    warnings::Vector{String} = String[]
    status::Symbol = :ok
end

"""
    ExpansionResult

Paths and hit counts from the MMseqs2 and HMMER expansion stage.

# Fields

- `run_dir::String`: Directory for this expansion run.
- `seed_stockholm::String`: Seed MSA used for expansion.
- `seed_fasta::Union{Nothing, String} = nothing`: Seed MSA in FASTA format.
  `nothing` means no FASTA copy is available.
- `hits_fasta::String`: Raw sequence hits found by MMseqs2.
- `full_stockholm::String`: Expanded MSA with all columns.
- `match_stockholm::String`: Expanded MSA restricted to match columns.
- `a3m_path::String`: Expanded MSA in A3M format.
- `s_exon_blocks_tsv::Union{Nothing, String} = nothing`: Table that maps
  expanded MSA columns to s-exons. `nothing` means no table is available.
- `db_dir::String`: Directory with temporary MMseqs2 databases.
- `hmm_dir::String`: Directory with HMMER files.
- `logs_dir::String`: Directory with expansion logs.
- `n_hits::Int = 0`: Number of sequence hits found.
- `n_new_hits::Int = 0`: Number of hits not already present in the seed.
- `status::Symbol = :ok`: Stage status, such as `:ok`, `:warn`, or `:error`.
- `workdir::Union{Nothing, String} = nothing`: Work directory used for paths
  stored in this result. `nothing` means paths are not tied to a work directory.
"""
Base.@kwdef struct ExpansionResult
    run_dir::String
    seed_stockholm::String
    seed_fasta::Union{Nothing, String} = nothing
    hits_fasta::String
    full_stockholm::String
    match_stockholm::String
    a3m_path::String
    s_exon_blocks_tsv::Union{Nothing, String} = nothing
    db_dir::String
    hmm_dir::String
    logs_dir::String
    n_hits::Int = 0
    n_new_hits::Int = 0
    status::Symbol = :ok
    workdir::Union{Nothing, String} = nothing
end

"""
    ValidationResult

Statistics for the seed MSA, the expanded MSA, and the UniProt comparison.

# Fields

- `stats_path::String`: CSV file with the validation statistics.
- `query_name::Union{Nothing, String} = nothing`: Sequence name used for the
  query in the MSA. `nothing` means no query sequence was found.
- `query_vs_uniprot_path::Union{Nothing, String} = nothing`: Alignment report
  comparing the query with UniProt. `nothing` means no report is available.
- `seed_nseq::Union{Nothing, Int} = nothing`: Number of sequences in the seed MSA.
- `seed_ncol::Union{Nothing, Int} = nothing`: Number of columns in the seed MSA.
- `seed_clusters62::Union{Nothing, Int} = nothing`: Number of Hobohm clusters at
  62% identity in the seed MSA.
- `seed_neff80::Union{Nothing, Float64} = nothing`: Effective sequence count at
  80% identity in the seed MSA.
- `expanded_nseq::Union{Nothing, Int} = nothing`: Number of sequences in the
  expanded MSA. `nothing` means no expanded MSA statistics are available.
- `expanded_ncol::Union{Nothing, Int} = nothing`: Number of columns in the
  expanded MSA.
- `expanded_clusters62::Union{Nothing, Int} = nothing`: Number of Hobohm clusters
  at 62% identity in the expanded MSA.
- `expanded_neff80::Union{Nothing, Float64} = nothing`: Effective sequence count
  at 80% identity in the expanded MSA.
- `aln_identical::Union{Nothing, Bool} = nothing`: Whether the final query
  sequence matched the UniProt sequence. `nothing` means it was not checked.
- `aln_mismatches::Union{Nothing, Int} = nothing`: Number of mismatched residues
  in the UniProt comparison.
- `aln_insertions::Union{Nothing, Int} = nothing`: Number of insertions in the
  UniProt comparison.
- `aln_deletions::Union{Nothing, Int} = nothing`: Number of deletions in the
  UniProt comparison.
- `warnings::Vector{String} = String[]`: Non-fatal validation problems. An empty
  vector means no warnings were recorded.
- `status::Symbol = :ok`: Validation status, such as `:ok`, `:warn`, or `:error`.
"""
Base.@kwdef struct ValidationResult
    stats_path::String
    query_name::Union{Nothing, String} = nothing
    query_vs_uniprot_path::Union{Nothing, String} = nothing
    seed_nseq::Union{Nothing, Int} = nothing
    seed_ncol::Union{Nothing, Int} = nothing
    seed_clusters62::Union{Nothing, Int} = nothing
    seed_neff80::Union{Nothing, Float64} = nothing
    expanded_nseq::Union{Nothing, Int} = nothing
    expanded_ncol::Union{Nothing, Int} = nothing
    expanded_clusters62::Union{Nothing, Int} = nothing
    expanded_neff80::Union{Nothing, Float64} = nothing
    aln_identical::Union{Nothing, Bool} = nothing
    aln_mismatches::Union{Nothing, Int} = nothing
    aln_insertions::Union{Nothing, Int} = nothing
    aln_deletions::Union{Nothing, Int} = nothing
    warnings::Vector{String} = String[]
    status::Symbol = :ok
end

"""
    IdunaResult

Top-level result returned by [`Iduna.iduna`](@ref).

# Fields

- `input_id::String`: ID given by the user.
- `workdir::String`: Absolute path to the work directory.
- `target::ResolvedTarget`: Resolved target IDs and sequence checks.
- `thoraxe_msa::ThorAxeMSAResult`: ThorAxe MSA outputs and selected seeds.
- `expansions::Vector{Union{Missing, ExpansionResult}}`: Expansion results,
  indexed like `thoraxe_msa.seeds`.
- `validations::Vector{ValidationResult}`: Validation results for each seed.
- `stages::Vector{Any} = Any[]`: Compact summaries of completed or failed
  pipeline stages.
- `warnings::Vector{String} = String[]`: Non-fatal problems from the full run.
  An empty vector means no warnings were recorded.
- `status::Symbol = :ok`: Overall status, such as `:ok`, `:warn`, or `:error`.

`expansions` is empty when `no_expansion=true`. Otherwise it may contain
`missing` for a seed where expansion did not finish.
"""
Base.@kwdef struct IdunaResult
    input_id::String
    workdir::String
    target::ResolvedTarget
    thoraxe_msa::ThorAxeMSAResult
    expansions::Vector{Union{Missing, ExpansionResult}}
    validations::Vector{ValidationResult}
    stages::Vector{Any} = Any[]
    warnings::Vector{String} = String[]
    status::Symbol = :ok
end

struct _IdunaSummarySegment
    text::String
    color::Union{Nothing, Symbol}
end

_iduna_segment(text; color = nothing) = _IdunaSummarySegment(string(text), color)

function _iduna_printstyled(io::IO, text::AbstractString; color = nothing, bold = false)
    if !get(io, :color, false) || (color === nothing && !bold)
        print(io, text)
    elseif color === nothing
        printstyled(io, text; bold)
    else
        printstyled(io, text; color, bold)
    end
    return nothing
end

function _iduna_print_segments(io::IO, segments)
    for segment in segments
        _iduna_printstyled(io, segment.text; color = segment.color)
    end
    return nothing
end

_iduna_status_word(status) = lowercase(string(status))
_iduna_status_text(status::Symbol) = string(':', String(status))

function _iduna_semantic_color(value)
    word = _iduna_status_word(value)
    if word in ("ok", "done", "completed")
        return :green
    elseif word in ("warn", "warning", "warnings")
        return :yellow
    elseif word in ("error", "failed")
        return :red
    end
    return :light_black
end

function _iduna_status_segment(status)
    return _iduna_segment(_iduna_status_text(status); color = _iduna_semantic_color(status))
end

function _iduna_count_segment(count::Integer, singular::AbstractString;
        plural::AbstractString = string(singular, "s"), color = nothing, zero_color = nothing)
    label = count == 1 ? singular : plural
    segment_color = count == 0 ? zero_color : color
    return _iduna_segment("$(count) $(label)"; color = segment_color)
end

function _iduna_push_summary_part!(segments, segment)
    push!(segments, _iduna_segment(", "))
    push!(segments, segment)
    return segments
end

function _iduna_unknown_if_empty(value)
    text = string(value)
    isempty(text) && return _iduna_segment("unknown"; color = :light_black)
    return _iduna_segment(text)
end

function _iduna_target_summary(target::ResolvedTarget)
    return [
        _iduna_segment("ResolvedTarget"),
        _iduna_segment(", "),
        _iduna_unknown_if_empty(target.ensembl_gene_id),
        _iduna_segment(", "),
        _iduna_unknown_if_empty(target.transcript_id)
    ]
end

function _iduna_dimension_label(count::Integer, singular::AbstractString,
        plural::AbstractString)
    return count == 1 ? singular : plural
end

function _iduna_msa_dimensions(nseq::Union{Nothing, Integer},
        ncol::Union{Nothing, Integer}; compact::Bool = false)
    (nseq === nothing || ncol === nothing) && return nothing
    if compact
        return "($(nseq), $(ncol))"
    end
    seq_label = _iduna_dimension_label(nseq, "seq", "seqs")
    col_label = _iduna_dimension_label(ncol, "col", "cols")
    return "($(nseq) $(seq_label), $(ncol) $(col_label))"
end

function _iduna_validation_at(validations::AbstractVector, index::Integer)
    1 <= index <= length(validations) ? validations[index] : nothing
end

function _iduna_seed_dimensions(validations::AbstractVector, index::Integer;
        compact::Bool = false)
    validation = _iduna_validation_at(validations, index)
    validation === nothing && return nothing
    return _iduna_msa_dimensions(validation.seed_nseq, validation.seed_ncol; compact)
end

function _iduna_expanded_dimensions(expansions::AbstractVector,
        validations::AbstractVector, index::Integer; compact::Bool = false)
    1 <= index <= length(expansions) || return nothing
    ismissing(expansions[index]) && return nothing
    validation = _iduna_validation_at(validations, index)
    validation === nothing && return nothing
    return _iduna_msa_dimensions(validation.expanded_nseq, validation.expanded_ncol;
        compact)
end

function _iduna_seed_pid_dimensions(seeds::AbstractVector{SeedSelection},
        validations::AbstractVector, index::Integer; compact::Bool = false)
    text = format_pid(seeds[index].pid)
    dimensions = _iduna_seed_dimensions(validations, index; compact)
    dimensions === nothing && return text
    return "$(text) $(dimensions)"
end

function _iduna_thoraxe_msa_summary(thoraxe_msa::ThorAxeMSAResult,
        validations::AbstractVector)
    n_seeds = length(thoraxe_msa.seeds)
    pid_label = n_seeds == 1 ? "selected PID " : "selected PIDs "
    segments = [
        _iduna_segment("ThorAxeMSAResult"),
        _iduna_segment(", "),
        _iduna_count_segment(n_seeds, "seed"),
        _iduna_segment(", "),
        _iduna_segment(pid_label)
    ]
    if isempty(thoraxe_msa.seeds)
        push!(segments, _iduna_segment("unknown"; color = :light_black))
    else
        compact = length(thoraxe_msa.seeds) > 1
        has_dimensions = any(
            index -> _iduna_seed_dimensions(validations, index;
                compact) !== nothing,
            eachindex(thoraxe_msa.seeds))
        compact && has_dimensions && push!(segments, _iduna_segment("(seqs, cols) "))
        push!(segments,
            _iduna_segment(join(
                (_iduna_seed_pid_dimensions(thoraxe_msa.seeds, validations, index;
                     compact)
                for index in eachindex(thoraxe_msa.seeds)), ", ")))
    end
    return segments
end

function _iduna_classify_status(status)
    word = _iduna_status_word(status)
    if word in ("ok", "done", "completed")
        return :ok
    elseif word in ("warn", "warning", "warnings")
        return :warn
    elseif word in ("error", "failed")
        return :failed
    end
    return :unknown
end

function _iduna_expansion_status(expansion)
    ismissing(expansion) && return :missing
    _iduna_classify_status(expansion.status)
end

function _iduna_expansion_hits(expansion)
    ismissing(expansion) && return 0
    return expansion.n_hits
end

function _iduna_expansion_dimension_entries(seeds::AbstractVector,
        expansions::AbstractVector, validations::AbstractVector)
    entries = String[]
    for index in eachindex(expansions)
        dimensions = _iduna_expanded_dimensions(expansions, validations, index;
            compact = true)
        dimensions === nothing && continue
        label = 1 <= index <= length(seeds) ? format_pid(seeds[index].pid) : string(index)
        push!(entries, "$(label) $(dimensions)")
    end
    return entries
end

function _iduna_append_expansion_dimensions!(segments, expansions::AbstractVector,
        seeds::AbstractVector, validations::AbstractVector)
    n_expansions = length(expansions)
    if n_expansions == 1
        dimensions = _iduna_expanded_dimensions(expansions, validations, 1)
        dimensions === nothing ||
            _iduna_push_summary_part!(segments, _iduna_segment("selected MSA $(dimensions)"))
    elseif n_expansions > 1
        entries = _iduna_expansion_dimension_entries(seeds, expansions, validations)
        isempty(entries) ||
            _iduna_push_summary_part!(segments,
                _iduna_segment("selected MSAs (seqs, cols) $(join(entries, ", "))"))
    end
    return segments
end

function _iduna_expansions_summary(expansions::AbstractVector,
        seeds::AbstractVector, validations::AbstractVector)
    n_expansions = length(expansions)
    completed = count(expansion -> _iduna_expansion_status(expansion) === :ok, expansions)
    warned = count(expansion -> _iduna_expansion_status(expansion) === :warn, expansions)
    failed = count(expansion -> _iduna_expansion_status(expansion) === :failed, expansions)
    unknown = count(expansion -> _iduna_expansion_status(expansion) === :unknown, expansions)
    missing = count(expansion -> _iduna_expansion_status(expansion) === :missing, expansions)
    n_hits = sum(_iduna_expansion_hits, expansions; init = 0)

    segments = [_iduna_count_segment(n_expansions, "expansion")]
    if n_expansions == 0
        _iduna_push_summary_part!(segments, _iduna_segment("empty"; color = :light_black))
    else
        completed > 0 &&
            _iduna_push_summary_part!(segments,
                _iduna_count_segment(completed, "completed", plural = "completed",
                    color = :green))
        warned > 0 &&
            _iduna_push_summary_part!(segments,
                _iduna_count_segment(warned, "warn", plural = "warn", color = :yellow))
        failed > 0 &&
            _iduna_push_summary_part!(segments,
                _iduna_count_segment(failed, "failed", plural = "failed", color = :red))
        unknown > 0 &&
            _iduna_push_summary_part!(segments,
                _iduna_count_segment(unknown, "unknown", plural = "unknown",
                    color = :light_black))
        missing > 0 &&
            _iduna_push_summary_part!(segments,
                _iduna_count_segment(missing, "missing", plural = "missing",
                    color = :light_black))
    end
    _iduna_push_summary_part!(segments,
        _iduna_count_segment(n_hits, "hit"; color = nothing, zero_color = :light_black))
    _iduna_append_expansion_dimensions!(segments, expansions, seeds, validations)
    return segments
end

function _iduna_push_status_count!(segments, count::Integer, word::AbstractString,
        color::Symbol)
    count == 0 && return segments
    text = count == 1 ? word : "$(count) $(word)"
    return _iduna_push_summary_part!(segments, _iduna_segment(text; color))
end

function _iduna_validations_summary(validations::AbstractVector)
    ok = count(validation -> _iduna_classify_status(validation.status) === :ok, validations)
    warned = count(
        validation -> _iduna_classify_status(validation.status) === :warn, validations)
    failed = count(
        validation -> _iduna_classify_status(validation.status) === :failed, validations)
    unknown = length(validations) - ok - warned - failed

    segments = [_iduna_count_segment(length(validations), "result")]
    if isempty(validations)
        _iduna_push_summary_part!(segments, _iduna_segment("empty"; color = :light_black))
    else
        _iduna_push_status_count!(segments, ok, "ok", :green)
        _iduna_push_status_count!(segments, warned, "warn", :yellow)
        _iduna_push_status_count!(segments, failed, "failed", :red)
        _iduna_push_status_count!(segments, unknown, "unknown", :light_black)
    end
    return segments
end

function _iduna_stage_status(stage)
    if stage isa AbstractDict
        return get(stage, "status", nothing)
    elseif hasproperty(stage, :status)
        return getproperty(stage, :status)
    end
    return nothing
end

function _iduna_stages_summary(stages::AbstractVector)
    done = 0
    warned = 0
    failed = 0
    unknown = 0
    for stage in stages
        classification = _iduna_classify_status(_iduna_stage_status(stage))
        if classification === :ok
            done += 1
        elseif classification === :warn
            warned += 1
        elseif classification === :failed
            failed += 1
        else
            unknown += 1
        end
    end

    segments = [_iduna_count_segment(length(stages), "stage")]
    if isempty(stages)
        _iduna_push_summary_part!(segments, _iduna_segment("empty"; color = :light_black))
    else
        done > 0 &&
            _iduna_push_summary_part!(segments,
                _iduna_count_segment(done, "done", plural = "done", color = :green))
        warned > 0 &&
            _iduna_push_summary_part!(segments,
                _iduna_count_segment(warned, "warn", plural = "warn", color = :yellow))
        failed > 0 &&
            _iduna_push_summary_part!(segments,
                _iduna_count_segment(failed, "failed", plural = "failed", color = :red))
        unknown > 0 &&
            _iduna_push_summary_part!(segments,
                _iduna_count_segment(unknown, "unknown", plural = "unknown",
                    color = :light_black))
    end
    return segments
end

function _iduna_result_rows(result::IdunaResult)
    return [
        ("input_id", [_iduna_segment(result.input_id)]),
        ("workdir", [_iduna_segment(result.workdir)]),
        ("target", _iduna_target_summary(result.target)),
        ("thoraxe_msa", _iduna_thoraxe_msa_summary(result.thoraxe_msa,
            result.validations)),
        ("expansions",
            _iduna_expansions_summary(result.expansions,
                result.thoraxe_msa.seeds, result.validations)),
        ("validations", _iduna_validations_summary(result.validations)),
        ("stages", _iduna_stages_summary(result.stages)),
        ("warnings",
            [_iduna_count_segment(length(result.warnings), "warning";
                color = :yellow, zero_color = :light_black)]),
        ("status", [_iduna_status_segment(result.status)])
    ]
end

function Base.show(io::IO, ::MIME"text/plain", result::IdunaResult)
    print(io, "IdunaResult ", result.input_id, " [")
    _iduna_print_segments(io, [_iduna_status_segment(result.status)])
    println(io, "]")
    println(io)

    rows = _iduna_result_rows(result)
    field_width = maximum(length(label) for (label, _) in rows)
    _iduna_printstyled(io, rpad("FIELD", field_width); bold = true)
    print(io, "  ")
    _iduna_printstyled(io, "SUMMARY"; bold = true)
    for (label, segments) in rows
        println(io)
        print(io, rpad(label, field_width), "  ")
        _iduna_print_segments(io, segments)
    end
    return nothing
end

@def_pprint mime_types="text/plain" base_show=true ResolvedTarget
@def_pprint mime_types="text/plain" base_show=true SeedSelection
@def_pprint mime_types="text/plain" base_show=true ThorAxeMSAResult
@def_pprint mime_types="text/plain" base_show=true ExpansionResult
@def_pprint mime_types="text/plain" base_show=true ValidationResult
