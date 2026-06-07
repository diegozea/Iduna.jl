using AutoPrettyPrinting: @def_pprint

"""
    ResolvedTarget

Resolved identifiers and sequence-validation metadata for one input ID.
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

The ThorAxe PID seed chosen for expansion, plus the summary values used to pick it.
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

Paths and metadata produced by the ThorAxe MSA-building stage.
`sampling_strategy` records how species samples were prepared for PID
selection.
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

Paths and hit counts produced by the MMseqs2/HMMER expansion stage.
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

Seed and expanded-MSA statistics plus optional query-vs-UniProt checks.
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

Top-level result returned by `Iduna.iduna`. It keeps the full pipeline status
and links to the per-stage result objects. `expansions` is empty when the run
was requested with `no_expansion=true`; otherwise it is indexed like
`thoraxe_msa.seeds` and can contain `nothing` for missing expansion results.
"""
Base.@kwdef struct IdunaResult
    input_id::String
    workdir::String
    target::ResolvedTarget
    thoraxe_msa::ThorAxeMSAResult
    expansions::Vector{Union{Nothing, ExpansionResult}}
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

function _iduna_thoraxe_msa_summary(thoraxe_msa::ThorAxeMSAResult)
    segments = [
        _iduna_segment("ThorAxeMSAResult"),
        _iduna_segment(", "),
        _iduna_count_segment(length(thoraxe_msa.seeds), "seed"),
        _iduna_segment(", selected PIDs ")
    ]
    if isempty(thoraxe_msa.seeds)
        push!(segments, _iduna_segment("unknown"; color = :light_black))
    else
        push!(segments,
            _iduna_segment(join((format_pid(seed.pid) for seed in thoraxe_msa.seeds),
                ", ")))
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
    expansion === nothing && return :missing
    _iduna_classify_status(expansion.status)
end

function _iduna_expansion_hits(expansion)
    expansion === nothing && return 0
    return expansion.n_hits
end

function _iduna_expansions_summary(expansions::AbstractVector)
    n_slots = length(expansions)
    completed = count(expansion -> _iduna_expansion_status(expansion) === :ok, expansions)
    warned = count(expansion -> _iduna_expansion_status(expansion) === :warn, expansions)
    failed = count(expansion -> _iduna_expansion_status(expansion) === :failed, expansions)
    unknown = count(expansion -> _iduna_expansion_status(expansion) === :unknown, expansions)
    missing = count(expansion -> _iduna_expansion_status(expansion) === :missing, expansions)
    n_hits = sum(_iduna_expansion_hits, expansions; init = 0)

    segments = [_iduna_count_segment(n_slots, "slot")]
    if n_slots == 0
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
                _iduna_count_segment(missing, "nothing", plural = "nothing",
                    color = :light_black))
    end
    _iduna_push_summary_part!(segments,
        _iduna_count_segment(n_hits, "hit"; color = nothing, zero_color = :light_black))
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
        ("thoraxe_msa", _iduna_thoraxe_msa_summary(result.thoraxe_msa)),
        ("expansions", _iduna_expansions_summary(result.expansions)),
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
