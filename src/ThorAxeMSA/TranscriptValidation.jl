# Transcript Validation
# ---------------------

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
        msa::AbstractMultipleSequenceAlignment;
        workdir::Union{Nothing, AbstractString} = target.workdir)
    target.uniprot_sequence_path === nothing && return String[]
    artifact_workdir = target.workdir === nothing ? workdir : target.workdir
    uniprot_sequence_path = artifact_workdir === nothing ? target.uniprot_sequence_path :
                            _resolve_artifact_path(
        target.uniprot_sequence_path, artifact_workdir)
    isfile(uniprot_sequence_path) || return String[
        "UniProt sequence file is missing; skipped ThorAxe transcript validation."]

    query_name,
    query_seq = _extract_reference_sequence(msa,
        target.ensembl_gene_id, target.transcript_id)
    reference_seq = _read_single_fasta_sequence(uniprot_sequence_path)
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

function _hhsuite_query_indices(seq::AbstractString, start::Integer, stop::Integer)
    return EPLI._hhsuite_query_indices(seq, start, stop)
end

function _parse_hhsuite_query_segment(line::AbstractString)
    return EPLI._parse_hhsuite_query_segment(line)
end

function _append_hhsuite_code_line!(positions, codes, line::AbstractString, cols, indices)
    return EPLI._append_hhsuite_code_line!(positions, codes, line, cols, indices)
end

function _get_codes(output::AbstractString)
    return EPLI._get_codes(output)
end

function _identity_from_codes(positions::Vector{Int}, codes::Vector{Char})
    counts = EPLI._identity_counts_from_codes(positions, codes)
    score = (;
        matched_positions = counts.matched_positions,
        comparable_positions = counts.comparable_positions)
    return EPLI.comparable_positions_normalization(score).normalized_score
end

"""
    compute_identity_against_reference(reference_fasta, sample_fasta;
                                       logs_dir=nothing, label="seed")

Compare a sampled MSA with its full reference MSA using HHsuite.

# Arguments

- `reference_fasta::AbstractString`: full reference MSA in FASTA format.
- `sample_fasta::AbstractString`: sampled MSA in FASTA format.

# Keywords

- `logs_dir = nothing`: directory where the HHalign output should be copied. When
  `nothing`, the temporary HHalign output is not copied.
- `label::AbstractString = "seed"`: label used in the copied log file name.

# Returns

- Percent identity as a `Float64`.
"""
function compute_identity_against_reference(reference_fasta::AbstractString,
        sample_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        label::AbstractString = "seed")
    score = EPLI.hhsuite_identity_score(reference_fasta, sample_fasta;
        logs_dir, label)
    return EPLI.comparable_positions_normalization(score).normalized_score
end
