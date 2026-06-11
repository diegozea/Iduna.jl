"""
    strip_ensembl_version(id) -> String

Remove the version suffix from an Ensembl ID.

# Arguments

- `id::AbstractString`: Ensembl ID, with or without a version suffix.
"""
strip_ensembl_version(id::AbstractString)::String = String(split(
    String(id), '.'; limit = 2)[1])

function _strip_numeric_prefix(id::AbstractString)
    parts = split(String(id), '|'; limit = 2)
    return length(parts) == 2 && all(isdigit, parts[1]) ? String(parts[2]) : String(id)
end

"""
    sequence_name_variants(id) -> Vector{String}

Return common sequence-name forms for an ID, including version-stripped Ensembl
names and names without numeric prefixes.

# Arguments

- `id::AbstractString`: sequence name or biological identifier to normalize.
"""
function sequence_name_variants(id::AbstractString)
    variants = String[]
    for base in (String(id), _strip_numeric_prefix(id))
        push!(variants, base)
        core = strip_ensembl_version(base)
        core == base || push!(variants, core)
    end
    return unique(variants)
end

"""
    resolve_sequence_name(msa, ids; fallback=false)

Find the first sequence name in `msa` that matches one of the requested IDs.

# Arguments

- `msa`: MSA whose sequence names should be searched.
- `ids`: one ID or an iterable of IDs to match against sequence-name variants.

# Keywords

- `fallback::Bool = false`: return the first sequence name when no ID matches.

# Returns

- The matching sequence name, or `nothing` when no match is found.
"""
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

"""
    protein_alignment_stats(query_seq, reference_seq; include_alignment=false)

Align two protein sequences and count mismatches, insertions, and deletions.

# Arguments

- `query_seq::AbstractString`: query protein sequence.
- `reference_seq::AbstractString`: reference protein sequence.

# Keywords

- `include_alignment::Bool = false`: include the BioAlignments alignment object
  in the returned value.

# Returns

- A named tuple with `identical`, `mismatches`, `insertions`, and `deletions`.
"""
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

# ID Classification
# -----------------

"""
    is_ensembl_transcript_id(id) -> Bool

Return `true` when `id` looks like an Ensembl transcript ID.

# Arguments

- `id::AbstractString`: identifier to classify.
"""
is_ensembl_transcript_id(id::AbstractString)::Bool = occursin(
    r"^ENS[A-Z]*T[0-9]+(\.[0-9]+)?$", String(id))

# UniProt documents accessions as either [OPQ][0-9][A-Z0-9]{3}[0-9] or
# [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}. The optional suffix accepts
# UniProt isoform IDs, which are written as accession-number plus "-<number>".
"""
    is_uniprot_id(id) -> Bool

Return `true` when `id` looks like a UniProt accession or isoform accession.

# Arguments

- `id::AbstractString`: identifier to classify.
"""
is_uniprot_id(id::AbstractString)::Bool = occursin(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9])(-[0-9]+)?$",
    String(id))

"""
    id_kind(id) -> Symbol

Classify an input ID as `:uniprot` or `:ensembl_transcript`.

# Arguments

- `id::AbstractString`: identifier to classify.

# Throws

- `ErrorException`: if the ID does not look like either supported ID type.
"""
function id_kind(id::AbstractString)::Symbol
    is_ensembl_transcript_id(id) && return :ensembl_transcript
    is_uniprot_id(id) && return :uniprot
    error("Input ID $(id) is not recognized as a UniProt accession or Ensembl transcript ID.")
end
