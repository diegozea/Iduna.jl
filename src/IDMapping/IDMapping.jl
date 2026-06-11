"""
    IDMapping

Resolve UniProt and Ensembl IDs into the target information needed by Iduna.
"""
module IDMapping

import JSON

using ..Utils: ResolvedTarget, _http_get_request, _http_get_with_retries,
               _is_transient_http_status,
               decode_body, fasta_sequence, format_fasta, id_kind,
               is_ensembl_transcript_id, protein_alignment_stats,
               strip_ensembl_version, write_text

export EnsemblCandidate,
       fetch_uniprot_entry,
       resolve_target,
       resolve_transcript_gene,
       sequences_match

const _UNIPROT_BASE = "https://rest.uniprot.org"
const _ENSEMBL_REST_BASE = "https://rest.ensembl.org"
const _JSON_HEADERS = ["Accept" => "application/json"]
const _FASTA_HEADERS = ["Accept" => "text/x-fasta"]

# Result Types
# ------------

"""
    EnsemblCandidate

One possible Ensembl transcript match for a UniProt entry.

# Fields

- `transcript_id::String`: Ensembl transcript ID.
- `ensembl_gene_id::Union{Nothing, String} = nothing`: Ensembl gene ID for the
  transcript. `nothing` means no gene ID was supplied or resolved.
- `ensembl_protein_id::Union{Nothing, String} = nothing`: Ensembl protein ID for
  the transcript. `nothing` means no protein ID was supplied or resolved.
- `isoform_id::Union{Nothing, String} = nothing`: UniProt isoform ID linked to
  this transcript. `nothing` means no isoform ID was supplied or resolved.
- `species::Union{Nothing, String} = nothing`: Species name from UniProt or
  Ensembl. `nothing` means no species was supplied or resolved.
- `sequence_validated::Bool = false`: Whether the protein sequence matched
  UniProt. `false` means it did not match or was not checked.
- `mapping_confirmed::Union{Nothing, Bool} = nothing`: Whether the mapping was
  confirmed by sequence or metadata. `nothing` means it was not checked.
"""
Base.@kwdef struct EnsemblCandidate
    transcript_id::String
    ensembl_gene_id::Union{Nothing, String} = nothing
    ensembl_protein_id::Union{Nothing, String} = nothing
    isoform_id::Union{Nothing, String} = nothing
    species::Union{Nothing, String} = nothing
    sequence_validated::Bool = false
    mapping_confirmed::Union{Nothing, Bool} = nothing
end

Base.@kwdef struct _EnsemblTranscriptLookup
    transcript_id::String
    ensembl_gene_id::String
    species::Union{Nothing, String} = nothing
end

Base.@kwdef struct _UniProtEntry
    id::String
    species::Union{Nothing, String}
    gene_ids::Vector{String}
    transcript_ids::Vector{String}
    transcript_to_protein::Dict{String, String}
    transcript_to_gene::Dict{String, String}
    transcript_to_isoform::Dict{String, String}
    protein_sequence::Union{Nothing, String}
end

# HTTP Fetching
# -------------

function _http_get(url::AbstractString,
        headers = _JSON_HEADERS;
        retries::Int = 4,
        sleep_seconds::Real = 1.5,
        http_get::Function = _http_get_request)
    resp = _http_get_with_retries(url, headers; retries, sleep_seconds, http_get)
    if resp.status == 200
        return resp
    elseif _is_transient_http_status(resp.status)
        @warn "HTTP request failed after retries." url status=resp.status
    end
    return nothing
end

# UniProt Parsing
# ---------------

"""
    fetch_uniprot_entry(uniprot_id)

Fetch UniProt metadata and linked Ensembl transcript IDs for one accession.

# Arguments

- `uniprot_id::AbstractString`: UniProt accession to fetch.

# Returns

- An internal record with species, sequence, gene IDs, and transcript mappings.

# Throws

- `ErrorException`: if UniProt metadata cannot be fetched.
"""
function fetch_uniprot_entry(uniprot_id::AbstractString;
        _http_get_fn::Function = _http_get)::_UniProtEntry
    url = "$(_UNIPROT_BASE)/uniprotkb/$(strip(String(uniprot_id))).json"
    resp = _http_get_fn(url)
    resp === nothing && error("Could not fetch UniProt metadata for $(uniprot_id).")
    data = JSON.parse(decode_body(resp))
    gene_ids, transcripts, transcript_to_protein,
    transcript_to_gene, transcript_to_isoform = _parse_xrefs(data)
    return _UniProtEntry(;
        id = String(uniprot_id),
        species = _get_species(data),
        gene_ids,
        transcript_ids = transcripts,
        transcript_to_protein,
        transcript_to_gene,
        transcript_to_isoform,
        protein_sequence = _extract_uniprot_sequence(data)
    )
end

function _get_species(data)::Union{Nothing, String}
    organism = get(data, "organism", nothing)
    organism === nothing && return nothing
    scientific = get(organism, "scientificName", nothing)
    return scientific isa AbstractString ? String(scientific) : nothing
end

function _extract_uniprot_sequence(data)::Union{Nothing, String}
    seq = get(data, "sequence", nothing)
    seq === nothing && return nothing
    value = get(seq, "value", nothing)
    return value isa AbstractString ? uppercase(String(value)) : nothing
end

function _new_xref_data()
    return (;
        gene_ids = String[],
        transcripts = String[],
        transcript_to_protein = Dict{String, String}(),
        transcript_to_gene = Dict{String, String}(),
        transcript_to_isoform = Dict{String, String}()
    )
end

function _record_ensembl_transcript!(xrefs, ref)::Union{Nothing, String}
    transcript = get(ref, "id", nothing)
    transcript isa AbstractString || return nothing

    transcript = String(transcript)
    push!(xrefs.transcripts, transcript)
    isoform = get(ref, "isoformId", nothing)
    isoform isa AbstractString &&
        (xrefs.transcript_to_isoform[transcript] = String(isoform))
    return transcript
end

function _record_ensembl_property!(xrefs, transcript::Union{Nothing, String}, prop)
    key = get(prop, "key", "")
    value = get(prop, "value", nothing)
    if key == "GeneId" && value isa AbstractString
        gene_id = String(value)
        push!(xrefs.gene_ids, gene_id)
        transcript !== nothing && (xrefs.transcript_to_gene[transcript] = gene_id)
    elseif key == "ProteinId" && value isa AbstractString && transcript !== nothing
        xrefs.transcript_to_protein[transcript] = String(value)
    end
    return nothing
end

function _fill_default_transcript_genes!(transcript_to_gene, transcripts, gene_ids)
    isempty(gene_ids) && return transcript_to_gene

    default_gene = first(gene_ids)
    for transcript in transcripts
        # Some UniProt records give one gene ID but omit it from individual transcripts.
        haskey(transcript_to_gene, transcript) ||
            (transcript_to_gene[transcript] = default_gene)
    end
    return transcript_to_gene
end

function _parse_xrefs(data)
    xrefs = _new_xref_data()

    # UniProt stores Ensembl links as cross-references with extra properties.
    for ref in get(data, "uniProtKBCrossReferences", Any[])
        get(ref, "database", "") == "Ensembl" || continue
        transcript = _record_ensembl_transcript!(xrefs, ref)
        for prop in get(ref, "properties", Any[])
            _record_ensembl_property!(xrefs, transcript, prop)
        end
    end

    gene_ids = unique(xrefs.gene_ids)
    transcripts = unique(xrefs.transcripts)
    _fill_default_transcript_genes!(xrefs.transcript_to_gene, transcripts, gene_ids)
    return gene_ids, transcripts, xrefs.transcript_to_protein,
    xrefs.transcript_to_gene,
    xrefs.transcript_to_isoform
end

# Sequence Fetching and Validation
# --------------------------------

function _fetch_ensembl_protein_sequence(protein_id::AbstractString;
        _http_get_fn::Function = _http_get)::Union{Nothing, String}
    core = strip_ensembl_version(protein_id)
    url = "$(_ENSEMBL_REST_BASE)/sequence/id/$(core)?type=protein"
    resp = _http_get_fn(url, _FASTA_HEADERS)
    resp === nothing && return nothing
    return fasta_sequence(decode_body(resp))
end

function _fetch_uniprot_fasta_sequence(uniprot_id::AbstractString;
        _http_get_fn::Function = _http_get)::String
    url = "$(_UNIPROT_BASE)/uniprotkb/$(strip(String(uniprot_id))).fasta"
    resp = _http_get_fn(url, _FASTA_HEADERS; retries = 5, sleep_seconds = 2.0)
    resp === nothing && error("Could not fetch UniProt FASTA for $(uniprot_id).")
    seq = fasta_sequence(decode_body(resp))
    seq === nothing && error("UniProt FASTA for $(uniprot_id) did not contain a sequence.")
    return seq
end

"""
    sequences_match(uniprot_seq, ensembl_seq) -> Bool

Return `true` when two protein sequences are the same after FASTA parsing and
case normalization.

# Arguments

- `uniprot_seq::AbstractString`: UniProt protein sequence.
- `ensembl_seq::AbstractString`: Ensembl protein sequence.
"""
function sequences_match(uniprot_seq::AbstractString, ensembl_seq::AbstractString)::Bool
    protein_alignment_stats(ensembl_seq, uniprot_seq).identical
end

function _validate_candidates(entry::_UniProtEntry, sequence_dir::AbstractString;
        _uniprot_fasta_fetcher::Function = _fetch_uniprot_fasta_sequence,
        _ensembl_protein_fetcher::Function = _fetch_ensembl_protein_sequence)
    # Keep only Ensembl proteins that exactly match the UniProt sequence.
    uniprot_seq = entry.protein_sequence === nothing ?
                  _uniprot_fasta_fetcher(entry.id) : entry.protein_sequence
    uniprot_path = joinpath(sequence_dir, "uniprot", "$(entry.id).fasta")
    write_text(uniprot_path, format_fasta(entry.id, uniprot_seq))

    candidates = EnsemblCandidate[]
    for transcript in entry.transcript_ids
        protein_id = get(entry.transcript_to_protein, transcript, nothing)
        protein_id === nothing && continue
        protein_seq = _ensembl_protein_fetcher(protein_id)
        protein_seq === nothing && continue
        protein_path = joinpath(sequence_dir, "ensembl_proteins", "$(protein_id).fasta")
        write_text(protein_path, format_fasta(protein_id, protein_seq))
        if sequences_match(uniprot_seq, protein_seq)
            push!(candidates,
                EnsemblCandidate(;
                    transcript_id = transcript,
                    ensembl_gene_id = get(entry.transcript_to_gene, transcript, nothing),
                    ensembl_protein_id = protein_id,
                    isoform_id = get(entry.transcript_to_isoform, transcript, nothing),
                    species = entry.species,
                    sequence_validated = true,
                    # The candidate came from UniProt's own Ensembl cross-references.
                    mapping_confirmed = true
                ))
        end
    end
    return candidates, uniprot_path
end

# Transcript Metadata
# -------------------

function _choose_candidate(candidates::Vector{EnsemblCandidate}, transcript_id::Union{
        Nothing, AbstractString})
    isempty(candidates) &&
        error("No Ensembl transcript/protein candidates passed sequence validation.")
    if transcript_id !== nothing
        # A user-supplied transcript must match either the full or core ID.
        wanted = String(transcript_id)
        idx = findfirst(
            c -> c.transcript_id == wanted ||
                 strip_ensembl_version(c.transcript_id) == strip_ensembl_version(wanted),
            candidates)
        idx === nothing &&
            error("Requested transcript_id $(wanted) was not among the validated candidates.")
        return candidates[idx], String[]
    end
    warnings = String[]
    if length(candidates) > 1
        push!(warnings,
            "Multiple validated Ensembl transcripts were found; using $(first(candidates).transcript_id).")
    end
    return first(candidates), warnings
end

function _parse_transcript_lookup(data, transcript_id::AbstractString)::_EnsemblTranscriptLookup
    parent = get(data, "Parent", nothing)
    parent isa AbstractString ||
        error("Ensembl lookup for $(transcript_id) did not include a parent gene ID.")
    species = get(data, "species", nothing)
    return _EnsemblTranscriptLookup(;
        transcript_id = String(transcript_id),
        ensembl_gene_id = String(parent),
        species = species isa AbstractString ? String(species) : nothing
    )
end

function _resolve_transcript_metadata(transcript_id::AbstractString;
        _http_get_fn::Function = _http_get)::_EnsemblTranscriptLookup
    core = strip_ensembl_version(transcript_id)
    url = "$(_ENSEMBL_REST_BASE)/lookup/id/$(core)?expand=0"
    resp = _http_get_fn(url)
    resp === nothing &&
        error("Could not resolve Ensembl metadata for transcript $(transcript_id).")
    data = JSON.parse(decode_body(resp))
    return _parse_transcript_lookup(data, transcript_id)
end

"""
    resolve_transcript_gene(transcript_id) -> String

Look up the parent Ensembl gene ID for an Ensembl transcript.

# Arguments

- `transcript_id::AbstractString`: Ensembl transcript ID to look up.

# Throws

- `ErrorException`: if Ensembl metadata cannot be fetched or parsed.
"""
function resolve_transcript_gene(transcript_id::AbstractString)::String
    return _resolve_transcript_metadata(transcript_id).ensembl_gene_id
end

# Target Resolution
# -----------------

function _resolve_uniprot_target(input_id::AbstractString,
        sequence_dir::AbstractString,
        workdir::AbstractString;
        ensembl_gene_id::Union{Nothing, AbstractString},
        ensembl_protein_id::Union{Nothing, AbstractString},
        transcript_id::Union{Nothing, AbstractString},
        species::Union{Nothing, AbstractString},
        _uniprot_entry_fetcher::Function,
        _candidate_validator::Function)
    entry = _uniprot_entry_fetcher(input_id)
    candidates, uniprot_path = _candidate_validator(entry, sequence_dir)
    chosen, warnings = _choose_candidate(candidates, transcript_id)
    gene = ensembl_gene_id === nothing ? chosen.ensembl_gene_id : String(ensembl_gene_id)
    gene === nothing && error("Could not resolve an Ensembl gene ID for $(input_id).")
    protein = ensembl_protein_id === nothing ? chosen.ensembl_protein_id :
              String(ensembl_protein_id)
    protein_path = protein === nothing ? nothing :
                   joinpath(sequence_dir, "ensembl_proteins", "$(protein).fasta")
    return ResolvedTarget(;
        input_id = String(input_id),
        input_kind = :uniprot,
        uniprot_id = String(input_id),
        ensembl_gene_id = String(gene),
        transcript_id = chosen.transcript_id,
        ensembl_protein_id = protein,
        species = species === nothing ? chosen.species : String(species),
        uniprot_sequence_path = uniprot_path,
        ensembl_protein_sequence_path = protein_path,
        sequence_validated = chosen.sequence_validated,
        mapping_confirmed = chosen.mapping_confirmed,
        workdir = String(workdir),
        warnings
    )
end

function _resolve_transcript_target(input_id::AbstractString,
        workdir::AbstractString;
        uniprot_id::Union{Nothing, AbstractString},
        ensembl_gene_id::Union{Nothing, AbstractString},
        ensembl_protein_id::Union{Nothing, AbstractString},
        transcript_id::Union{Nothing, AbstractString},
        species::Union{Nothing, AbstractString},
        _transcript_metadata_resolver::Function)
    # Transcript input is intentionally lighter: ThorAxe needs the parent gene.
    tx = transcript_id === nothing ? String(input_id) : String(transcript_id)
    is_ensembl_transcript_id(tx) ||
        error("Transcript input $(tx) is not an Ensembl transcript ID.")
    metadata = (ensembl_gene_id === nothing || species === nothing) ?
               _transcript_metadata_resolver(tx) : nothing
    gene = ensembl_gene_id === nothing ? metadata.ensembl_gene_id :
           String(ensembl_gene_id)
    resolved_species = species === nothing ? metadata.species : String(species)
    return ResolvedTarget(;
        input_id = String(input_id),
        input_kind = :ensembl_transcript,
        uniprot_id = uniprot_id === nothing ? nothing : String(uniprot_id),
        ensembl_gene_id = gene,
        transcript_id = tx,
        ensembl_protein_id = ensembl_protein_id === nothing ? nothing :
                             String(ensembl_protein_id),
        species = resolved_species,
        sequence_validated = nothing,
        mapping_confirmed = nothing,
        workdir = String(workdir)
    )
end

"""
    resolve_target(input_id; workdir, uniprot_id=nothing, ensembl_gene_id=nothing,
                   ensembl_protein_id=nothing, transcript_id=nothing,
                   species=nothing) -> ResolvedTarget

Resolve one UniProt accession or Ensembl transcript ID into the IDs Iduna needs
for ThorAxe and validation.

# Arguments

- `input_id::AbstractString`: UniProt accession or Ensembl transcript ID to
  resolve.

# Keywords

- `workdir::AbstractString`: Iduna work directory where fetched sequences are
  saved.
- `uniprot_id = nothing`: UniProt accession to record when the main input is a
  transcript. When `nothing`, no UniProt accession is attached unless the input
  itself is UniProt.
- `ensembl_gene_id = nothing`: Ensembl gene ID to use instead of a fetched value.
  When `nothing`, Iduna fetches or infers the gene ID.
- `ensembl_protein_id = nothing`: Ensembl protein ID to use instead of a fetched
  value. When `nothing`, Iduna uses fetched metadata when available.
- `transcript_id = nothing`: Ensembl transcript ID to prefer when resolving a
  UniProt ID. When `nothing`, Iduna chooses from validated candidates.
- `species = nothing`: Species name to use instead of a fetched value. When
  `nothing`, Iduna uses fetched metadata when available.

# Returns

- A [`ResolvedTarget`](@ref) with resolved IDs, saved sequence paths, and warnings.
"""
function resolve_target(input_id::AbstractString;
        workdir::AbstractString,
        uniprot_id::Union{Nothing, AbstractString} = nothing,
        ensembl_gene_id::Union{Nothing, AbstractString} = nothing,
        ensembl_protein_id::Union{Nothing, AbstractString} = nothing,
        transcript_id::Union{Nothing, AbstractString} = nothing,
        species::Union{Nothing, AbstractString} = nothing,
        _transcript_metadata_resolver::Function = _resolve_transcript_metadata,
        _uniprot_entry_fetcher::Function = fetch_uniprot_entry,
        _candidate_validator::Function = _validate_candidates)
    kind = id_kind(input_id)
    sequence_dir = joinpath(workdir, "sequences")
    kind === :uniprot &&
        return _resolve_uniprot_target(input_id, sequence_dir, workdir;
            ensembl_gene_id, ensembl_protein_id, transcript_id, species,
            _uniprot_entry_fetcher, _candidate_validator)
    return _resolve_transcript_target(input_id, workdir;
        uniprot_id, ensembl_gene_id, ensembl_protein_id, transcript_id, species,
        _transcript_metadata_resolver)
end

end
