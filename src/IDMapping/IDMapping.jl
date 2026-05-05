module IDMapping

using HTTP
using JSON3

using ..Utils: ResolvedTarget, decode_body, fasta_sequence, format_fasta, id_kind,
               is_ensembl_transcript_id, protein_alignment_stats, strip_ensembl_version,
               write_text

export EnsemblCandidate,
       fetch_uniprot_entry,
       resolve_target,
       resolve_transcript_gene,
       sequences_match

const _UNIPROT_BASE = "https://rest.uniprot.org"
const _ENSEMBL_REST_BASE = "https://rest.ensembl.org"
const _JSON_HEADERS = ["Accept" => "application/json"]
const _FASTA_HEADERS = ["Accept" => "text/x-fasta"]

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

function _http_get(url::AbstractString, headers = _JSON_HEADERS; retries::Int = 4, sleep_seconds::Real = 1.5)
    last_status = nothing
    for attempt in 1:max(retries, 1)
        resp = HTTP.get(url; headers = headers, retry = false, status_exception = false)
        if resp.status == 200
            return resp
        elseif resp.status in (429, 500, 502, 503, 504)
            last_status = resp.status
            sleep(sleep_seconds * attempt)
        else
            return nothing
        end
    end
    @warn "HTTP request failed after retries." url status=last_status
    return nothing
end

function fetch_uniprot_entry(uniprot_id::AbstractString)::_UniProtEntry
    url = "$(_UNIPROT_BASE)/uniprotkb/$(strip(String(uniprot_id))).json"
    resp = _http_get(url)
    resp === nothing && error("Could not fetch UniProt metadata for $(uniprot_id).")
    data = JSON3.read(decode_body(resp))
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

function _parse_xrefs(data)
    gene_ids = String[]
    transcripts = String[]
    transcript_to_protein = Dict{String, String}()
    transcript_to_gene = Dict{String, String}()
    transcript_to_isoform = Dict{String, String}()

    for ref in get(data, "uniProtKBCrossReferences", Any[])
        get(ref, "database", "") == "Ensembl" || continue
        transcript = get(ref, "id", nothing)
        if transcript isa AbstractString
            transcript = String(transcript)
            push!(transcripts, transcript)
            isoform = get(ref, "isoformId", nothing)
            isoform isa AbstractString &&
                (transcript_to_isoform[transcript] = String(isoform))
        end

        for prop in get(ref, "properties", Any[])
            key = get(prop, "key", "")
            value = get(prop, "value", nothing)
            if key == "GeneId" && value isa AbstractString
                push!(gene_ids, String(value))
                transcript isa AbstractString &&
                    (transcript_to_gene[String(transcript)] = String(value))
            elseif key == "ProteinId" && value isa AbstractString &&
                   transcript isa AbstractString
                transcript_to_protein[String(transcript)] = String(value)
            end
        end
    end

    gene_ids = unique(gene_ids)
    transcripts = unique(transcripts)
    if !isempty(gene_ids)
        default_gene = first(gene_ids)
        for transcript in transcripts
            haskey(transcript_to_gene, transcript) ||
                (transcript_to_gene[transcript] = default_gene)
        end
    end
    return gene_ids, transcripts, transcript_to_protein, transcript_to_gene,
    transcript_to_isoform
end

function _fetch_ensembl_protein_sequence(protein_id::AbstractString)::Union{Nothing, String}
    core = strip_ensembl_version(protein_id)
    url = "$(_ENSEMBL_REST_BASE)/sequence/id/$(core)?type=protein"
    resp = _http_get(url, _FASTA_HEADERS)
    resp === nothing && return nothing
    return fasta_sequence(decode_body(resp))
end

function _fetch_uniprot_fasta_sequence(uniprot_id::AbstractString)::String
    url = "$(_UNIPROT_BASE)/uniprotkb/$(strip(String(uniprot_id))).fasta"
    resp = _http_get(url, _FASTA_HEADERS; retries = 5, sleep_seconds = 2.0)
    resp === nothing && error("Could not fetch UniProt FASTA for $(uniprot_id).")
    seq = fasta_sequence(decode_body(resp))
    seq === nothing && error("UniProt FASTA for $(uniprot_id) did not contain a sequence.")
    return seq
end

function sequences_match(uniprot_seq::AbstractString, ensembl_seq::AbstractString)::Bool
    protein_alignment_stats(ensembl_seq, uniprot_seq).identical
end

function _validate_candidates(entry::_UniProtEntry, sequence_dir::AbstractString)
    uniprot_seq = entry.protein_sequence === nothing ?
                  _fetch_uniprot_fasta_sequence(entry.id) : entry.protein_sequence
    uniprot_path = joinpath(sequence_dir, "uniprot", "$(entry.id).fasta")
    write_text(uniprot_path, format_fasta(entry.id, uniprot_seq))

    candidates = EnsemblCandidate[]
    for transcript in entry.transcript_ids
        protein_id = get(entry.transcript_to_protein, transcript, nothing)
        protein_id === nothing && continue
        protein_seq = _fetch_ensembl_protein_sequence(protein_id)
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

function _choose_candidate(candidates::Vector{EnsemblCandidate}, transcript_id::Union{
        Nothing, AbstractString})
    isempty(candidates) &&
        error("No Ensembl transcript/protein candidates passed sequence validation.")
    if transcript_id !== nothing
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

function _resolve_transcript_metadata(transcript_id::AbstractString)::_EnsemblTranscriptLookup
    core = strip_ensembl_version(transcript_id)
    url = "$(_ENSEMBL_REST_BASE)/lookup/id/$(core)?expand=0"
    resp = _http_get(url)
    resp === nothing &&
        error("Could not resolve Ensembl metadata for transcript $(transcript_id).")
    data = JSON3.read(decode_body(resp))
    return _parse_transcript_lookup(data, transcript_id)
end

function resolve_transcript_gene(transcript_id::AbstractString)::String
    return _resolve_transcript_metadata(transcript_id).ensembl_gene_id
end

function resolve_target(input_id::AbstractString;
        workdir::AbstractString,
        uniprot_id::Union{Nothing, AbstractString} = nothing,
        ensembl_gene_id::Union{Nothing, AbstractString} = nothing,
        ensembl_protein_id::Union{Nothing, AbstractString} = nothing,
        transcript_id::Union{Nothing, AbstractString} = nothing,
        species::Union{Nothing, AbstractString} = nothing,
        _transcript_metadata_resolver::Function = _resolve_transcript_metadata)
    kind = id_kind(input_id)
    sequence_dir = joinpath(workdir, "sequences")

    if kind === :uniprot
        entry = fetch_uniprot_entry(input_id)
        candidates, uniprot_path = _validate_candidates(entry, sequence_dir)
        chosen, warnings = _choose_candidate(candidates, transcript_id)
        gene = ensembl_gene_id === nothing ? chosen.ensembl_gene_id :
               String(ensembl_gene_id)
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
            warnings
        )
    end

    # Transcript input is intentionally lighter: ThorAxe needs the parent gene,
    # and UniProt mapping is not required for this path.
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
        mapping_confirmed = nothing
    )
end

end
