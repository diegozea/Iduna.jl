# Transcript MSA Assembly
# -----------------------

function _join_msas_consistently(
        lhs::MSAType, rhs::MSAType) where {MSAType <:
                                           AbstractMultipleSequenceAlignment}
    common = intersect(sequencenames(lhs), sequencenames(rhs))
    isempty(common) && error("Cannot join s-exon MSAs that share no sequence names.")
    return join_msas(lhs, rhs)
end

function _transcript_path_from_table(thoraxe_dir::AbstractString, transcript_id::AbstractString)
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
    return String(path_table.Path[only(matches)])
end

function _transcript_exon_ids(transcript_path::AbstractString)
    exon_tokens = split(transcript_path, "/")
    return [String(exon) for exon in exon_tokens if !(exon in ("", "start", "stop"))]
end

function _transcript_exon_file(thoraxe_dir::AbstractString, exon_id::AbstractString)
    return joinpath(thoraxe_dir, "msa", "msa_s_exon_$(exon_id).fasta")
end

function _row_matches_transcript(row, transcript_id::AbstractString)
    transcript_core = strip_ensembl_version(transcript_id)
    ids = split(String(row.TranscriptIDCluster), "/")
    return transcript_core in ids || String(transcript_id) in ids
end

function _s_exon_table(thoraxe_dir::AbstractString)
    path = joinpath(thoraxe_dir, "s_exon_table.csv")
    isfile(path) || error("Missing ThorAxe s_exon_table.csv at $(path).")
    return DataFrame(CSV.File(path))
end

function _s_exon_sequence_value(value)
    ismissing(value) && return ""
    sequence = replace(String(value), "*" => "")
    return isempty(strip(sequence)) ? "" : sequence
end

function _s_exon_sequence_from_table(thoraxe_dir::AbstractString,
        exon_id::AbstractString,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    table = _s_exon_table(thoraxe_dir)
    gene_core = strip_ensembl_version(gene_id)
    matches = filter(eachrow(table)) do row
        String(row.S_exonID) == String(exon_id) &&
            strip_ensembl_version(String(row.GeneID)) == gene_core &&
            _row_matches_transcript(row, transcript_id)
    end
    isempty(matches) &&
        error("Could not find sequence for $(exon_id) in ThorAxe s_exon_table.csv.")
    sequences = unique(_s_exon_sequence_value(row.S_exon_Sequence) for row in matches)
    length(sequences) == 1 ||
        error("ThorAxe s_exon_table.csv has inconsistent sequences for $(exon_id).")
    return only(sequences)
end

function _read_single_sequence_msa(name::AbstractString, sequence::AbstractString)
    return mktemp() do path, io
        write(io, ">", name, "\n", sequence, "\n")
        close(io)
        read_file(path, FASTA)
    end
end

function _transcript_exon_msa(thoraxe_dir::AbstractString,
        exon_id::AbstractString,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    exon_file = _transcript_exon_file(thoraxe_dir, exon_id)
    isfile(exon_file) && return read_file(exon_file, FASTA)
    startswith(exon_id, "0_") ||
        error("Expected ThorAxe s-exon MSA is missing: $(exon_file).")
    # ThorAxe can impute query-only s-exons as 0_* entries without separate MSA files.
    sequence = _s_exon_sequence_from_table(thoraxe_dir, exon_id, gene_id, transcript_id)
    isempty(sequence) && return nothing
    return _read_single_sequence_msa(strip_ensembl_version(gene_id), sequence)
end

function _transcript_exon_segments(thoraxe_dir::AbstractString,
        exon_ids::AbstractVector{<:AbstractString},
        gene_id::AbstractString,
        transcript_id::AbstractString)
    segments = Tuple{String, AbstractMultipleSequenceAlignment}[]
    for exon_id in exon_ids
        exon_msa = _transcript_exon_msa(thoraxe_dir, exon_id, gene_id, transcript_id)
        exon_msa === nothing && continue
        push!(segments, (String(exon_id), exon_msa))
    end
    return segments
end

# S-Exon Provenance
# -----------------

function _phylosofs_dir(thoraxe_dir::AbstractString)
    return joinpath(thoraxe_dir, "phylosofs")
end

function _phylosofs_transcripts_pir(thoraxe_dir::AbstractString)
    return joinpath(_phylosofs_dir(thoraxe_dir), "transcripts.pir")
end

function _phylosofs_s_exons_tsv(thoraxe_dir::AbstractString)
    return joinpath(_phylosofs_dir(thoraxe_dir), "s_exons.tsv")
end

function _has_phylosofs_outputs(thoraxe_dir::AbstractString)
    return isfile(_phylosofs_transcripts_pir(thoraxe_dir)) &&
           isfile(_phylosofs_s_exons_tsv(thoraxe_dir))
end

function _phylosofs_s_exon_code_map(thoraxe_dir::AbstractString)
    path = _phylosofs_s_exons_tsv(thoraxe_dir)
    isfile(path) || error("ThorAxe PhyloSofS s-exon map is missing: $(path).")
    code_map = Pair{Char, String}[]
    for line in eachline(path)
        isempty(strip(line)) && continue
        # PhyloSofS writes the stable s-exon ID and its one-character alignment code.
        fields = split(line, '\t'; limit = 2)
        length(fields) == 2 ||
            error("Invalid PhyloSofS s-exon map line in $(path): $(repr(line)).")
        symbol = only(String(fields[2]))
        push!(code_map, symbol => String(fields[1]))
    end
    isempty(code_map) && error("ThorAxe PhyloSofS s-exon map is empty: $(path).")
    return code_map
end

function _phylosofs_transcript_tokens(seq)
    tokens = String[]
    for token in split(sequence_id(seq))
        for clustered_token in split(String(token), '/')
            transcript = strip_ensembl_version(clustered_token)
            isempty(transcript) || push!(tokens, transcript)
        end
    end
    return tokens
end

function _phylosofs_transcript_sequence(thoraxe_dir::AbstractString,
        transcript_id::AbstractString)
    path = _phylosofs_transcripts_pir(thoraxe_dir)
    isfile(path) || error("ThorAxe PhyloSofS transcript PIR is missing: $(path).")
    transcript_core = strip_ensembl_version(transcript_id)
    matches = filter(read_file(path, PIRSequences)) do seq
        transcript_core in _phylosofs_transcript_tokens(seq)
    end
    isempty(matches) &&
        error("Transcript $(transcript_id) was not found in ThorAxe PhyloSofS PIR $(path).")
    length(matches) == 1 ||
        error("Transcript $(transcript_id) matched more than one entry in ThorAxe PhyloSofS PIR $(path).")
    return only(matches)
end

function _phylosofs_transcript_symbols(thoraxe_dir::AbstractString,
        transcript_id::AbstractString)
    seq = _phylosofs_transcript_sequence(thoraxe_dir, transcript_id)
    symbols = getannotsequence(seq, "Title")
    length(symbols) == length(seq) ||
        error("ThorAxe PhyloSofS annotation length for $(transcript_id) does not match its sequence length.")
    return String(symbols)
end

function _s_exon_symbol_lookup(code_map::AbstractVector{<:Pair})
    return Dict(String(s_exon_id) => symbol for (symbol, s_exon_id) in code_map)
end

function _complete_s_exon_code_map(code_map::AbstractVector{<:Pair},
        exon_ids::AbstractVector{<:AbstractString})
    complete_map = Pair{Char, String}[symbol => String(s_exon_id)
                                      for (symbol, s_exon_id) in code_map]
    mapped_exons = Set(String(s_exon_id) for (_, s_exon_id) in complete_map)
    used_symbols = Set(symbol for (symbol, _) in complete_map)
    for exon_id in exon_ids
        exon = String(exon_id)
        exon in mapped_exons && continue
        startswith(exon, "0_") ||
            error("ThorAxe PhyloSofS s-exon map has no symbol for $(exon).")
        # Some reference-only s-exons are absent from PhyloSofS, so assign a safe symbol.
        symbol_idx = findfirst(symbol -> !(symbol in used_symbols), _PHYLOSOFS_FALLBACK_SYMBOLS)
        symbol_idx === nothing &&
            error("No unused PhyloSofS-compatible symbol is available for $(exon).")
        symbol = _PHYLOSOFS_FALLBACK_SYMBOLS[symbol_idx]
        push!(complete_map, symbol => exon)
        push!(mapped_exons, exon)
        push!(used_symbols, symbol)
    end
    return complete_map
end

function _s_exon_symbol_for_reference_residue(residue::Char, exon_symbol::Char)
    return residue == '.' ? '.' : exon_symbol
end

function _consume_phylosofs_symbol(symbol_state,
        symbols::AbstractString,
        exon_id::AbstractString,
        exon_by_symbol::AbstractDict,
        transcript_id::AbstractString)
    symbol_state === nothing &&
        error("ThorAxe PhyloSofS annotation for $(transcript_id) is shorter than the transcript MSA reference sequence.")
    observed_symbol, state = symbol_state
    pir_exon = get(exon_by_symbol, observed_symbol, nothing)
    (pir_exon === nothing || pir_exon == String(exon_id)) ||
        error("ThorAxe PhyloSofS annotation for $(transcript_id) does not match s-exon $(exon_id).")
    return observed_symbol, iterate(symbols, state)
end

function _write_projected_exon_symbols!(io::IO,
        exon_msa,
        exon_id::AbstractString,
        gene_id::AbstractString,
        symbols::AbstractString,
        exon_symbol::Char,
        exon_by_symbol::AbstractDict,
        symbol_state,
        transcript_id::AbstractString)
    reference = resolve_sequence_name(exon_msa, gene_id)
    reference === nothing &&
        error("Could not find $(gene_id) in the s-exon MSA for $(exon_id).")
    for residue in stringsequence(exon_msa, reference)
        if residue == '.' || residue == '-'
            write(io, _s_exon_symbol_for_reference_residue(residue, exon_symbol))
        else
            _,
            symbol_state = _consume_phylosofs_symbol(
                symbol_state, symbols, exon_id, exon_by_symbol, transcript_id)
            write(io, exon_symbol)
        end
    end
    return symbol_state
end

function _project_phylosofs_symbols(exon_msas::AbstractVector,
        exon_ids::AbstractVector{<:AbstractString},
        gene_id::AbstractString,
        symbols::AbstractString,
        code_map::AbstractVector{<:Pair},
        transcript_id::AbstractString)
    symbol_by_exon = _s_exon_symbol_lookup(code_map)
    exon_by_symbol = Dict(symbol => String(s_exon_id) for (symbol, s_exon_id) in code_map)
    symbol_state = iterate(symbols)
    io = IOBuffer()
    for (exon_id, exon_msa) in zip(exon_ids, exon_msas)
        exon_symbol = get(symbol_by_exon, String(exon_id), nothing)
        exon_symbol === nothing &&
            error("ThorAxe PhyloSofS s-exon map has no symbol for $(exon_id).")
        # Walk through exon MSAs in transcript order and copy one s-exon symbol per column.
        symbol_state = _write_projected_exon_symbols!(
            io, exon_msa, exon_id, gene_id, symbols, exon_symbol, exon_by_symbol,
            symbol_state, transcript_id)
    end
    symbol_state === nothing ||
        error("ThorAxe PhyloSofS annotation for $(transcript_id) is longer than the transcript MSA reference sequence.")
    return String(take!(io))
end

function _maybe_set_s_exon_annotations!(transcript_msa::AbstractMultipleSequenceAlignment,
        exon_msas::AbstractVector,
        exon_ids::AbstractVector{<:AbstractString},
        thoraxe_dir::AbstractString,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    _has_phylosofs_outputs(thoraxe_dir) || return transcript_msa
    code_map = _complete_s_exon_code_map(_phylosofs_s_exon_code_map(thoraxe_dir), exon_ids)
    codes = _project_phylosofs_symbols(
        exon_msas, exon_ids, gene_id,
        _phylosofs_transcript_symbols(thoraxe_dir, transcript_id),
        code_map, transcript_id)
    set_s_exon_annotations!(transcript_msa, codes, code_map)
    return transcript_msa
end

function _transcript_msa_species(thoraxe_dir::AbstractString,
        transcript_msa::AbstractMultipleSequenceAlignment)
    s_exon_path = joinpath(thoraxe_dir, "s_exon_table.csv")
    isfile(s_exon_path) || return fill("unknown", nsequences(transcript_msa))

    table = DataFrame(CSV.File(s_exon_path))
    lookup = Dict(String(row.GeneID) => String(row.Species) for row in eachrow(table))
    return [get(lookup, String(name), "unknown") for name in sequencenames(transcript_msa)]
end

# Public MSA Assembly
# -------------------

"""
    assemble_transcript_msa(thoraxe_dir, gene_id, transcript_id)

Rebuild the transcript-level MSA for one transcript from ThorAxe s-exon outputs.

# Arguments

- `thoraxe_dir::AbstractString`: ThorAxe output directory.
- `gene_id::AbstractString`: Ensembl gene ID.
- `transcript_id::AbstractString`: Ensembl transcript ID to rebuild.

# Returns

- A tuple `(msa, species)`, where `msa` is the reconstructed MSA and `species`
  gives one species label per sequence.

# Throws

- `ErrorException`: if the transcript path, s-exons, or reference sequence cannot
  be found.
"""
function assemble_transcript_msa(thoraxe_dir::AbstractString,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    transcript_path = _transcript_path_from_table(thoraxe_dir, transcript_id)
    path_exon_ids = _transcript_exon_ids(transcript_path)
    isempty(path_exon_ids) && error("No s-exons were found for $(transcript_id).")

    exon_segments = _transcript_exon_segments(
        thoraxe_dir, path_exon_ids, gene_id, transcript_id)
    isempty(exon_segments) &&
        error("No non-empty s-exons were found for $(transcript_id).")
    exon_ids = first.(exon_segments)
    exon_msas = last.(exon_segments)
    transcript_msa = reduce(_join_msas_consistently, exon_msas)

    reference = resolve_sequence_name(transcript_msa, gene_id)
    reference === nothing &&
        error("Could not find $(gene_id) in the reconstructed transcript MSA.")
    setreference!(transcript_msa, reference)
    _maybe_set_s_exon_annotations!(
        transcript_msa, exon_msas, exon_ids, thoraxe_dir, gene_id, transcript_id)

    return transcript_msa, _transcript_msa_species(thoraxe_dir, transcript_msa)
end
