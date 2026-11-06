# HHsuite Alignment Parsing
# -------------------------

function _hhsuite_residue_indices(seq::AbstractString, start::Integer, stop::Integer)
    start >= 1 ||
        throw(ArgumentError("HHsuite residue positions must start at 1 or later."))
    indices = Int[]
    current = Int(start)
    for residue in seq
        if residue in ('-', '.')
            push!(indices, 0)
        else
            push!(indices, current)
            current += 1
        end
    end
    current - 1 == stop ||
        throw(ArgumentError("HHsuite residue positions disagree with the aligned sequence."))
    return indices
end

function _parse_hhsuite_sequence_row(line::AbstractString)
    m = match(r"^[QT]\s+\S+\s+(-?\d+)\s+(\S+)\s+(-?\d+)(?:\s|$)", line)
    m === nothing && return nothing

    start = parse(Int, m.captures[1])
    seq = m.captures[2]
    stop = parse(Int, m.captures[3])
    # PrintHHR uses %4i for coordinates, but keeps the comparison row fixed when
    # a coordinate needs more digits. Its symbols still begin at the original column.
    comparison_start = m.offsets[2] - max(0, ndigits(start) - 4)
    return (;
        cols = comparison_start:(comparison_start + ncodeunits(seq) - 1),
        indices = _hhsuite_residue_indices(seq, start, stop))
end

function _validated_hhsuite_sequence(lines, row_indices, label)
    first_sequence = nothing
    for index in row_indices
        line = lines[index]
        # Secondary-structure annotations have no residue coordinates. A malformed
        # sequence row must not be skipped in favour of another sequence or consensus.
        occursin(r"^[QT]\s+(?:ss_pred|ss_conf|ss_dssp|sa_dssp)\S*\s+\S+\s*$", line) &&
            continue
        parsed = _parse_hhsuite_sequence_row(line)
        parsed === nothing &&
            throw(ArgumentError("Malformed HHsuite $(label) sequence row: $(line)"))
        if first_sequence === nothing
            first_sequence = parsed
        elseif length(parsed.indices) != length(first_sequence.indices)
            throw(ArgumentError("HHsuite $(label) sequence rows have different alignment widths."))
        end
    end
    first_sequence === nothing &&
        throw(ArgumentError("HHsuite alignment has no recognized $(label) sequence row."))
    return first_sequence
end

function _append_hhsuite_block!(positions, codes, lines; summary_pending = false)
    query_rows = findall(line -> startswith(line, "Q "), lines)
    target_rows = findall(line -> startswith(line, "T "), lines)
    if isempty(query_rows) && isempty(target_rows)
        return false
    end
    summary_pending && throw(ArgumentError("HHsuite hit has no summary header."))
    query = _validated_hhsuite_sequence(lines, query_rows, "query")
    target = _validated_hhsuite_sequence(lines, target_rows, "target")

    # Count positions from the first displayed query sequence, but read comparison
    # symbols after the last query row (which can be the consensus).
    comparison_index = last(query_rows) + 1
    comparison_index < first(target_rows) ||
        throw(ArgumentError("HHsuite alignment has no comparison row."))
    length(query.indices) == length(target.indices) ||
        throw(ArgumentError("HHsuite query and target alignment widths differ."))
    comparison_line = lines[comparison_index]
    if !startswith(comparison_line, " ") || lastindex(comparison_line) < last(query.cols)
        throw(ArgumentError("HHsuite comparison row is malformed or truncated."))
    end
    symbols = comparison_line[query.cols]
    all(c -> c in " |+.-=", symbols) ||
        throw(ArgumentError("HHsuite comparison row contains unrecognized symbols."))

    append!(positions, query.indices)
    append!(codes, symbols)
    return true
end

function _require_hhsuite_block_separator(block)
    any(row -> startswith(row, "Q ") || startswith(row, "T "), block) &&
        throw(ArgumentError("HHsuite alignment blocks must be separated."))
    return nothing
end

function _require_complete_hhsuite_hit(hit_recognized, summary_pending = false)
    summary_pending && throw(ArgumentError("HHsuite hit has no summary header."))
    hit_recognized || throw(ArgumentError("HHsuite hit has no complete alignment block."))
    return nothing
end

function _get_codes(output::AbstractString)
    lines = split(output, '\n')
    start = findfirst(line -> startswith(line, "Probab=") || occursin(r"^No \d+$", line), lines)
    start === nothing && throw(ArgumentError("HHsuite output has no recognized alignment."))
    positions = Int[]
    codes = Char[]
    block = SubString{String}[]
    hit_recognized = false
    # No N and its following Probab line belong to the same hit header.
    summary_pending = startswith(lines[start], "No ")
    # HHsuite 3.3 PrintHHR separates blocks with empty lines; Confidence is optional.
    # Preserve all-space comparison rows: only truly empty lines end a block.
    for line in @view(lines[(start + 1):end])
        if occursin(r"^No \d+$", line)
            _require_hhsuite_block_separator(block)
            _require_complete_hhsuite_hit(hit_recognized, summary_pending)
            hit_recognized = false
            summary_pending = true
            empty!(block)
        elseif startswith(line, "Probab=")
            _require_hhsuite_block_separator(block)
            summary_pending || _require_complete_hhsuite_hit(hit_recognized)
            hit_recognized = false
            summary_pending = false
            empty!(block)
        elseif isempty(line) || startswith(line, "Confidence")
            hit_recognized |= _append_hhsuite_block!(positions, codes, block; summary_pending)
            empty!(block)
        else
            push!(block, line)
        end
    end
    hit_recognized |= _append_hhsuite_block!(positions, codes, block; summary_pending)
    _require_complete_hhsuite_hit(hit_recognized, summary_pending)
    return positions, codes
end

function _identity_counts_from_codes(positions::Vector{Int}, codes::Vector{Char})
    length(positions) == length(codes) ||
        throw(ArgumentError("HHsuite positions and comparison symbols have different lengths."))
    all(>=(0), positions) || throw(ArgumentError("HHsuite positions cannot be negative."))
    seen = Dict{Int, Bool}()
    for (pos, code) in zip(positions, codes)
        pos == 0 && continue
        # A reference position counts as matched if any compared row marks it with '|'.
        seen[pos] = get(seen, pos, false) || code == '|'
    end
    matched = count(identity, values(seen))
    return IdentityCounts(matched, length(seen))
end

function _parse_hhalign_float(line::AbstractString, key::AbstractString)
    pattern = Regex("(?:^|\\s)$(key)=(\\S+)")
    m = match(pattern, line)
    m === nothing && return missing
    value = tryparse(Float64, m.captures[1])
    return value === nothing ? missing : value
end

function _hhalign_summary(output::AbstractString)
    lines = split(output, '\n')
    idx = findfirst(l -> startswith(l, "Probab="), lines)
    idx === nothing &&
        return (; probability = missing, e_value = missing, hhalign_score = missing)
    line = lines[idx]
    return (;
        probability = _parse_hhalign_float(line, "Probab"),
        e_value = _parse_hhalign_float(line, "E-value"),
        hhalign_score = _parse_hhalign_float(line, "Score"))
end

# HHsuite Execution
# -----------------

function _run_hhsuite_alignment(reference_fasta::AbstractString,
        sample_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        label::AbstractString = "sample")
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
        return output
    end
end

# Public Score Functions
# ----------------------

"""
    hhsuite_identity_score(reference_msa_fasta, sample_msa_fasta;
                           logs_dir=nothing, label="sample")

Compare two MSA FASTA files with HHsuite and return identity-count fields for EPLI.

# Arguments

- `reference_msa_fasta::AbstractString`: full reference MSA in FASTA format.
- `sample_msa_fasta::AbstractString`: sampled MSA in FASTA format.

# Keywords

- `logs_dir = nothing`: optional directory where the HHalign output is copied.
- `label::AbstractString = "sample"`: label used in the copied log file name.

# Returns

An `IdentityScore` with validated counts and separate tool diagnostics.
"""
function hhsuite_identity_score(reference_msa_fasta::AbstractString,
        sample_msa_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        label::AbstractString = "sample")
    output = _run_hhsuite_alignment(reference_msa_fasta, sample_msa_fasta;
        logs_dir, label)
    positions, codes = _get_codes(output)
    counts = _identity_counts_from_codes(positions, codes)
    summary = _hhalign_summary(output)
    return IdentityScore(counts; diagnostics = merge(summary, (;
        score_name = "hhsuite_identity_score")))
end

"""
    hhsuite_profile_score(reference_msa_fasta, sample_msa_fasta;
                          logs_dir=nothing, label="sample")

Compare two MSA FASTA files with HHsuite and return the HHalign profile score.
"""
function hhsuite_profile_score(reference_msa_fasta::AbstractString,
        sample_msa_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        label::AbstractString = "sample")
    output = _run_hhsuite_alignment(reference_msa_fasta, sample_msa_fasta;
        logs_dir, label)
    summary = _hhalign_summary(output)
    ismissing(summary.hhalign_score) &&
        error("Could not parse HHalign profile score.")
    return RawScore(summary.hhalign_score;
        diagnostics = merge(summary, (;
            score_name = "hhsuite_profile_score")))
end
