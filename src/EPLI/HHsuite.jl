# HHsuite Alignment Parsing
# -------------------------

function _hhsuite_query_indices(seq::AbstractString, start::Integer, stop::Integer)
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
    current - 1 == stop || error("Could not parse HHsuite alignment positions.")
    return indices
end

function _parse_hhsuite_query_segment(line::AbstractString)
    m = match(r"\S+\s+(\d+)\s+(\S+)\s+(\d+)", line)
    m === nothing && return nothing

    start = parse(Int, m.captures[1])
    seq = m.captures[2]
    stop = parse(Int, m.captures[3])
    return (;
        cols = m.offsets[2]:(m.offsets[2] + ncodeunits(seq) - 1),
        indices = _hhsuite_query_indices(seq, start, stop))
end

function _append_hhsuite_code_line!(positions, codes, line::AbstractString, cols, indices)
    occursin(r"^\s", line) || return nothing
    isempty(indices) && return nothing
    (cols === nothing || isempty(cols)) && return nothing
    append!(positions, indices)
    append!(codes, collect(line[cols]))
    return nothing
end

function _hhsuite_alignment_row(lines, prefix)
    idx = findfirst(
        line -> startswith(line, prefix) &&
                occursin(r"^\S+\s+\S+\s+\d+\s+\S+\s+\d+", line), lines)
    idx === nothing &&
        throw(ArgumentError("HHsuite alignment has no recognized $(prefix)row."))
    parsed = _parse_hhsuite_query_segment(lines[idx])
    return idx, parsed
end

function _hhsuite_block_codes!(positions, codes, lines)
    _query_idx, query = _hhsuite_alignment_row(lines, "Q ")
    query_end = findlast(line -> startswith(line, "Q "), lines)
    _target_idx, target = _hhsuite_alignment_row(lines, "T ")
    query_end + 1 < findfirst(line -> startswith(line, "T "), lines) ||
        throw(ArgumentError("HHsuite alignment has no comparison row."))
    length(query.indices) == length(target.indices) ||
        throw(ArgumentError("HHsuite query and target alignment widths differ."))
    code_line = lines[query_end + 1]
    startswith(code_line, " ") && lastindex(code_line) >= last(query.cols) ||
        throw(ArgumentError("HHsuite comparison row is malformed or truncated."))
    symbols = code_line[query.cols]
    all(c -> c in " |+.-=", symbols) ||
        throw(ArgumentError("HHsuite comparison row contains unrecognized symbols."))
    _append_hhsuite_code_line!(positions, codes, code_line, query.cols, query.indices)
    return nothing
end

function _append_hhsuite_block!(positions, codes, lines)
    any(line -> startswith(line, "Q ") || startswith(line, "T "), lines) || return false
    _hhsuite_block_codes!(positions, codes, lines)
    return true
end

function _get_codes(output::AbstractString)
    lines = split(output, '\n')
    start = findfirst(line -> startswith(line, "Probab="), lines)
    start === nothing && throw(ArgumentError("HHsuite output has no recognized alignment."))
    positions = Int[]
    codes = Char[]
    block = SubString{String}[]
    recognized = false
    # HHsuite 3.3 PrintHHR separates blocks with empty lines; Confidence is optional.
    # Preserve all-space comparison rows: only truly empty lines end a block.
    for line in @view(lines[(start + 1):end])
        if isempty(line) || startswith(line, "Confidence")
            recognized |= _append_hhsuite_block!(positions, codes, block)
            empty!(block)
        else
            push!(block, line)
        end
    end
    recognized |= _append_hhsuite_block!(positions, codes, block)
    recognized || throw(ArgumentError("HHsuite output has no complete alignment block."))
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
    pattern = Regex("$(key)=([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)")
    m = match(pattern, line)
    m === nothing && return missing
    return parse(Float64, m.captures[1])
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
