function _hhsuite_query_indices(seq::AbstractString, start::Integer, stop::Integer)
    indices = Int[]
    current = Int(start)
    for residue in seq
        if residue == '-'
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
        cols = findfirst(seq, line),
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

function _get_codes(output::AbstractString)
    in_alignment = false
    query_line = true
    cols = 1:0
    indices = Int[]
    positions = Int[]
    codes = Char[]
    for line in split(output, '\n')
        if startswith(line, "Probab=")
            in_alignment = true
            continue
        end
        in_alignment || continue
        isempty(line) && continue
        # HHalign prints one query row followed by rows of symbols for each block.
        if query_line
            parsed = _parse_hhsuite_query_segment(line)
            if parsed !== nothing
                cols = parsed.cols
                indices = parsed.indices
            end
            query_line = false
        elseif startswith(line, "Confidence")
            query_line = true
        else
            _append_hhsuite_code_line!(positions, codes, line, cols, indices)
        end
    end
    return positions, codes
end

function _identity_counts_from_codes(positions::Vector{Int}, codes::Vector{Char})
    seen = Dict{Int, Bool}()
    for (pos, code) in zip(positions, codes)
        pos == 0 && continue
        # A reference position counts as matched if any compared row marks it with '|'.
        seen[pos] = get(seen, pos, false) || code == '|'
    end
    matched = count(identity, values(seen))
    return (;
        matched_positions = matched,
        comparable_positions = length(seen))
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

A named tuple with `raw_score`, `matched_positions`, and `comparable_positions`.
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
    return merge(summary,
        (;
            score_name = "hhsuite_identity_score",
            raw_score = Float64(counts.matched_positions),
            matched_positions = counts.matched_positions,
            comparable_positions = counts.comparable_positions))
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
    return merge(summary, (;
        score_name = "hhsuite_profile_score",
        raw_score = Float64(summary.hhalign_score)))
end
