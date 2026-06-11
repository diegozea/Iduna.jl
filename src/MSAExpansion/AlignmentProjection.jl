# Alignment Projection
# --------------------

function _reorder_alignment(msa::AbstractMultipleSequenceAlignment, priority::Vector{String})
    name_to_idx = Dict(String(name) => idx for (idx, name) in enumerate(sequencenames(msa)))
    priority_idx = Int[]
    seen = Set{Int}()
    for name in priority
        idx = get(name_to_idx, String(name), nothing)
        idx === nothing && continue
        push!(priority_idx, idx)
        push!(seen, idx)
    end
    remaining = [idx for idx in 1:nsequences(msa) if !(idx in seen)]
    return msa[vcat(priority_idx, remaining), :]
end

function _is_rf_match(rf_char::Char)
    return !(rf_char == '.' || rf_char == '-' || rf_char == ' ')
end

function _rf_match_state_mask(rf::AbstractString, ncols::Integer)
    if isempty(rf)
        return nothing
    end
    length(rf) == ncols ||
        error("RF annotation has $(length(rf)) characters, but the MSA has $(ncols) columns.")
    return [_is_rf_match(rf_char) for rf_char in rf]
end

function _aligned_match_state_mask(aligned::AbstractString, ncols::Integer)
    if isempty(aligned)
        return nothing
    end
    length(aligned) == ncols ||
        error("Aligned annotation has $(length(aligned)) characters, but the MSA has $(ncols) columns.")
    return [aligned_char == '1' ? true :
            aligned_char == '0' ? false :
            error("Aligned annotation contains $(repr(aligned_char)); expected '0' or '1'.")
            for aligned_char in aligned]
end

function _match_state_mask(msa::AbstractMultipleSequenceAlignment;
        default_aligned::Bool)
    ncols = ncolumns(msa)
    aligned_mask = _aligned_match_state_mask(getannotcolumn(msa, "Aligned", ""), ncols)
    aligned_mask === nothing || return aligned_mask
    rf_mask = _rf_match_state_mask(getannotcolumn(msa, "RF", ""), ncols)
    rf_mask === nothing || return rf_mask
    # If no annotation says which columns are seed columns, use the caller's default.
    return fill(default_aligned, ncols)
end

function _project_s_exon_codes(msa::AbstractMultipleSequenceAlignment,
        seed_match_codes::AbstractString)
    match_mask = _match_state_mask(msa; default_aligned = false)
    seed_match_chars = collect(seed_match_codes)
    io = IOBuffer(; sizehint = ncolumns(msa))
    seed_idx = 0
    for is_match in match_mask
        if is_match
            # Match columns inherit s-exon labels from the seed; inserted columns do not.
            seed_idx += 1
            write(io, seed_idx <= length(seed_match_chars) ? seed_match_chars[seed_idx] :
                      '.')
        else
            write(io, '.')
        end
    end
    return String(take!(io))
end

function _seed_match_s_exon_codes(seed_codes::AbstractString,
        annotated_seed::AbstractMultipleSequenceAlignment)
    match_mask = _match_state_mask(annotated_seed; default_aligned = true)
    length(match_mask) == length(seed_codes) ||
        error("Annotated seed match mask has $(length(match_mask)) characters, but SExonCode has $(length(seed_codes)) characters.")
    io = IOBuffer(; sizehint = length(seed_codes))
    for (is_match, code) in zip(match_mask, seed_codes)
        is_match && write(io, code)
    end
    return String(take!(io))
end

function _with_seed_match_s_exon_codes(archived, annotated_seed_path::AbstractString)
    annotated_seed = read_file(annotated_seed_path, Stockholm; keepinserts = true)
    return merge(archived,
        (;
            seed_match_s_exon_codes = _seed_match_s_exon_codes(
            archived.seed_s_exon_codes, annotated_seed)))
end

function _restore_s_exon_annotations!(msa::AbstractMultipleSequenceAlignment, archived)
    seed_match_codes = :seed_match_s_exon_codes in propertynames(archived) ?
                       archived.seed_match_s_exon_codes : archived.seed_s_exon_codes
    codes = _project_s_exon_codes(msa, seed_match_codes)
    set_s_exon_annotations!(msa, codes, archived.seed_s_exon_code_map)
    return msa
end

function _write_expansion_s_exon_blocks(path::AbstractString,
        full_alignment::AbstractMultipleSequenceAlignment,
        match_alignment::AbstractMultipleSequenceAlignment,
        pid::Real)
    write_s_exon_blocks_tsv(path, match_alignment;
        alignment = "expanded_match",
        pid = Float64(pid))
    write_s_exon_blocks_tsv(path, full_alignment;
        alignment = "expanded_full",
        pid = Float64(pid),
        append = true)
    return path
end

function _ensure_alignment_s_exon_blocks(path::AbstractString,
        full_stockholm::AbstractString,
        match_stockholm::AbstractString,
        pid::Real;
        full_label::AbstractString,
        match_label::AbstractString)
    isfile(path) && return path
    (isfile(full_stockholm) && isfile(match_stockholm)) || return path
    full_alignment = read_file(full_stockholm, Stockholm; keepinserts = true)
    match_alignment = read_file(match_stockholm, Stockholm; keepinserts = true)
    write_s_exon_blocks_tsv(path, match_alignment;
        alignment = match_label,
        pid = Float64(pid))
    write_s_exon_blocks_tsv(path, full_alignment;
        alignment = full_label,
        pid = Float64(pid),
        append = true)
    return path
end

function _ensure_expansion_s_exon_blocks(outputs::NamedTuple, pid::Real)
    _ensure_alignment_s_exon_blocks(outputs.s_exon_blocks_tsv,
        outputs.full_stockholm,
        outputs.match_stockholm,
        pid;
        full_label = "expanded_full",
        match_label = "expanded_match")
    if :centroid_s_exon_blocks_tsv in propertynames(outputs)
        _ensure_alignment_s_exon_blocks(outputs.centroid_s_exon_blocks_tsv,
            outputs.centroid_full_stockholm,
            outputs.centroid_match_stockholm,
            pid;
            full_label = "centroid_full",
            match_label = "centroid_match")
    end
    return outputs.s_exon_blocks_tsv
end
