"""
    s_exon_blocks_path(stockholm_path) -> String

Return the default TSV path for s-exon column blocks next to a Stockholm MSA.

# Arguments

- `stockholm_path::AbstractString`: Stockholm MSA path.
"""
function s_exon_blocks_path(stockholm_path::AbstractString)
    root, _ = splitext(String(stockholm_path))
    return string(root, "_s_exon_blocks.tsv")
end

"""
    s_exon_codes(msa; feature=S_EXON_CODE_FEATURE)

Read the per-column s-exon codes from an MSA.

When the annotation is absent, this returns one [`S_EXON_MISSING_CODE`](@ref)
for each MSA column.

# Arguments

- `msa`: MSA to inspect.

# Keywords

- `feature::AbstractString = S_EXON_CODE_FEATURE`: Stockholm column annotation
  name used to read s-exon codes.

# Throws

- `ErrorException`: if the requested annotation exists but has a different length
  from the MSA.
"""
function s_exon_codes(msa::AbstractMultipleSequenceAlignment;
        feature::AbstractString = S_EXON_CODE_FEATURE)
    codes = getannotcolumn(msa, String(feature), "")
    if isempty(codes)
        return repeat(string(S_EXON_MISSING_CODE), ncolumns(msa))
    end
    length(codes) == ncolumns(msa) ||
        error("$(feature) has $(length(codes)) characters, but the MSA has $(ncolumns(msa)) columns.")
    return String(codes)
end

"""
    has_s_exon_annotations(msa; feature=S_EXON_CODE_FEATURE) -> Bool

Return `true` when an MSA has Iduna s-exon column annotations.

# Arguments

- `msa`: MSA to inspect.

# Keywords

- `feature::AbstractString = S_EXON_CODE_FEATURE`: Stockholm column annotation
  name used to detect s-exon codes.
"""
function has_s_exon_annotations(msa::AbstractMultipleSequenceAlignment;
        feature::AbstractString = S_EXON_CODE_FEATURE)
    codes = getannotcolumn(msa, String(feature), "")
    return !isempty(codes) && length(codes) == ncolumns(msa)
end

"""
    s_exon_code_map(msa; feature=S_EXON_CODE_MAP_FEATURE)

Read the map from short s-exon codes to ThorAxe s-exon IDs.

When the annotation is absent, this returns an empty dictionary.

# Arguments

- `msa`: MSA to inspect.

# Keywords

- `feature::AbstractString = S_EXON_CODE_MAP_FEATURE`: Stockholm file annotation
  name used to read the s-exon code map.
"""
function s_exon_code_map(msa::AbstractMultipleSequenceAlignment;
        feature::AbstractString = S_EXON_CODE_MAP_FEATURE)
    raw = getannotfile(msa, String(feature), "")
    map = Dict{Char, String}()
    isempty(raw) && return map
    pos = firstindex(raw)
    # The map is stored as compact quoted pairs, for example "A"=>"exon_1".
    for match in eachmatch(r"(\"(?:\\.|[^\"])*\")=>(\"(?:\\.|[^\"])*\")", raw)
        match_start = match.offset
        separator = pos == match_start ? "" : raw[pos:prevind(raw, match_start)]
        all(==(','), separator) ||
            error("Invalid $(feature) entry $(repr(raw)); expected quoted code=>s_exon_id pairs.")
        code = only(String(Meta.parse(match.captures[1])))
        code == S_EXON_MISSING_CODE &&
            error("$(S_EXON_MISSING_CODE) is reserved for columns without s-exon provenance.")
        map[code] = String(Meta.parse(match.captures[2]))
        pos = nextind(raw, match.offset + ncodeunits(match.match) - 1)
    end
    pos > lastindex(raw) ||
        error("Invalid $(feature) entry $(repr(raw)); expected quoted code=>s_exon_id pairs.")
    return map
end

function _format_s_exon_code_map(code_map)
    entries = ["$(repr(string(code)))=>$(repr(String(id)))" for (code, id) in code_map]
    return join(entries, ',')
end

"""
    set_s_exon_annotations!(msa, codes, code_map;
                            code_feature=S_EXON_CODE_FEATURE,
                            map_feature=S_EXON_CODE_MAP_FEATURE)

Store Iduna s-exon annotations on an MSA.

# Arguments

- `msa`: MSA to annotate.
- `codes::AbstractString`: one short code for each MSA column.
- `code_map`: map from each short code to its ThorAxe s-exon ID.

# Keywords

- `code_feature::AbstractString = S_EXON_CODE_FEATURE`: Stockholm column
  annotation name used to store s-exon codes.
- `map_feature::AbstractString = S_EXON_CODE_MAP_FEATURE`: Stockholm file
  annotation name used to store the s-exon code map.
"""
function set_s_exon_annotations!(msa::AbstractMultipleSequenceAlignment,
        codes::AbstractString,
        code_map;
        code_feature::AbstractString = S_EXON_CODE_FEATURE,
        map_feature::AbstractString = S_EXON_CODE_MAP_FEATURE)
    length(codes) == ncolumns(msa) ||
        error("$(code_feature) has $(length(codes)) characters, but the MSA has $(ncolumns(msa)) columns.")
    setannotcolumn!(msa, String(code_feature), String(codes))
    setannotfile!(msa, String(map_feature), _format_s_exon_code_map(code_map))
    return msa
end

function _write_s_exon_block_line(io,
        alignment::AbstractString,
        pid_value::AbstractString,
        code::Char,
        s_exon_id::AbstractString,
        start_col::Integer,
        end_col::Integer)
    println(io,
        alignment, '\t',
        pid_value, '\t',
        code, '\t',
        s_exon_id, '\t',
        start_col, '\t',
        end_col, '\t',
        end_col - start_col + 1)
    return nothing
end

function _write_s_exon_block_rows(io,
        codes::AbstractString,
        code_map,
        alignment::AbstractString,
        pid_value::AbstractString)
    start_col = 0
    current = S_EXON_MISSING_CODE
    for (idx, code) in enumerate(codes)
        if code != current
            if current != S_EXON_MISSING_CODE
                # Finish the previous run of neighboring columns with the same s-exon.
                _write_s_exon_block_line(io, alignment, pid_value, current,
                    get(code_map, current, ""), start_col, idx - 1)
            end
            current = code
            start_col = idx
        end
    end
    current == S_EXON_MISSING_CODE && return nothing
    _write_s_exon_block_line(io, alignment, pid_value, current,
        get(code_map, current, ""), start_col, length(codes))
    return nothing
end

"""
    write_s_exon_blocks_tsv(path, msa; alignment, pid=nothing, append=false,
                            code_feature=S_EXON_CODE_FEATURE,
                            map_feature=S_EXON_CODE_MAP_FEATURE)

Write a TSV table that groups neighboring MSA columns with the same s-exon.

# Arguments

- `path::AbstractString`: TSV output path.
- `msa`: MSA to summarize. When s-exon annotations are absent, only the header is
  written.

# Keywords

- `alignment::AbstractString`: label for the MSA being written.
- `pid = nothing`: percent identity threshold to record in the table. When
  `nothing`, the PID column is left empty for those rows.
- `append::Bool = false`: append rows to an existing table.
- `code_feature::AbstractString = S_EXON_CODE_FEATURE`: Stockholm column
  annotation name used to read s-exon codes.
- `map_feature::AbstractString = S_EXON_CODE_MAP_FEATURE`: Stockholm file
  annotation name used to read the s-exon code map.
"""
function write_s_exon_blocks_tsv(path::AbstractString,
        msa::AbstractMultipleSequenceAlignment;
        alignment::AbstractString,
        pid = nothing,
        append::Bool = false,
        code_feature::AbstractString = S_EXON_CODE_FEATURE,
        map_feature::AbstractString = S_EXON_CODE_MAP_FEATURE)
    codes = s_exon_codes(msa; feature = code_feature)
    code_map = s_exon_code_map(msa; feature = map_feature)
    mkpath(dirname(path))
    pid_value = pid === nothing ? "" : string(pid)
    open(path, append ? "a" : "w") do io
        append ||
            println(io, "alignment\tpid\tcode\ts_exon_id\tstart_col\tend_col\tn_columns")
        _write_s_exon_block_rows(io, codes, code_map, alignment, pid_value)
    end
    return path
end
