# Stockholm Preparation
# ---------------------

function _is_s_exon_provenance_stockholm_line(line::AbstractString)
    return startswith(line, "#=GC SExonCode") ||
           startswith(line, "#=GF SExonCodeMap")
end

"""
    prepare_stockholm_for_mmseqs(source, dest) -> String

Copy a Stockholm MSA to a form that MMseqs2 can read.

Iduna removes its s-exon provenance annotations from this temporary copy because
MMseqs2 does not need them. The original file is not changed.

# Arguments

- `source::AbstractString`: source Stockholm MSA path.
- `dest::AbstractString`: destination path for the sanitized copy.
"""
function prepare_stockholm_for_mmseqs(source::AbstractString, dest::AbstractString)
    mkpath(dirname(dest))
    lines = readlines(source)
    has_header = !isempty(lines) && startswith(strip(lines[1]), "# STOCKHOLM")
    open(dest, "w") do io
        # MMseqs expects a complete Stockholm file with header and terminator.
        has_header || println(io, "# STOCKHOLM 1.0")
        for line in lines
            _is_s_exon_provenance_stockholm_line(line) && continue
            println(io, line)
        end
        if isempty(lines) || strip(lines[end]) != "//"
            println(io, "//")
        end
    end
    return dest
end

function _stockholm_annotation_data()
    # Keep each Stockholm record type separate so split lines can be joined safely.
    return (;
        comments = String[],
        gf_order = String[],
        gf_data = Dict{String, Vector{String}}(),
        gs_order = Tuple{String, String}[],
        gs_data = Dict{Tuple{String, String}, Vector{String}}(),
        gc_order = String[],
        gc_data = Dict{String, String}(),
        gr_order = Tuple{String, String}[],
        gr_data = Dict{Tuple{String, String}, String}(),
        seq_order = String[],
        seq_data = Dict{String, String}()
    )
end

function _append_stockholm_vector!(order, data, key, value)
    haskey(data, key) || (push!(order, key); data[key] = String[])
    push!(data[key], value)
    return nothing
end

function _append_stockholm_string!(order, data, key, value)
    if haskey(data, key)
        data[key] = string(data[key], value)
    else
        push!(order, key)
        data[key] = value
    end
    return nothing
end

function _record_gf_line!(records, line::AbstractString)
    parts = split(line; limit = 3)
    length(parts) < 3 && return nothing
    _append_stockholm_vector!(records.gf_order, records.gf_data, parts[2], parts[3])
    return nothing
end

function _record_gs_line!(records, line::AbstractString)
    parts = split(line; limit = 4)
    length(parts) < 4 && return nothing
    key = (parts[2], parts[3])
    _append_stockholm_vector!(records.gs_order, records.gs_data, key, parts[4])
    return nothing
end

function _record_gc_line!(records, line::AbstractString)
    parts = split(line; limit = 3)
    length(parts) < 3 && return nothing
    _append_stockholm_string!(
        records.gc_order, records.gc_data, parts[2], replace(parts[3], ' ' => ""))
    return nothing
end

function _record_gr_line!(records, line::AbstractString)
    parts = split(line; limit = 4)
    length(parts) < 4 && return nothing
    key = (parts[2], parts[3])
    _append_stockholm_string!(
        records.gr_order, records.gr_data, key, replace(parts[4], ' ' => ""))
    return nothing
end

function _record_sequence_line!(records, line::AbstractString)
    parts = split(line; limit = 2)
    length(parts) < 2 && return nothing
    name = strip(parts[1])
    fragment = replace(strip(parts[2]), ' ' => "")
    _append_stockholm_string!(records.seq_order, records.seq_data, name, fragment)
    return nothing
end

function _record_stockholm_line!(records, line::AbstractString)
    stripped = strip(line)
    isempty(stripped) && return nothing
    (startswith(stripped, "# STOCKHOLM") || startswith(stripped, "//")) &&
        return nothing
    startswith(line, "#=GF") && return _record_gf_line!(records, line)
    startswith(line, "#=GS") && return _record_gs_line!(records, line)
    startswith(line, "#=GC") && return _record_gc_line!(records, line)
    startswith(line, "#=GR") && return _record_gr_line!(records, line)
    startswith(line, '#') && (push!(records.comments, stripped); return nothing)
    return _record_sequence_line!(records, line)
end

function _write_stockholm_metadata(io, records)
    foreach(line -> println(io, line), records.comments)
    for feature in records.gf_order, value in records.gf_data[feature]

        println(io, "#=GF ", feature, ' ', value)
    end
    for key in records.gs_order, value in records.gs_data[key]

        println(io, "#=GS ", key[1], ' ', key[2], ' ', value)
    end
    for feature in records.gc_order
        println(io, "#=GC ", feature, ' ', records.gc_data[feature])
    end
    return nothing
end

function _write_stockholm_sequences(io, records)
    for name in records.seq_order
        println(io, name, '\t', records.seq_data[name])
        for key in records.gr_order
            key[1] == name || continue
            println(io, "#=GR ", key[1], ' ', key[2], ' ', records.gr_data[key])
        end
    end
    return nothing
end

function _write_normalized_stockholm(path::AbstractString, records)
    open(path, "w") do io
        println(io, "# STOCKHOLM 1.0")
        _write_stockholm_metadata(io, records)
        _write_stockholm_sequences(io, records)
        println(io, "//")
    end
    return path
end

"""
    normalize_stockholm_annotations!(path) -> String

Rewrite a Stockholm file so split annotations and sequence fragments are stored
in a stable order.

# Arguments

- `path::AbstractString`: Stockholm file to rewrite in place.
"""
function normalize_stockholm_annotations!(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return path

    records = _stockholm_annotation_data()
    # Merge split Stockholm records while keeping the original output order.
    for line in lines
        _record_stockholm_line!(records, line)
    end
    return _write_normalized_stockholm(path, records)
end
