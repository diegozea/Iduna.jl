"""
    ensure_mmseqs_db(db) -> String

Check that an MMseqs2 database prefix has the files Iduna needs.

# Arguments

- `db::AbstractString`: MMseqs2 database prefix.

# Throws

- `ErrorException`: if the database looks incomplete.
"""
function ensure_mmseqs_db(db::AbstractString)
    for db_path in (String(db), string(db, "_aln"), string(db, "_seq"))
        dbtype = string(db_path, ".dbtype")
        isfile(dbtype) ||
            error("MMseqs2 database at $(db_path) looks incomplete; missing $(dbtype).")
    end
    return String(db)
end

# Immutable Database Snapshots
# ----------------------------

function _mmseqs_snapshot_files(db::AbstractString, manifest::AbstractString)
    root, prefix = dirname(abspath(db)), basename(db)
    return sort!([name
                  for name in readdir(root)
                  if (name == prefix || startswith(name, prefix * ".") ||
                      startswith(name, prefix * "_")) &&
                         isfile(joinpath(root, name)) &&
                         abspath(joinpath(root, name)) != abspath(manifest)])
end

function _require_mmseqs_components(files, prefix::AbstractString)
    for suffix in ("", "_seq", "_aln", "_h", "_seq_h")
        component = prefix * suffix
        for tail in (".dbtype", ".index")
            component * tail in files || error("missing component $(component * tail)")
        end
        # MMseqs split databases begin with part zero.
        component in files || component * ".0" in files ||
            error("missing data for $(component)")
    end
    return nothing
end

function _read_mmseqs_snapshot(db::AbstractString, manifest::AbstractString)
    bytes = read(manifest, String)
    value = JSON.parse(bytes)
    value isa AbstractDict || error("snapshot manifest must be an object")
    get(value, "schema_version", nothing) === 1 || error("unsupported snapshot schema")
    get(value, "immutable", nothing) === true ||
        error("snapshot must declare immutable=true")
    snapshot = get(value, "snapshot_id", nothing)
    snapshot isa AbstractString && !isempty(strip(snapshot)) || error("missing snapshot_id")
    components = get(value, "components", nothing)
    components isa AbstractDict || error("missing component checksums")
    files = _mmseqs_snapshot_files(db, manifest)
    _require_mmseqs_components(files, basename(db))
    Set(keys(components)) == Set(files) ||
        error("component inventory does not match database files")
    all(hash -> hash isa AbstractString && occursin(r"^[0-9a-f]{64}$", hash),
        values(components)) || error("components require SHA-256 checksums")
    # These checks bind a prepared immutable snapshot declaration, not mutable DB bytes.
    return (; snapshot_id = String(snapshot),
        manifest_sha256 = bytes2hex(SHA.sha256(bytes)), components)
end

function _mmseqs_snapshot_identity(db::AbstractString,
        manifest::Union{Nothing, AbstractString})
    if manifest === nothing
        @warn "MMseqs expansion reuse disabled: no immutable snapshot manifest was supplied."
        return nothing
    end
    try
        return _read_mmseqs_snapshot(db, manifest)
    catch err
        err isa InterruptException && rethrow()
        @warn "MMseqs expansion reuse disabled: invalid snapshot manifest." reason=sprint(showerror, err)
        return nothing
    end
end
