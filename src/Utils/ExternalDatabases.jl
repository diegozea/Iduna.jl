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
