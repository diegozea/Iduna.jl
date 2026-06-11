# Species Lists
# -------------

function _normalized_specieslist(specieslist::Union{Nothing, AbstractString})
    specieslist === nothing && return nothing
    stripped = strip(String(specieslist))
    return isempty(stripped) ? nothing : stripped
end

function _resolve_specieslist_preset(specieslist::AbstractString)
    stripped = strip(String(specieslist))
    if stripped == "ases"
        return (;
            specieslist = _ASES_DEFAULT_SPECIESLIST,
            mode = :ases,
            input = stripped)
    elseif isempty(stripped) || stripped == "all"
        return (;
            specieslist = nothing,
            mode = :all,
            input = stripped)
    end
    return (;
        specieslist = stripped,
        mode = :explicit,
        input = stripped)
end

function _log_specieslist_choice(resolved)
    specieslist = resolved.specieslist
    if resolved.mode === :ases
        @info "Using Ases default ThorAxe species list." n_species=length(_ASES_DEFAULT_SPECIES) specieslist
    elseif resolved.mode === :all
        @info "Using unrestricted ThorAxe species selection." specieslist=resolved.input
    elseif occursin(',', specieslist)
        summary = _specieslist_log_summary(specieslist)
        @info "Using explicit ThorAxe species list." n_species=summary.n_species specieslist_preview=summary.specieslist_preview
    elseif isfile(specieslist)
        summary = _specieslist_log_summary(specieslist)
        @info "Using ThorAxe species list file." specieslist_path=specieslist n_species=summary.n_species specieslist_preview=summary.specieslist_preview
    else
        species = only(_parse_specieslist(specieslist))
        @info "Using explicit ThorAxe species." species
    end
    return nothing
end

function _orthology_relationships(orthology::AbstractString)
    value = strip(String(orthology))
    if value == "1:1"
        return ["ortholog_one2one"]
    elseif value == "1:n"
        return ["ortholog_one2one", "ortholog_one2many"]
    elseif value == "m:n"
        return ["ortholog_one2one", "ortholog_one2many", "ortholog_many2many"]
    end
    error("Orthology must be one of 1:1, 1:n, or m:n; got $(orthology).")
end

function _unique_nonempty_species(items)
    seen = Set{String}()
    species = String[]
    for item in items
        item === nothing && continue
        value = String(item)
        isempty(value) && continue
        value in seen && continue
        push!(seen, value)
        push!(species, value)
    end
    return species
end

function _parse_specieslist(specieslist::Union{Nothing, AbstractString})
    normalized = _normalized_specieslist(specieslist)
    normalized === nothing && return nothing

    raw_species = String[]
    # The same option accepts a comma list, a file path, or one species name.
    if occursin(',', normalized)
        append!(raw_species, split(normalized, ','))
    elseif isfile(normalized)
        append!(raw_species, readlines(normalized))
    else
        push!(raw_species, normalized)
    end

    return _unique_nonempty_species(
        _normalize_species_name(species) for species in raw_species)
end

function _specieslist_string(species::AbstractVector{<:AbstractString})
    isempty(species) && return nothing
    return join(species, ",")
end

function _specieslist_log_summary(specieslist::Union{Nothing, AbstractString})
    species = _parse_specieslist(specieslist)
    species === nothing && return (; n_species = nothing, specieslist_preview = nothing)
    preview_count = min(length(species), 5)
    preview = isempty(species) ? "" : join(species[1:preview_count], ",")
    length(species) > preview_count && (preview = string(preview, ",..."))
    return (; n_species = length(species), specieslist_preview = preview)
end

function _specieslist_log_details(specieslist::Union{Nothing, AbstractString})
    summary = _specieslist_log_summary(specieslist)
    summary.n_species === nothing && return (;)
    return (; n_species = summary.n_species,
        specieslist_preview = summary.specieslist_preview)
end

function _pid_thresholds_preview(pid_thresholds::AbstractVector{<:Real})
    pid_values = Float64.(pid_thresholds)
    preview_count = min(length(pid_values), 5)
    preview = isempty(pid_values) ? "" : join(pid_values[1:preview_count], ",")
    length(pid_values) > preview_count && (preview = string(preview, ",..."))
    return preview
end

# BioMart Species Filtering
# -------------------------

function _fetch_biomart_datasets_text(;
        url::AbstractString = _BIOMART_DATASETS_URL,
        retries::Integer = 4,
        http_get::Function = _http_get_request)
    # BioMart can return temporary 5xx/429 responses, so fetch through retry.
    resp = _http_get_with_retries(
        url, _BIOMART_TEXT_HEADERS; retries, sleep_seconds = 1.0, http_get)
    resp.status == 200 && return decode_body(resp)
    error("BioMart datasets metadata request failed with HTTP status $(resp.status).")
end

function _parse_biomart_gene_datasets(text::AbstractString)
    datasets = Set{String}()
    for line in split(String(text), '\n')
        fields = split(strip(line), '\t')
        length(fields) >= 2 || continue
        fields[1] == "TableSet" || continue
        dataset = strip(fields[2])
        endswith(dataset, "_gene_ensembl") || continue
        push!(datasets, dataset)
    end
    return datasets
end

function _read_biomart_cache_date(metadata_path::AbstractString)
    isfile(metadata_path) || return nothing
    try
        metadata = JSON.parse(read(metadata_path, String))
        date = get(metadata, "download_date", nothing)
        date === nothing && return nothing
        return String(date)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _write_biomart_datasets_cache!(cache_dir::AbstractString,
        text::AbstractString;
        today::Dates.Date = Dates.today(),
        url::AbstractString = _BIOMART_DATASETS_URL)
    mkpath(cache_dir)
    datasets_path = joinpath(cache_dir, _BIOMART_DATASETS_FILE)
    metadata_path = joinpath(cache_dir, _BIOMART_DATASETS_METADATA_FILE)
    tmp_datasets = string(datasets_path, ".tmp")
    write_text(tmp_datasets, text)
    mv(tmp_datasets, datasets_path; force = true)
    metadata = (;
        download_date = string(today),
        download_time = string(Dates.now()),
        url,
        status = 200
    )
    tmp_metadata = string(metadata_path, ".tmp")
    write_json(tmp_metadata, metadata)
    mv(tmp_metadata, metadata_path; force = true)
    return datasets_path
end

function _read_cached_biomart_datasets(cache_dir::AbstractString)
    datasets_path = joinpath(cache_dir, _BIOMART_DATASETS_FILE)
    isfile(datasets_path) || return nothing
    datasets = _parse_biomart_gene_datasets(read(datasets_path, String))
    isempty(datasets) && return nothing
    return datasets
end

function _load_biomart_gene_datasets(;
        cache_dir::AbstractString = _biomart_cache_dir(),
        today::Dates.Date = Dates.today(),
        fetcher::Function = _fetch_biomart_datasets_text)
    datasets_path = joinpath(cache_dir, _BIOMART_DATASETS_FILE)
    metadata_path = joinpath(cache_dir, _BIOMART_DATASETS_METADATA_FILE)
    cached_date = _read_biomart_cache_date(metadata_path)
    cached_is_current = isfile(datasets_path) && cached_date == string(today)
    # Use one fresh BioMart dataset list per day to avoid repeated network hits.
    if cached_is_current
        cached = _read_cached_biomart_datasets(cache_dir)
        cached !== nothing && return (datasets = cached, warnings = String[])
    end

    try
        text = fetcher()
        datasets = _parse_biomart_gene_datasets(text)
        isempty(datasets) &&
            error("BioMart datasets metadata did not contain Ensembl Gene datasets.")
        _write_biomart_datasets_cache!(cache_dir, text; today)
        return (datasets = datasets, warnings = String[])
    catch err
        err isa InterruptException && rethrow()
        cached = _read_cached_biomart_datasets(cache_dir)
        if cached !== nothing
            stale_date = cached_date === nothing ? "unknown date" : cached_date
            warning = "BioMart datasets metadata refresh failed; using stale cache from $(stale_date). $(sprint(showerror, err))"
            return (datasets = cached, warnings = [warning])
        end
        warning = "BioMart datasets metadata refresh failed; using the unfiltered specieslist. $(sprint(showerror, err))"
        return (datasets = nothing, warnings = [warning])
    end
end

function _biomart_gene_dataset_for_species(species::AbstractString)
    normalized = _normalize_species_name(species)
    normalized === nothing && return nothing
    parts = split(normalized, '_')
    length(parts) == 2 || return nothing
    isempty(parts[1]) && return nothing
    isempty(parts[2]) && return nothing
    return string(first(parts[1]), parts[2], "_gene_ensembl")
end

function _classify_biomart_species!(
        kept::Vector{String},
        removed::Vector{String},
        unchecked::Vector{String},
        species::AbstractString,
        datasets,
        query_species::Union{Nothing, String})
    dataset = _biomart_gene_dataset_for_species(species)
    if dataset === nothing
        push!(kept, species)
        push!(unchecked, species)
        return false
    elseif dataset in datasets
        push!(kept, species)
        return false
    elseif species == query_species
        push!(kept, species)
        return true
    end
    push!(removed, species)
    return false
end

function _resolve_biomart_datasets_specieslist(target::ResolvedTarget,
        specieslist::Union{Nothing, AbstractString};
        dataset_loader::Function = _load_biomart_gene_datasets)
    requested_species = _parse_specieslist(specieslist)
    requested_species === nothing &&
        return (specieslist = _normalized_specieslist(specieslist), warnings = String[])

    loaded = dataset_loader()
    warnings = String.(loaded.warnings)
    datasets = loaded.datasets
    datasets === nothing &&
        return (specieslist = _normalized_specieslist(specieslist), warnings = warnings)

    query_species = _normalize_species_name(target.species)
    kept = String[]
    removed = String[]
    unchecked = String[]
    missing_query_dataset = false
    # Drop species that BioMart cannot serve, but keep aliases we cannot prove.
    for species in requested_species
        missing_query_dataset |= _classify_biomart_species!(
            kept, removed, unchecked, species, datasets, query_species)
    end

    if !isempty(removed)
        push!(warnings,
            "BioMart datasets filter removed species without matching Ensembl Gene datasets or recognized aliases: $(join(removed, ", ")).")
    end
    if !isempty(unchecked)
        push!(warnings,
            "BioMart datasets filter could not derive Ensembl Gene dataset names for possible species aliases: $(join(unchecked, ", ")); keeping them unchanged.")
    end
    if missing_query_dataset && query_species !== nothing
        push!(warnings,
            "BioMart datasets filter did not find an Ensembl Gene dataset for query species $(query_species); keeping it because transcript_query requires the query species.")
    end
    isempty(kept) && error("BioMart datasets filter removed all requested species.")
    return (specieslist = _specieslist_string(kept), warnings = warnings)
end

# Ensembl Homology Filtering
# --------------------------

function _homology_species(data, orthology::AbstractString)
    wanted = Set(_orthology_relationships(orthology))
    species = String[]
    for item in get(data, "data", Any[])
        for homology in get(item, "homologies", Any[])
            type = get(homology, "type", nothing)
            type isa AbstractString && String(type) in wanted || continue
            target = get(homology, "target", nothing)
            target === nothing && continue
            target_species = get(target, "species", nothing)
            target_species isa AbstractString || continue
            push!(species, _normalize_species_name(target_species))
        end
    end
    return _unique_nonempty_species(species)
end

function _fetch_ensembl_homology_data(species::AbstractString,
        gene_id::AbstractString;
        retries::Integer = 4,
        sleep_seconds::Real = 1.5,
        http_get::Function = _http_get_request)
    gene_core = strip_ensembl_version(gene_id)
    url = "$(_ENSEMBL_REST_BASE)/homology/id/$(species)/$(gene_core)?type=orthologues;sequence=none"
    resp = _http_get_with_retries(
        url, _ENSEMBL_JSON_HEADERS; retries, sleep_seconds, http_get)
    resp.status == 200 && return JSON.parse(decode_body(resp))
    error("Ensembl homology specieslist filter failed for $(gene_core) in $(species) with HTTP status $(resp.status).")
end

function _fetch_ortholog_species(target::ResolvedTarget, orthology::AbstractString;
        homology_data_fetcher::Function = _fetch_ensembl_homology_data)
    species = _normalize_species_name(target.species)
    species === nothing &&
        error("Cannot run Ensembl specieslist filter because the target species is unknown.")
    data = homology_data_fetcher(species, target.ensembl_gene_id)
    return _homology_species(data, orthology)
end

function _prepend_query_species(species::AbstractVector{<:AbstractString},
        query_species::Union{Nothing, AbstractString})
    query_species === nothing && return String.(species)
    return _unique_nonempty_species(Iterators.flatten(([String(query_species)], species)))
end

function _resolve_effective_specieslist(target::ResolvedTarget,
        specieslist::Union{Nothing, AbstractString},
        orthology::AbstractString;
        homology_species_fetcher::Function = _fetch_ortholog_species)
    _orthology_relationships(orthology)
    query_species = _normalize_species_name(target.species)
    # If Ensembl homology is unavailable, fall back to the user's species list.
    ortholog_species = try
        homology_species_fetcher(target, orthology)
    catch err
        err isa InterruptException && rethrow()
        warning = "Ensembl specieslist filter failed; using the unfiltered specieslist. $(sprint(showerror, err))"
        return (specieslist = _normalized_specieslist(specieslist), warnings = [warning])
    end

    ortholog_species = _unique_nonempty_species(ortholog_species)
    if isempty(ortholog_species)
        error("Ensembl specieslist filter found no $(orthology) ortholog species for $(target.ensembl_gene_id).")
    end

    requested_species = _parse_specieslist(specieslist)
    if requested_species === nothing
        effective_species = _prepend_query_species(ortholog_species, query_species)
        return (specieslist = _specieslist_string(effective_species), warnings = String[])
    end

    allowed_targets = Set(ortholog_species)
    matching_targets = [species
                        for species in requested_species if species in allowed_targets]
    if isempty(matching_targets)
        error("Ensembl specieslist filter found no requested species with $(orthology) orthologs for $(target.ensembl_gene_id).")
    end

    removed = [species
               for species in requested_species
               if species != query_species && !(species in allowed_targets)]
    warnings = String[]
    if !isempty(removed)
        removed_text = join(removed, ", ")
        push!(warnings,
            "Removed species without $(orthology) Ensembl orthologs for $(target.ensembl_gene_id): $(removed_text).")
    end

    effective_species = _prepend_query_species(matching_targets, query_species)
    return (specieslist = _specieslist_string(effective_species), warnings = warnings)
end
