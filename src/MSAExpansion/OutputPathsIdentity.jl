# Output Paths and Identity
# -------------------------

_normalize_id(s::AbstractString) = String(split(String(s))[1])

const _STEP_STATE_FILE = "stage_state.json"

function _expansion_root(workdir::AbstractString)
    joinpath(workdir, "expansion")
end

function _expansion_output_paths(run_dir::AbstractString,
        transcript_id::AbstractString;
        centroids::Bool = false)
    expanded_dir = joinpath(run_dir, "expanded_msa")
    outputs = (;
        full_stockholm = joinpath(expanded_dir, "$(transcript_id)_full.sto"),
        match_stockholm = joinpath(expanded_dir, "$(transcript_id)_matchonly.sto"),
        a3m_path = joinpath(expanded_dir, "$(transcript_id)_expanded.a3m"),
        s_exon_blocks_tsv = joinpath(expanded_dir, "$(transcript_id)_s_exon_blocks.tsv"),
        hits_fasta = joinpath(expanded_dir, "$(transcript_id)_hits_raw.fasta")
    )
    if !centroids
        return outputs
    end
    centroid_dir = joinpath(run_dir, "centroid_msa")
    return merge(outputs,
        (;
            centroid_full_stockholm = joinpath(
                centroid_dir, "$(transcript_id)_centroids_full.sto"),
            centroid_match_stockholm = joinpath(
                centroid_dir, "$(transcript_id)_centroids_matchonly.sto"),
            centroid_a3m_path = joinpath(centroid_dir, "$(transcript_id)_centroids.a3m"),
            centroid_s_exon_blocks_tsv = joinpath(
                centroid_dir, "$(transcript_id)_centroids_s_exon_blocks.tsv"),
            centroid_hits_fasta = joinpath(
                centroid_dir, "$(transcript_id)_centroid_hits_raw.fasta")
        ))
end

_step_state_path(run_dir::AbstractString) = _stage_state_path(run_dir)

function _required_expansion_outputs(outputs::NamedTuple)
    # The s-exon block tables can be rebuilt from the Stockholm alignments.
    return Dict(String(name) => path
    for (name, path) in pairs(outputs)
    if !endswith(String(name), "s_exon_blocks_tsv"))
end

function _expansion_identity(target::ResolvedTarget,
        seed::SeedSelection,
        seed_stockholm::AbstractString,
        seed_fasta::Union{Nothing, AbstractString},
        mmseqs_db::AbstractString;
        match_mode::Integer,
        match_ratio::Union{Nothing, Real},
        hmmbuild_symfrac::Real,
        centroids::Bool)
    return (;
        target = (;
            input_id = target.input_id,
            input_kind = String(target.input_kind),
            uniprot_id = target.uniprot_id,
            ensembl_gene_id = target.ensembl_gene_id,
            transcript_id = target.transcript_id,
            ensembl_protein_id = target.ensembl_protein_id,
            species = target.species
        ),
        seed = (;
            pid = Float64(seed.pid),
            stockholm_sha256 = _file_sha256(seed_stockholm),
            fasta_sha256 = seed_fasta === nothing || !isfile(seed_fasta) ? nothing :
                           _file_sha256(seed_fasta)
        ),
        expansion = (;
            mmseqs_db = abspath(String(mmseqs_db)),
            match_mode = Int(match_mode),
            match_ratio = match_ratio === nothing ? nothing : Float64(match_ratio),
            hmmbuild_symfrac = Float64(hmmbuild_symfrac),
            centroids = Bool(centroids)
        )
    )
end
