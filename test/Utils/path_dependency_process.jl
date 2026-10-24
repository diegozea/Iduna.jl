using Iduna, JSON

const root, phase = ARGS
const utils = Iduna.Utils
target = Iduna.ResolvedTarget(; input_id = "P20963", input_kind = :uniprot,
    uniprot_id = "P20963", ensembl_gene_id = "seq1", transcript_id = "ENST000001",
    uniprot_sequence_path = joinpath(root, "uniprot.fasta"))
seed = Iduna.SeedSelection(; pid = 10.0, epli = 1.0,
    stockholm_path = joinpath(root, "seed.sto"),
    summary_path = joinpath(root, "seed_summary.csv"))
identity = utils._implementation_identity()
@assert identity !== nothing
dependency = only(filter(record -> record.name == "BioAlignments", identity.dependencies))
@assert dependency.source_sha256 !== nothing
cached = Iduna.ResultsValidation.validate_results(target, seed, nothing,
    joinpath(root, "shared"))
state = utils._read_stage_state(dirname(cached.stats_path))
fresh = Iduna.ResultsValidation.validate_results(target, seed, nothing,
    joinpath(root, "fresh_" * phase))
function result_fields(value)
    (; status = value.status, identical = value.aln_identical,
        insertions = value.aln_insertions, deletions = value.aln_deletions)
end
result = (; source_sha256 = dependency.source_sha256,
    implementation_hash = utils._identity_hash(identity),
    package_path = pkgdir(utils.BioAlignments),
    blosum62_aa = utils.BLOSUM62[utils.BioSequences.AA_A, utils.BioSequences.AA_A],
    action = state["action"], cached = result_fields(cached), fresh = result_fields(fresh))
write(joinpath(root, phase * ".json"), JSON.json(result))
println(JSON.json(result))
