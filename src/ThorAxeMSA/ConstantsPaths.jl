# Constants and Paths
# -------------------

const _PHYLOSOFS_RESERVED_SYMBOLS = Set([
    ' ', '\t', '\n', '\r', '\v', '\f', '\\', '*', '>', '"', '\'', ',', '-', '_',
    '/', ';', '#', '$', '.', '&', '!', '@', '[', ']'
])
const _PHYLOSOFS_FALLBACK_SYMBOLS = [c
                                     for c in vcat(collect(Char(33):Char(126)),
                                             collect(Char(0x00BC):Char(0x017F)))
                                     if !(c in _PHYLOSOFS_RESERVED_SYMBOLS)]

const _REQUIRED_ENSEMBL_FILES = (
    "sequences.fasta",
    "exonstable.tsv",
    "tree.nh",
    "ensembl_version.csv",
    "tsl.csv"
)
const _ENSEMBL_REST_BASE = "https://rest.ensembl.org"
const _ENSEMBL_JSON_HEADERS = [
    "Accept" => "application/json",
    "Content-Type" => "application/json"
]
const _BIOMART_DATASETS_URL = "https://www.ensembl.org/biomart/martservice?type=datasets&mart=ENSEMBL_MART_ENSEMBL"
const _BIOMART_TEXT_HEADERS = ["Accept" => "text/plain"]
const _BIOMART_DATASETS_FILE = "ENSEMBL_MART_ENSEMBL_datasets.tsv"
const _BIOMART_DATASETS_METADATA_FILE = "ENSEMBL_MART_ENSEMBL_datasets.json"
const _TRANSCRIPT_QUERY_METADATA_FILE = "iduna_transcript_query.json"
const _SAMPLING_STRATEGIES = Set([:independent, :common, :input])
const _LOW_COMMON_SPECIES_THRESHOLD = 6
const _TRANSCRIPT_QUERY_SPINNER_INTERVAL_SECONDS = 1 / 3
const _ASES_DEFAULT_SPECIES = [
    "homo_sapiens",
    "gorilla_gorilla",
    "macaca_mulatta",
    "monodelphis_domestica",
    "rattus_norvegicus",
    "mus_musculus",
    "bos_taurus",
    "sus_scrofa",
    "ornithorhynchus_anatinus",
    "xenopus_tropicalis",
    "danio_rerio",
    "caenorhabditis_elegans"
]
const _ASES_DEFAULT_SPECIESLIST = join(_ASES_DEFAULT_SPECIES, ",")

_normalize_species_name(species::Nothing) = nothing
function _normalize_species_name(species::AbstractString)
    lowercase(replace(strip(String(species)), ' ' => '_'))
end

_thoraxe_input_dir(workdir::AbstractString) = joinpath(workdir, "thoraxe_input")
_thoraxe_msa_dir(workdir::AbstractString) = joinpath(workdir, "thoraxe_msa")
function _thoraxe_candidates_dir(workdir::AbstractString)
    joinpath(_thoraxe_msa_dir(workdir), "candidates")
end
_thoraxe_logs_dir(workdir::AbstractString) = joinpath(workdir, "logs", "thoraxe")
_thoraxe_pid_runs_dir(workdir::AbstractString) = joinpath(_thoraxe_msa_dir(workdir), "runs")
function _thoraxe_sample_species_dir(workdir::AbstractString)
    joinpath(_thoraxe_msa_dir(workdir),
        "samples", "species")
end
function _thoraxe_input_stage_dir(workdir::AbstractString)
    _pipeline_stage_dir(workdir, "thoraxe_input")
end
function _thoraxe_msa_stage_dir(workdir::AbstractString)
    _pipeline_stage_dir(workdir, "thoraxe_msa")
end
