"""
    DEFAULT_PID_THRESHOLDS

Default ThorAxe percent identity (PID) thresholds tested during seed selection:
`[10.0, 20.0, 30.0, 60.0, 80.0]`.
"""
const DEFAULT_PID_THRESHOLDS = Float64[10, 20, 30, 60, 80]

"""
    S_EXON_CODE_FEATURE

Stockholm column annotation name used to store one short s-exon code per MSA
column. The value is `"SExonCode"`.
"""
const S_EXON_CODE_FEATURE = "SExonCode"

"""
    S_EXON_CODE_MAP_FEATURE

Stockholm file annotation name used to map short s-exon codes back to ThorAxe
s-exon IDs. The value is `"SExonCodeMap"`.
"""
const S_EXON_CODE_MAP_FEATURE = "SExonCodeMap"

# ThorAxe excludes '.' from PhyloSofS one-character s-exon codes and imputes
# missing s-exon IDs as "0_<number>", so this marker should not collide with
# ThorAxe-assigned s-exon provenance.
"""
    S_EXON_MISSING_CODE

Column code used when an MSA column has no assigned s-exon. The value is `'.'`.
"""
const S_EXON_MISSING_CODE = '.'

const _STAGE_STATE_SCHEMA_VERSION = 1
const _STAGE_STATE_FILE = "stage_state.json"
const _PROTEIN_ALIGNMENT_SCORE_MODEL = AffineGapScoreModel(
    BLOSUM62, gap_open = -10, gap_extend = -1)
const _TRANSIENT_HTTP_STATUSES = Set([429, 500, 502, 503, 504])

struct _RetryableHTTPStatus <: Exception
    response::HTTP.Response
end
