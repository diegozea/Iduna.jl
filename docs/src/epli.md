# EPLI Score

EPLI means Estimated Profile-Level Identities. It measures how stable an aligner
is when the input sequence set is subsampled.

## Quick Example

This example creates four fake protein sequences with MIToS, aligns them with
ProGraphMSA, and computes EPLI with the default sampling settings.

```jldoctest; output = false
using Iduna
using MIToS.MSA

sequences = [
    AnnotatedSequence("ref", "MKTAYIAKQRQISFVKSHFSRQDILD"),
    AnnotatedSequence("seq2", "MKTAYIAKQRQISFVKSHFSRQNILD"),
    AnnotatedSequence("seq3", "MKTAYIAKQKQISFVKSHFSRQDILD"),
    AnnotatedSequence("seq4", "MKTAYIAKQRQISFVKSHFTRQDILD"),
]

workdir = mktempdir()

result = Iduna.EPLI.epli_score(
    sequences,
    workdir,
    Iduna.EPLI.prographmsa_aligner,
)

(result.n_samples, result.input_source, isfile(result.summary_path))

# output

(45, "mitos_sequences", true)
```

The call writes a full input FASTA under `workdir`, runs ProGraphMSA for the full
sequence set and for each sampled subset, compares each sampled MSA to the full
MSA with HHsuite, and stores the median normalized score in `result.epli`.

The main output files are:

- `result.summary_path`: one-row CSV with the run settings and final EPLI.
- `result.scores_path`: per-sample scores used to compute the median.
- `result.reference_msa_fasta`: the full ProGraphMSA alignment.

## Input Types

`epli_score` accepts these sequence sources:

- An unaligned FASTA path.
- A `Vector` of MIToS `AbstractSequence` objects, such as `AnnotatedSequence`.
- A single MIToS `AbstractSequence`.
- A `MIToS.MSA.AbstractMultipleSequenceAlignment`.

MIToS sequence and MSA inputs are converted to unaligned FASTA before the aligner
runs. Gap characters `-` and `.` are removed from each sequence.

For aligned MSA input:

```julia
using Iduna
using MIToS.MSA

msa = read_file("aligned_sequences.fasta", FASTA)

result = Iduna.EPLI.epli_score(
    msa,
    "msa_epli",
    Iduna.EPLI.prographmsa_aligner;
    reference_sequence = "seq1",
)
```

The returned result and `summary.csv` use `input_fasta = missing` for MIToS
object inputs.

## Wrap An Aligner

An aligner wrapper is a Julia function with this required contract:

```julia
function my_aligner(input_fasta, output_fasta; logs_dir = nothing,
        run_label = "run", aligner_args = Cmd(String[]))
    run(`my-aligner --input $input_fasta --output $output_fasta $aligner_args`)
    return output_fasta
end
```

The wrapper receives an unaligned FASTA and must write an aligned FASTA. It can
return the output path, return `nothing`, or return a named tuple with a
`fasta_path` field.

EPLI sets these arguments when it calls the wrapper:

- `input_fasta`: unaligned FASTA for this run. It can contain the full input set
  or one sampled subset.
- `output_fasta`: path where EPLI expects the aligned FASTA for this run.
- `logs_dir`: optional directory for logs.
- `run_label`: stable label for the current run, such as `full` or
  `sequence_subset_001`.
- `aligner_args`: extra command arguments supplied by the `epli_score` caller.

`aligner_args` is intentionally a `Cmd`, not open-ended Julia keywords. This
keeps EPLI options separate from command-line options for the external aligner.
Wrappers must accept `run_label` and `aligner_args`, even when callers use the
default empty command.

## Optional JLL Aligners

Iduna also provides unexported wrappers for FAMSA, Kalign, MAFFT, Clustal
Omega, and MUSCLE through Julia package extensions. These wrappers are
available only after loading the matching JLL package.

All bundled aligner wrappers return aligned FASTA records in the same order as
the input FASTA. This keeps EPLI's reference row stable across aligners.
Wrappers that reorder aligner output use MIToS FASTA parsing and writing, so
input and aligned output sequences should use the MIToS residue alphabet.

For MAFFT:

```julia
using Iduna
using MAFFT_jll

result = Iduna.EPLI.epli_score(
    "sequences.fasta",
    "mafft_epli",
    Iduna.EPLI.mafft_aligner,
    aligner_args = `--thread 4 --localpair --maxiterate 1000`,
)
```

For FAMSA:

```julia
using Iduna
using FAMSA_jll

result = Iduna.EPLI.epli_score(
    "sequences.fasta",
    "famsa_epli",
    Iduna.EPLI.famsa_aligner,
    aligner_args = `-t 4`,
)
```

For Kalign:

```julia
using Iduna
using kalign_jll

result = Iduna.EPLI.epli_score(
    "sequences.fasta",
    "kalign_epli",
    Iduna.EPLI.kalign_aligner,
    aligner_args = `-n 4`,
)
```

The Kalign wrapper fixes FASTA output with `--format fasta`. It does not force
`--type protein`, because Kalign rejects DNA-like protein alphabets such as
`AAAA` before alignment when that flag is set. Pass `--type` in `aligner_args`
only when you want Kalign to validate a specific sequence type.

For Clustal Omega:

```julia
using Iduna
using ClustalO_jll

result = Iduna.EPLI.epli_score(
    "sequences.fasta",
    "clustalo_epli",
    Iduna.EPLI.clustalo_aligner,
    aligner_args = `--threads=4 --auto`,
)
```

For MUSCLE:

```julia
using Iduna
using MUSCLE_jll

result = Iduna.EPLI.epli_score(
    "sequences.fasta",
    "muscle_epli",
    Iduna.EPLI.muscle_aligner,
)
```

Use the qualified names `Iduna.EPLI.famsa_aligner`,
`Iduna.EPLI.kalign_aligner`, `Iduna.EPLI.mafft_aligner`,
`Iduna.EPLI.clustalo_aligner`, and `Iduna.EPLI.muscle_aligner`; they are
intentionally not exported.

## Alternative Scores

By default, EPLI uses HHsuite identity counts and comparable-position
normalization:

```text
100 * matched_positions / comparable_positions
```

For profile-score experiments, pass a different score and normalization function:

```julia
using Iduna
using Iduna.EPLI

result = epli_score(
    "sequences.fasta",
    "epli_work",
    Iduna.EPLI.prographmsa_aligner;
    score_fn = hhsuite_profile_score,
    normalization_fn = self_reference_normalization,
)
```

This computes each sample's raw HHalign profile score and normalizes it by the
raw score from aligning the full reference MSA against itself.

### Score values and callback contracts

The built-in score functions return immutable domain values: `IdentityScore`
contains `IdentityCounts`, while `RawScore` holds a finite profile or custom score.
`IdentityCounts` requires integer counts with `0 ≤ matched ≤ comparable`. Its raw
identity score is derived from the matched count. Use `raw_score(score)` to read
raw values. Optional tool information lives in `score.diagnostics`; unavailable
or non-finite diagnostics do not invalidate an unrelated scientific score.

All normalizers return `NormalizationResult`. Its constructor converts to Float64
and rejects non-finite components and non-finite supplied denominators. Negative
raw scores and normalized ratios above 100 remain valid. Self-reference normalization
requires a nonzero reference score; no normalization uses a missing denominator.

Custom score callbacks should return `RawScore` or `IdentityScore`. Named tuples
with `raw_score` are accepted at the input boundary and converted to `RawScore`;
complete, nonmissing identity counts are validated there as well. Normalization
callbacks receive this validated object, not the original tuple. Access additional
callback fields through `score.diagnostics`. A normalization callback may return
`NormalizationResult`, a real number, or a named tuple with `normalized_score` and
an optional `normalization_score`; the latter two are converted at the boundary.
Callable objects are supported. Extend `requires_reference_score` for a normalizer
that needs a reference self score, or set `compute_reference_self_score` explicitly.

The final median is an `EPLIScore`, whose constructor requires a finite Float64
value. Both `result.epli` and an available `seed.epli` expose that value as `.value`.
A seed without a score has `seed.epli === nothing`. Its optional `epli_reason`
metadata explains absence: sampling disabled, ineligible, or not recorded.
`require_epli` returns the number and throws if the score is absent. Scoring itself
returns a validated score or throws; an error never becomes `nothing`.

CSV stores absent scores as `missing`, and JSON uses `null`. No separate state flag
is needed. Eligible candidates in sampled selection must have scores before ranking
or cache reuse. Disabled sampling permits absent scores. Older JSON results with an
unavailable score retain `:not_recorded`, which does not assert that the score was
never computed. Absence reasons are optional; supplied reasons must agree with the
score and sampling metadata. Absence is never converted into zero.

### Zero comparable positions

The explicit convention `(matched=0, comparable=0)` gives zero, and the result
records `no_comparable_positions=true`. This differs from comparing positions and
finding no matches. The zero component still enters the median and normal seed
ranking: ties use MSA size, then PID order. It is not silently removed from the
sample set. Per-sample CSVs retain both counts and this diagnostic; summaries report
`no_comparable_samples`. Inspect these diagnostics before interpreting a zero EPLI
as measured zero identity. Disabling sampling leaves the score absent instead.

HHsuite output must contain complete, recognizable query/comparison/target blocks.
Empty, unrecognized, or truncated output raises an error rather than becoming a
zero-comparison result. Residue coordinates start at one; zero is reserved for gaps
in the parsed position vector. A successfully parsed block without any reference positions
can still produce valid `(0,0)` counts. Optional HHsuite diagnostics remain separate.

Finite constructors do not guarantee that arithmetic is numerically accurate.
Percentage normalization uses higher precision near the Float64 limits to avoid
intermediate overflow and premature underflow. The median also preserves
representable subnormal results when averaging its two central values. Truly
unrepresentable tiny results round to signed zero; true overflow raises an error.
The normalization and median formulas are unchanged, and invalid computed results
raise errors rather than dropping samples.

Iduna uses EPLI internally to choose the best seed generated by ThorAxe. That
adapter uses the same scoring layer, but its sampling step is different: ThorAxe
does not normally accept arbitrary sampled sequence rows as input, so Iduna
samples species lists instead of sequences for ThorAxe PID selection.

```@docs
Iduna.EPLI.epli_score
Iduna.EPLI.hhsuite_identity_score
Iduna.EPLI.hhsuite_profile_score
Iduna.EPLI.comparable_positions_normalization
Iduna.EPLI.self_reference_normalization
Iduna.EPLI.no_normalization
Iduna.EPLI.prographmsa_aligner
```
