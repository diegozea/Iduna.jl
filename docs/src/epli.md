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

### Read the score

`result.epli` is an `EPLIScore` containing the median of the sample scores.
Use `epli_value` to get its numerical value. This function also accepts a ThorAxe
seed and raises `ArgumentError` if its score is unavailable.

```julia
value = epli_value(result)
```

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

By convention, a valid comparison with no comparable positions contributes zero
to the median. In `scores.csv`, `no_comparable_positions` is `true` for this case,
`false` when positions were compared (even if none matched), and `missing` when
the scoring function did not supply counts. Summaries count samples with
no comparable positions in `no_comparable_samples` and samples lacking counts
in `samples_without_identity_counts`.
Malformed HHsuite output raises an error instead of being treated as zero.

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

Use `no_normalization` if you want the median of the raw scores in their original
units. Negative profile scores and ratios above 100 are allowed; a reference self
score of zero cannot be used as a divisor. Scores must be finite numbers:
`NaN`, `Inf`, and `-Inf` raise errors.

## Write a Custom Scoring Function

The built-in functions above handle the usual workflows. To experiment with a
different comparison method, pass your function as `score_fn` to `epli_score`.
It receives two aligned FASTA paths, for the full reference and the sample,
and must accept the keywords `logs_dir` and `label`. Return one of these types:

| What your function calculates | What to return |
| --- | --- |
| Matched and comparable position counts | `IdentityScore(IdentityCounts(matched, comparable))` |
| A profile score or another numerical measure | `RawScore(value)` |

`IdentityCounts` checks that the counts are integers and satisfy
`0 ≤ matched ≤ comparable`. The raw value of an `IdentityScore` is its matched
count. `RawScore` stores the finite number you supply, including negative values.
Use `raw_score` to read either kind of score:

```jldoctest score_values
using Iduna.EPLI

count_score = IdentityScore(IdentityCounts(3, 4));
profile_score = RawScore(-2.0);
(raw_score(count_score), raw_score(profile_score))

# output

(3.0, -2.0)
```

[`ScoreResult`](@ref Iduna.EPLI.ScoreResult) is the name for
`Union{IdentityScore,RawScore}`, not an abstract supertype. Custom scoring methods
return these same types; you implement the comparison function rather than a new
score subtype. Pair count-based scores with `comparable_positions_normalization`,
or use `self_reference_normalization` or `no_normalization` for other measures.

### Optional information from scoring tools

The built-in HHsuite functions also provide the alignment probability, E-value,
and HHalign profile score in `score.diagnostics`. These fields are copied to
`scores.csv` to help you inspect each comparison. You do not need to fill them in
when using the built-in functions.

A custom scoring function can attach extra information for its normalizer to
read. For example:

```jldoctest score_values
score = RawScore(-2.0; diagnostics=(; e_value=1e-6));
score.diagnostics.e_value

# output

1.0e-6
```

Adding this information does not change the raw score. Optional diagnostic values
may be missing or non-finite; the numerical score itself must remain finite.

### Write a normalization function

Pass a custom function as `normalization_fn` to change how each sample score is
scaled. It receives the score returned by `score_fn` and the keyword
`reference_score`, and returns `NormalizationResult(value, denominator)`.
`value` enters the final median; the optional `denominator` is recorded in
`scores.csv`. Both supplied numbers must be finite.

For example, report identity as a fraction from 0 to 1 when comparing it with
other measures on that scale. This reuses the built-in calculation and its
zero-comparison convention, then changes the units from percent to fraction:

```jldoctest identity_fraction
using Iduna.EPLI

function identity_fraction(score::IdentityScore; reference_score=nothing)
    percent = comparable_positions_normalization(score)
    return NormalizationResult(percent.normalized_score / 100,
        percent.normalization_score)
end

score = IdentityScore(IdentityCounts(3, 4));
identity_fraction(score).normalized_score

# output

0.75
```

Use this function with the default identity scorer by passing
`normalization_fn=identity_fraction` to `epli_score`.

`reference_score` is normally `nothing`. If your normalizer needs the score from
comparing the reference MSA with itself, also pass
`compute_reference_self_score=true`. Iduna then calls your `score_fn` with the
reference path twice and passes that result to the normalizer. The built-in
`self_reference_normalization` requests this automatically.

## EPLI for ThorAxe Seed Selection

`iduna` uses EPLI to choose among the alignments generated by ThorAxe at different
percent identity (PID) thresholds. It samples species lists because ThorAxe does
not normally accept arbitrary sampled sequence rows as input. Each selected
alignment is stored as a seed in `result.thoraxe_msa.seeds`.

Setting `pid_sample_count=0` keeps all eligible candidates without computing EPLI.
Their `seed.epli` is `nothing`. When reading a seed from an `iduna` result, check
for an unavailable score before using its value:

```julia
seed = first(result.thoraxe_msa.seeds)
if seed.epli === nothing
    println("EPLI unavailable: ", seed.epli_reason)
else
    println(epli_value(seed))
end
```

An optional `epli_reason` explains why a score is unavailable: sampling was
disabled, the candidate was ineligible, or no score was recorded. When analysing
saved results, an empty EPLI cell in a CSV (read as `missing` in Julia) or `null`
in JSON also means unavailable. Do not interpret it as zero. Scoring errors stop
the calculation; they do not produce an unavailable score to rank or reuse.

## Function Reference

```@docs
Iduna.EPLI.epli_score
Iduna.EPLI.hhsuite_identity_score
Iduna.EPLI.hhsuite_profile_score
Iduna.EPLI.comparable_positions_normalization
Iduna.EPLI.self_reference_normalization
Iduna.EPLI.no_normalization
Iduna.EPLI.prographmsa_aligner
```
