using Iduna, Test

struct InterruptingAlignmentModel end
Base.getproperty(::InterruptingAlignmentModel, ::Symbol) = throw(InterruptException())

@testset "Loaded protein alignment model identity" begin
    utils = Iduna.Utils
    fingerprint = utils._protein_alignment_score_identity
    model = utils._PROTEIN_ALIGNMENT_SCORE_MODEL
    original_matrix = model.submat
    original_open, original_extend = model.gap_open, model.gap_extend
    original = fingerprint()
    implementation = utils._implementation_identity()
    @test original !== nothing
    @test implementation.protein_alignment_sha256 == original
    @test fingerprint(nothing) === nothing
    @test_throws InterruptException fingerprint(InterruptingAlignmentModel())
    alphabet = utils.BioSequences.symbols(utils.BioSequences.AminoAcidAlphabet())
    @test Set(Char.(alphabet)) == Set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*-")
    # A hand-specified constant matrix fixes the serialization order and score coverage.
    constant_matrix = copy(original_matrix)
    for a in utils.LongAA("ABCDEFGHIJKLMNOPQRSTUVWXYZ*"),
        b in utils.LongAA("ABCDEFGHIJKLMNOPQRSTUVWXYZ*")

        constant_matrix[a, b] = 1
    end
    constant_model = utils.AffineGapScoreModel(constant_matrix; gap_open = -3, gap_extend = -2)
    @test fingerprint(constant_model) == utils._identity_hash((;
        model = string(typeof(constant_model)), matrix = string(typeof(constant_matrix)),
        alphabet = "ARNDCQEGHILKMFPSTWYVOUBJZX*", scores = fill(1, 27, 27),
        gap_open = -3, gap_extend = -2))
    try
        # A separately serialized Iduna model need not alias BioAlignments.BLOSUM62.
        model.submat = copy(original_matrix)
        @test model.submat !== utils.BLOSUM62
        @test fingerprint() == original
        @test utils._implementation_identity() == implementation
        for residue in utils.LongAA("ABCDEFGHIJKLMNOPQRSTUVWXYZ*")
            partner = first(utils.LongAA("A"))
            value = model.submat[residue, partner]
            model.submat[residue, partner] = value + 1
            @test fingerprint() != original
            @test utils.BLOSUM62[residue, partner] == value
            model.submat[residue, partner] = value
        end
        @test fingerprint() == original
        model.submat['A', 'X'] += 1
        forward = fingerprint()
        model.submat['A', 'X'] -= 1
        model.submat['X', 'A'] += 1
        @test fingerprint() != forward
        @test utils._implementation_identity() != implementation
        model.submat['X', 'A'] -= 1
        model.gap_open -= 1
        @test fingerprint() != original
        model.gap_open = original_open
        model.gap_extend -= 1
        @test fingerprint() != original
        model.gap_extend = original_extend
        model.submat = utils.BioAlignments.DichotomousSubstitutionMatrix(4, -1)
        @test fingerprint() === nothing
        @test_logs (:warn, r"Protein alignment model identity unavailable") @test utils._implementation_identity() ===
                                                                                  nothing
    finally
        model.submat = original_matrix
        model.gap_open, model.gap_extend = original_open, original_extend
    end
    @test fingerprint() == original
    @test utils._implementation_identity() == implementation
end
