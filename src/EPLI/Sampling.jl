function _sample_label(sample_idx::Integer)
    return "sequence_subset_$(lpad(string(sample_idx), 3, '0'))"
end

function _validate_sampling_options(sample_count::Integer, sample_fraction::Real)
    sample_count > 0 || error("sample_count must be positive.")
    0.0 < Float64(sample_fraction) <= 1.0 ||
        error("sample_fraction must be greater than 0 and at most 1.")
    return nothing
end

function _normalize_sample_seed(seed::Integer)::UInt64
    seed < 0 && error("sample_seed must be non-negative.")
    return UInt64(seed)
end

_ungapped_sequence(sequence) = replace(String(sequence), '-' => "", '.' => "")

function _sequence_records(path::AbstractString)
    sequences = read_file(path, FASTASequences)
    records = _sequence_records_from_sequences(sequences)
    isempty(records) && error("Input FASTA has no sequences: $(path).")
    names = [name for (name, _seq) in records]
    return (; names, records)
end

function _sequence_records(msa::AbstractMultipleSequenceAlignment)
    sequence_names = collect(sequencenames(msa))
    records = [(String(name), _ungapped_sequence(stringsequence(msa, name)))
               for name in sequence_names]
    isempty(records) && error("Input MSA has no sequences.")
    names = [name for (name, _seq) in records]
    return (; names, records)
end

function _sequence_records(seq::AbstractSequence)
    return _sequence_records([seq])
end

function _sequence_records(sequences::AbstractVector{<:AbstractSequence})
    records = _sequence_records_from_sequences(sequences)
    isempty(records) && error("Input sequence collection has no sequences.")
    names = [name for (name, _seq) in records]
    return (; names, records)
end

function _sequence_records_from_sequences(sequences)
    return [(String(sequence_id(seq)), _ungapped_sequence(stringsequence(seq)))
            for seq in sequences]
end

function _reference_index(names::AbstractVector{<:AbstractString}, reference_sequence)
    reference_sequence === nothing && return 1
    reference = String(reference_sequence)
    idx = findfirst(==(reference), String.(names))
    idx === nothing &&
        error("Reference sequence $(reference) was not found in the input sequences.")
    return idx
end
