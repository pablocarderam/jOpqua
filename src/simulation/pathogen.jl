using StaticArrays

struct Pathogen
    sequence::String
    fitness::Float64
    coefficients::SVector{NUM_COEFFICIENTS,Float64}
end

function sequenceCoefficients(class::Class, sequence::String)
    return [class.parameters.pathogen_coefficient_functions[evt_id](sequence) for evt_id in NUM_COEFFICIENTS]
end

function sequenceFitness(class::Class, sequence::String)
    return class.parameters.pathogenFitness(sequence)
end

function newPathogen!(class::Class, sequence::String)
    pathogen = Pathogen(
        sequence,
        sequenceFitness(class, sequence),
        sequenceCoefficients(class, sequence),
    )
    push!(class.pathogens, pathogen)
end
