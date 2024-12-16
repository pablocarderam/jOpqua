using StaticArrays

struct Pathogen
    sequence::String
    coefficients::SVector{NUM_COEFFICIENTS,Float64}
end

function sequenceCoefficients(sequence::String, class::Class)
    return [class.parameters.pathogen_coefficient_functions[evt_id](sequence) for evt_id in NUM_COEFFICIENTS]
end

function newPathogen!(sequence::String, class::Class)
    pathogen = Pathogen(sequence, sequenceCoefficients(sequence, class))
    push!(class.pathogens, pathogen)
end
