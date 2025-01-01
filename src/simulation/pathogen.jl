function pathogenSequenceCoefficients(sequence::String, class::Class)
    return [
        class.parameters.pathogen_coefficient_functions[evt_id](sequence)
        for evt_id in NUM_COEFFICIENTS
    ]
end

function newPathogen!(sequence::String, class::Class)
    class.pathogens[class.pathogen_count] = Pathogen(
        class.pathogen_count, sequence, pathogenSequenceCoefficients(sequence, class)
    )
    class.pathogen_count += 1
end
