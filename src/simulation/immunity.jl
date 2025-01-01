function immunityStaticCoefficients(
    imprinted_pathogen::String, matured_pathogen::String, type::String, class::Class)
    return [
        class.parameters.immunity_types[type].static_coefficient_functions[evt_id](
            imprinted_pathogen, matured_pathogen
        )
        for evt_id in NUM_COEFFICIENTS
    ]
end

function immunitySpecificCoefficient(pathogen::Pathogen, immunity::Immunity, class::Class, coefficient::Int64)
    return class.parameters.immunity_types[immunity.type].specific_coefficient_functions[coefficient](
        class.pathogens[immunity.imprinted_pathogen].sequence,
        class.pathogens[immunity.matured_pathogen].sequence,
        pathogen.sequence
    )
end

function newImmunity!(imprinted_pathogen::Pathogen, matured_pathogen::Pathogen, class::Class, type::String)
    class.immunities[class.immunity_count] = Immunity(
        class.immunity_count, imprinted_pathogen.id, matured_pathogen.id,
        immunityStaticCoefficients(imprinted_pathogen.sequence, matured_pathogen.sequence, type, class),
        type
    )
    class.immunity_count += 1
end
