function weightedResponseArithmeticMean(pathogen::Pathogen, host::Host, evt::Int64)
    # reactivity-weighted arithmetic mean of specific coefficients
    if length(host.responses) > 0
        reac_sum = 0.0
        numerator_sum = 0.0
        for response in host.responses
            reac_sum += response.type.reactivityCoefficient(
                response.imprinted_pathogen.sequence,
                response.matured_pathogen.sequence,
                pathogen.sequence
            )
            numerator_sum += response.type.reactivityCoefficient(
                response.imprinted_pathogen.sequence,
                response.matured_pathogen.sequence,
                pathogen.sequence
            ) * response.type.specific_coefficient_functions[evt](
                response.imprinted_pathogen.sequence,
                response.matured_pathogen.sequence,
                pathogen.sequence
            )
        end

        return numerator_sum / reac_sum
    else
        return 1.0
    end
end

function infectionProbabilityArithmeticMean(pathogen::Pathogen, host::Host)
    # reactivity-weighted arithmetic mean of infection coefficients
    if length(host.responses) > 0
        reac_sum = 0.0
        numerator_sum = 0.0
        for response in host.responses
            reac_sum += response.type.reactivityCoefficient(
                response.imprinted_pathogen.sequence,
                response.matured_pathogen.sequence,
                pathogen.sequence
            )
            numerator_sum += response.type.reactivityCoefficient(
                response.imprinted_pathogen.sequence,
                response.matured_pathogen.sequence,
                pathogen.sequence
            ) * response.type.infectionCoefficient(
                response.imprinted_pathogen.sequence,
                response.matured_pathogen.sequence,
                pathogen.sequence
            )
        end

        return numerator_sum / reac_sum
    else
        return 1.0
    end
end
