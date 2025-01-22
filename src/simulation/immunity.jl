function infectionProbability(pathogen::Pathogen, host::Host)
    # reactivity-weighted arithmetic mean of infection coefficients
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

    # # different algorithm reactivity-weighted arithmetic mean of infection coefficients
    # return sum([
    #     response.type.reactivityCoefficient(
    #         response.imprinted_pathogen.sequence,
    #         response.matured_pathogen.sequence,
    #         pathogen.sequence
    #     ) * response.type.infectionCoefficient(
    #         response.imprinted_pathogen.sequence,
    #         response.matured_pathogen.sequence,
    #         pathogen.sequence
    #     )
    #     for response in host.responses
    #     ]) / sum([
    #         response.type.reactivityCoefficient(
    #             response.imprinted_pathogen.sequence,
    #             response.matured_pathogen.sequence,
    #             pathogen.sequence
    #         )
    #         for response in host.responses
    #     ])

    # # reactivity-weighted harmonic mean of infection coefficients
    # reac_sum = 0.0
    # denom_sum = 0.0
    # for response in host.responses
    #     reac_sum += response.type.reactivityCoefficient(
    #         response.imprinted_pathogen.sequence,
    #         response.matured_pathogen.sequence,
    #         pathogen.sequence
    #     )
    #     denom_sum += response.type.reactivityCoefficient(
    #         response.imprinted_pathogen.sequence,
    #         response.matured_pathogen.sequence,
    #         pathogen.sequence
    #     ) / response.type.infectionCoefficient(
    #         response.imprinted_pathogen.sequence,
    #         response.matured_pathogen.sequence,
    #         pathogen.sequence
    #     )
    # end

    # return reac_sum / denom_sum

    # different algorithm reactivity-weighted harmonic mean of infection coefficients
    # return sum([ # reactivity-weighted harmonic mean of infection coefficients
    #     response.type.reactivityCoefficient(
    #         response.imprinted_pathogen.sequence,
    #         response.matured_pathogen.sequence,
    #         pathogen.sequence
    #     )
    #     for response in host.responses
    #     ]) / sum([
    #     response.type.reactivityCoefficient(
    #         response.imprinted_pathogen.sequence,
    #         response.matured_pathogen.sequence,
    #         pathogen.sequence
    #     ) / response.type.infectionCoefficient(
    #         response.imprinted_pathogen.sequence,
    #         response.matured_pathogen.sequence,
    #         pathogen.sequence
    #     )
    #     for response in host.responses
    #     ])
end
