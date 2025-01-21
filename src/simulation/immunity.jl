function immunityProbability(pathogen::Pathogen, host::Host)
    return max([
        response.type.immunityProbability(
            response.imprinted_pathogen.sequence,
            response.matured_pathogen.sequence,
            pathogen.sequence
        )
        for response in host.responses
        ])
end
