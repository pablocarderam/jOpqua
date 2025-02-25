function weightedResponseArithmeticMean(pathogen::Pathogen, host::Host, evt::Int64)
    # reactivity-weighted arithmetic mean of specific coefficients
    if length(host.responses) > 0
        reac_sum = 0.0
        numerator_sum = 0.0
        for response in host.responses
            reac_sum += response.type.reactivityCoefficient(
                host.sequence,
                isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
                isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
                pathogen.sequence
            )
            numerator_sum += response.type.reactivityCoefficient(
                host.sequence,
                isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
                isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
                pathogen.sequence
            ) * response.type.specific_coefficient_functions[evt](
                host.sequence,
                isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
                isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
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
                host.sequence,
                isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
                isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
                pathogen.sequence
            )
            numerator_sum += response.type.reactivityCoefficient(
                host.sequence,
                isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
                isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
                pathogen.sequence
            ) * response.type.infectionCoefficient(
                host.sequence,
                isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
                isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
                pathogen.sequence
            )
        end

        return numerator_sum / reac_sum
    else
        return 1.0
    end
end

function weightedResponseWinnerTakesAll(pathogen::Pathogen, host::Host, evt::Int64)
    if length(host.responses) > 0
        dominant_reaction = 1.0
        for response in host.responses
            reaction = response.type.reactivityCoefficient(
                host.sequence,
                isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
                isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
                pathogen.sequence
            ) * response.type.specific_coefficient_functions[evt](
                host.sequence,
                isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
                isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
                pathogen.sequence
            )
            if reaction < dominant_reaction
                dominant_reaction = reaction
            end
        end

        return dominant_reaction
    else
        return 1.0
    end
end

function infectionProbabilityWinnerTakesAll(pathogen::Pathogen, host::Host)
    if length(host.responses) > 0
        dominant_reaction = 1.0
        for response in host.responses
            reaction = response.type.reactivityCoefficient(
                host.sequence,
                isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
                isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
                pathogen.sequence
            ) * response.type.infectionCoefficient(
                host.sequence,
                isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
                isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
                pathogen.sequence
            )
            if reaction < dominant_reaction
                dominant_reaction = reaction
            end
        end

        return dominant_reaction
    else
        return 1.0
    end
end

function deNovoResponse(
        pathogen::Pathogen, host::Host,
        existing_responses::Dict{Tuple{String,String,String,String},Response},
        response_types::Dict{String,ResponseType}; response_type_id::String="Default")
    if haskey(existing_responses, (host.sequence, pathogen.sequence, "", response_type_id))
        return [existing_responses[(host.sequence, pathogen.sequence, "", response_type_id)]]
    else
        return [newResponse!(
            pathogen, nothing, host.sequence, existing_responses, response_types[response_type_id],
            parents=MVector{2,Union{Response,Nothing}}([nothing, nothing])
        )]
    end
end
