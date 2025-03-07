function weightedInteractionArithmeticMean(
        pathogen::Pathogen, host::Host, evt::Int64;
        coefficient=responseSpecificCoefficient)
    # reactivity-weighted arithmetic mean of specific coefficients
    if length(host.responses) > 0
        reac_sum = 0.0
        numerator_sum = 0.0
        reac = 0.0
        for response in host.responses
            reac = reactivityCoefficient(pathogen, response, host)
            reac_sum += reac
            numerator_sum += reac * coefficient(pathogen, response, host, evt)
        end

        return numerator_sum / reac_sum
    else
        return 1.0
    end
end

# function transmissionEfficiencyArithmeticMean(pathogen::Pathogen, host::Host)
#     # reactivity-weighted arithmetic mean of infection coefficients
#     if length(host.responses) > 0
#         reac_sum = 0.0
#         numerator_sum = 0.0
#         reac = 0.0
#         for response in host.responses
#             reac = reactivityCoefficient(pathogen, response, host)
#             reac_sum += reac
#             numerator_sum += reac * responseSpecificCoefficient(
#                 pathogen, response, host, TRANSMISSION_EFFICIENCY
#             )
#         end

#         return numerator_sum / reac_sum
#     else
#         return 1.0
#     end
# end

function weightedInteractionWinnerTakesAll(
        pathogen::Pathogen, host::Host, evt::Int64;
        coefficient=responseSpecificCoefficient)
    if length(host.responses) > 0
        dominant_reaction = 1.0
        dominant_reactivity = 0.0
        for response in host.responses
            reactivity = reactivityCoefficient(pathogen, response, host)
            if reactivity > dominant_reactivity
                dominant_reaction = reactivity * coefficient(pathogen, response, host, evt)
            end
        end

        return dominant_reaction
    else
        return 1.0
    end
end

# function transmissionEfficiencyWinnerTakesAll(pathogen::Pathogen, host::Host)
#     if length(host.responses) > 0
#         dominant_reaction = 1.0
#         dominant_reactivity = 0.0
#         for response in host.responses
#             reactivity = reactivityCoefficient(pathogen, response, host)
#             if reactivity > dominant_reactivity
#                 dominant_reaction = reactivity * responseSpecificCoefficient(
#                     pathogen, response, host, TRANSMISSION_EFFICIENCY
#                 )
#             end
#         end

#         return dominant_reaction
#     else
#         return 1.0
#     end
# end

function deNovoResponse(
    pathogen::Pathogen, host::Host,
    existing_responses::Dict{Tuple{String,String,String,String},Response},
    response_types::Dict{String,ResponseType}, birth_time::Float64; response_type_id::String="Default")
    if haskey(existing_responses, (host.sequence, pathogen.sequence, "", response_type_id))
        return [existing_responses[(host.sequence, pathogen.sequence, "", response_type_id)]]
    else
        return [newResponse!(
            pathogen, nothing, host.sequence, existing_responses, response_types[response_type_id],
            parents=MVector{2,Union{Response,Nothing}}([nothing, nothing]), birth_time=birth_time
        )]
    end
end
