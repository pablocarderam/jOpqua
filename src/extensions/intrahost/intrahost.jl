function pathogenFractionsWinnerTakesAll(host::Host, weightedResponse::FunctionWrapper{Float64,Tuple{Pathogen,Host,Int64}})
    # Currently, we assume pathogen population fraction (share in total fitness)
    # impacts all events equally regardless of event; this is not necessarilly the case, however.
    # We also assume only the most fit pathogen affects events.
    # If multiple pathogens are tied in fitness, the first one to have infected dominates.
    # This all changes in popgen Opqua.
    # # The values in the returned vector must sum to 1.0.
    fracs = Vector{Float64}(undef, length(host.pathogens))
    if length(host.pathogens) > 0
        if length(host.responses) > 0
            for p in 1:length(host.pathogens)
                fracs[p] =
                    host.pathogens[p].coefficients[INTRAHOST_FITNESS] * weightedResponse(
                        host.pathogens[p], host, INTRAHOST_FITNESS
                    )
            end
        else
            for p in 1:length(host.pathogens)
                fracs[p] = host.pathogens[p].coefficients[INTRAHOST_FITNESS]
            end
        end
        max_coef = argmax(fracs)
        fracs = zeros(Float64, length(host.pathogens))
        fracs[max_coef] = 1.0
    end

    return fracs
end

function pathogenFractionsProportionalFitness(host::Host, weightedResponse::FunctionWrapper{Float64,Tuple{Pathogen,Host,Int64}})
    fracs = Vector{Float64}(undef, length(host.pathogens))
    if length(host.pathogens) > 0
        if length(host.responses) > 0
            for p in 1:length(host.pathogens)
                fracs[p] =
                    host.pathogens[p].coefficients[INTRAHOST_FITNESS] * weightedResponse(
                        host.pathogens[p], host, INTRAHOST_FITNESS
                    )
            end
        else
            for p in 1:length(host.pathogens)
                fracs[p] = host.pathogens[p].coefficients[INTRAHOST_FITNESS]
            end
        end
        fracs = fracs ./ sum(fracs)
    end

    return fracs
end
