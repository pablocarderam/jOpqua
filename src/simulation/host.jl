mutable struct Host
    id::Int64

    pathogens::Vector{Pathogen} # size MAX_PATHOGENS
    immunities::Vector{Immunity} # size MAX_IMMUNITIES

    pathogen_fractions::Vector{Float64} # size MAX_PATHOGENS
    immunity_fractions::Vector{Float64} # size MAX_IMMUNITIES

    pathogen_weights::Matrix{Float64} # size NUM_PATHOGEN_EVENTS x MAX_PATHOGENS
    immunity_weights::Matrix{Float64} # size NUM_IMMUNITY_EVENTS x MAX_IMMUNITIES

    # We could handle receiving rates at this level, but the only one relevant
    # to entities within hosts is recombination, and that uses the same rates
    # already calculated above since it's symmetrical. We therefore don't define
    # any additional matrices here and we handle all receiving rates for hosts
    # and larger at the Class level.
end

function immunityFractions!(host::Host, class::Class)
    # Currently, we assume immunity priority impacts all events equally regardless of event;
    # this is not necessarilly the case, however.
    # We also assume only the most dominant immunity affects events.
    # If multiple immunities are tied in dominance, the first one to have infected dominates.
    # The values in the returned vector must sum to 1.0.
    max_coef = argmax(
        [class.parameters.immunityDominance(host.pathogens, host.immunities[i].sequence) for i in 1:length(host.immunities)]
    )
    fracs = zeros(Float64, length(host.immunities))
    fracs[max_coef] = 1.0

    host.immunity_fractions = fracs
end

function pathogenFractions!(host::Host, class::Class)
    # Currently, we assume pathogen population fraction (share in total fitness)
    # impacts all events equally regardless of event; this is not necessarilly the case, however.
    # We also assume only the most fit pathogen affects events.
    # If multiple pathogens are tied in fitness, the first one to have infected dominates.
    # This all changes in popgen Opqua.
    # # The values in the returned vector must sum to 1.0.
    fracs = Vector{Float64}(undef, length(host.pathogens))
    for p in 1:length(host.pathogens)
        fracs[p] =
            host.pathogens[p].coefficients[INTRAHOST_FITNESS] * sum([
                class.parameters.immunity_coefficient_effect_functions[INTRAHOST_FITNESS](host.pathogens[p].sequence, host.immunities[im].sequence) *
                host.immunity_fractions[im]
                for im in 1:length(host.immunities)
            ])
    end
    max_coef = argmax(fracs)
    fracs = zeros(Float64, length(host.pathogens))
    fracs[max_coef] = 1.0
    host.pathogen_fractions = fracs
end

function pathogenWeights!(host::Host, class::Class)
    # host.pathogen_weights = Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, length(host.pathogens))
    for evt in PATHOGEN_EVENTS
        for p in 1:length(host.pathogens)
            host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] =
                class.parameters.pathogen_coefficient_functions[evt](host.pathogens[p].sequence)
            if evt == CLEARANCE
                # Clearance likelihoods are not proportional to intrahost fraction
                host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] =
                    host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] /
                    length(host.pathogens)
            else
                # All other pathogen event likelihoods are proportional to intrahost fraction
                host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] =
                    host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] *
                    host.pathogen_fractions[p]
            end
            if length(host.immunities) > 0
                host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] =
                    host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] *
                    sum([
                        host.immunity_fractions[im] *
                        class.parameters.immunity_coefficient_effect_functions[evt](host.pathogens[p].sequence, host.immunities[im].sequence)
                        for im in 1:length(host.immunities)
                    ])
            end
        end
    end
end

function immunityWeights!(host::Host, class::Class)
    # host.immunity_weights = Matrix{Float64}(undef, NUM_IMMUNITY_EVENTS, length(host.immunities))
    for evt in IMMUNITY_EVENTS
        for im in 1:length(host.immunities)
            host.immunity_weights[evt-IMMUNITY_EVENTS[1]+1, im] =
            # No immunity fraction here, just presence of immunity is enough
                class.parameters.immunity_coefficient_functions[evt](host.immunities[im].sequence)
            if length(host.pathogens) > 0
                host.immunity_weights[evt-IMMUNITY_EVENTS[1]+1, im] =
                    host.immunity_weights[evt-IMMUNITY_EVENTS[1]+1, im] * sum([
                        # No pathogen fraction here, just presence of pathogen is enough
                        class.parameters.immunity_coefficient_effect_functions[evt](host.pathogens[p].sequence, host.immunities[im].sequence)
                        for p in 1:length(host.pathogens)
                    ])
            end
        end
    end
end

function setRates!(host::Host, class::Class)
    immunityFractions!(host, class)
    pathogenFractions!(host, class)
    immunityWeights!(host, class)
    pathogenWeights!(host, class)
end

function newHost!(class::Class)
    host = Host(
        length(class.hosts) + 1,
        Vector{Pathogen}(undef, 0), Vector{Immunity}(undef, 0),
        Vector{Float64}(undef, 0), Vector{Float64}(undef, 0),
        Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, 0),
        Matrix{Float64}(undef, NUM_IMMUNITY_EVENTS, 0),
    )
    push!(class.hosts, host)
end
