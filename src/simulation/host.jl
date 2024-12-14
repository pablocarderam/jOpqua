using StaticArrays
using Statistics

mutable struct Host
    id::Int64
    pathogens::Vector{Pathogen} # size MAX_PATHOGENS
    pathogen_fractions::Vector{Float64} # size MAX_PATHOGENS
    pathogen_weights::Matrix{Float64} # size NUM_PATHOGEN_EVENTS x MAX_PATHOGENS
    immunities::Vector{Immunity} # size MAX_IMMUNITIES
    immunity_fractions::Vector{Float64} # size MAX_IMMUNITIES
    immunity_weights::Matrix{Float64} # size NUM_IMMUNITY_EVENTS x MAX_IMMUNITIES
end

function immunityFractions!(host::Host, class::Class)
    # Currently, we assume immunity priority impacts all events equally regardless of event;
    # this is not necessarilly the case, however.
    # We also assume only the most dominant immunity affects events.
    # If multiple immunities are tied in dominance, the first one to have infected dominates.
    max_coef = argmax(
        [class.parameters.immunityDominance(host.immunities[i].sequence, host.pathogens) for i in 1:length(host.immunities)]
    )
    fracs = zeros(Float64, length(host.immunities))
    fracs[max_coef] = 1.0

    host.immunity_fractions = fracs

    # # We normalized based on immunodominance (maybe normalizing is wrong?):
    # weights = Vector{Float64}(undef, length(host.immunities))
    # total = 0.0
    # for i in 1:length(host.immunities)
    #     weights[i] = class.parameters.immunityDominance(host.immunities[i].sequence, host.pathogens)
    #     total += weights[i]
    # end
    # return weights ./ total
end

function pathogenFractions!(host::Host, class::Class)
    # Currently, we assume pathogen population fraction (share in total fitness)
    # impacts all events equally regardless of event; this is not necessarilly the case, however.
    # We also assume only the most fit pathogen affects events.
    # If multiple pathogens are tied in fitness, the first one to have infected dominates.
    # This all changes in popgen Opqua.
    fracs = Vector{Float64}(undef, length(host.pathogens))
    for p in 1:length(host.pathogens)
        fracs[p] =
            host.pathogens[p].fitness * sum([
                class.parameters.immunityEffectOnFitness(host.immunities[im].sequence, host.pathogens[p].sequence) *
                host.immunity_fractions[im]
                for im in 1:length(host.immunities)
            ])
    end
    max_coef = argmax(fracs)
    fracs = zeros(Float64, length(host.pathogens))
    fracs[max_coef] = 1.0
    host.pathogen_fractions = fracs
end

function pathogenWeights!(class::Class, host::Host)
    host.pathogen_weights = Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, length(host.pathogens))
    for evt in PATHOGEN_EVENTS
        for p in 1:length(host.pathogens)
            host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] =
                host.pathogen_fractions[p] *
                class.parameters.pathogen_coefficient_functions[evt](host.pathogens[p].sequence)
            if length(host.immunities) > 0
                host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] =
                    host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] *
                    sum([
                        host.immunity_fractions[im] *
                        class.parameters.immunity_coefficient_effect_functions[evt](host.immunities[im].sequence, host.pathogens[p].sequence)
                        for im in 1:length(host.immunities)
                    ])
            end
        end
    end
end

function immunityWeights!(class::Class, host::Host)
    host.immunity_weights = Matrix{Float64}(undef, NUM_IMMUNITY_EVENTS, length(host.immunities))
    for evt in IMMUNITY_EVENTS
        for im in 1:length(host.immunities)
            host.immunity_weights[evt-IMMUNITY_EVENTS[1]+1, im] =
            # No immunity fraction here, just presence of immunity is enough
                class.parameters.immunity_coefficient_functions[evt](host.immunities[im].sequence)
            if length(host.pathogens) > 0
                host.immunity_weights[evt-IMMUNITY_EVENTS[1]+1, im] =
                    host.immunity_weights[evt-IMMUNITY_EVENTS[1]+1, im] * sum([
                        # No pathogen fraction here, just presence of pathogen is enough
                        class.parameters.immunity_coefficient_effect_functions[evt](host.immunities[im].sequence, host.pathogens[p].sequence)
                        for p in 1:length(host.pathogens)
                    ])
            end
        end
    end
end

function newHost!(class::Class)
    host = Host(
        length(class.hosts) + 1,
        Vector{Pathogen}(undef, 0), Vector{Float64}(undef, 0), Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, 0),
        Vector{Immunity}(undef, 0), Vector{Float64}(undef, 0), Matrix{Float64}(undef, NUM_IMMUNITY_EVENTS, 0),
    )
    immunityFractions!(host, class)
    pathogenFractions!(host, class)
    immunityWeights!(class, host)
    pathogenWeights!(class, host)
    push!(class.hosts, host)
end
