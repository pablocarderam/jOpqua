function hostWeightsPathogen!(class::Class)
    for evt in PATHOGEN_EVENTS
        for h in 1:length(class.hosts)
            class.host_weights[evt, h] = sum(
                @views class.hosts[h].pathogen_weights[evt+PATHOGEN_EVENTS[1]-1, :]
            )
            # pathogen weights have already accounted for pathogen fraction
            # if length(class.hosts[h].immunities) > 0
            #     class.host_weights[evt, h] =
            #         class.host_weights[evt, h] *
            #         sum([
            #             class.hosts[h].immunity_fractions[im] *
            #             class.parameters.immunity_coefficient_functions[evt](
            #                 class.immunities[class.hosts[h].immunities[im]].sequence
            #             )
            #             for im in 1:length(class.hosts[h].immunities)
            #         ])
            # end
        end
    end
end

function hostWeightsImmunity!(class::Class)
    for evt in IMMUNITY_EVENTS
        for h in 1:length(class.hosts)
            class.host_weights[evt, h] = sum(
                @views class.hosts[h].immunity_weights[evt+IMMUNITY_EVENTS[1]-1, :]
            )
            # immunity weights have already accounted for immunity fraction
            # if length(class.hosts[h].pathogens) > 0
            #     class.host_weights[evt, h] =
            #         class.host_weights[evt, h] *
            #         sum([
            #             class.hosts[h].pathogen_fractions[p] *
            #             class.parameters.pathogen_coefficient_functions[evt](
            #                 class.pathogens[class.hosts[h].pathogens[p]].sequence
            #             )
            #             for p in 1:length(class.hosts[h].pathogens)
            #         ])
            # end
        end
    end
end

function hostWeightsHost!(class::Class)
    for evt in HOST_EVENTS
        for h in 1:length(class.hosts)
            class.host_weights[evt, h] = 1.0
            if length(class.hosts[h].pathogens) > 0
                class.host_weights[evt, h] =
                    class.host_weights[evt, h] *
                    sum([
                        class.hosts[h].pathogen_fractions[p] *
                        class.parameters.pathogen_coefficient_functions[evt](
                            class.pathogens[class.hosts[h].pathogens[p]].sequence
                        )
                        for p in 1:length(class.hosts[h].pathogens)
                    ])
                if length(class.hosts[h].immunities) > 0
                    class.host_weights[evt, h] =
                        class.host_weights[evt, h] *
                        sum([
                            class.hosts[h].pathogen_fractions[p] *
                            class.parameters.immunity_types[
                                class.immunities[im].type
                            ].specific_coefficient_functions[evt](
                                class.pathogens[class.hosts[h].pathogens[p]].sequence,
                                class.immunities[im].imprinted_sequence,
                                class.immunities[im].matured_sequence
                            )
                            for im in class.hosts[h].immunities for p in 1:length(class.hosts[h].pathogens)
                        ])
                end
            end
            if length(class.hosts[h].immunities) > 0
                class.host_weights[evt, h] =
                    class.host_weights[evt, h] *
                    sum([
                        # class.hosts[h].immunity_fractions[im] *
                        class.parameters.immunity_types[
                            class.immunities[im].type
                        ].static_coefficient_functions[evt](
                            class.immunities[im].imprinted_sequence,
                            class.immunities[im].matured_sequence
                        )
                        for im in class.hosts[h].immunities
                    ])
            end
        end
    end
end

function hostWeightsReceive!(class::Class)
    for evt in CHOICE_MODIFIERS[begin:end-1] # excludes intrahost fitness
        for h in 1:length(class.hosts)
            class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] = 1.0
            if length(class.hosts[h].pathogens) > 0
                class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] =
                    class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] *
                    sum([
                        class.hosts[h].pathogen_fractions[p] *
                        class.parameters.pathogen_coefficient_functions[evt](
                            class.pathogens[class.hosts[h].pathogens[p]].sequence
                        )
                        for p in 1:length(class.hosts[h].pathogens)
                    ])
                if length(class.hosts[h].immunities) > 0
                    class.host_weights[evt, h] =
                        class.host_weights[evt, h] *
                        sum([
                            class.hosts[h].pathogen_fractions[p] *
                            class.parameters.immunity_types[
                                class.immunities[im].type
                            ].specific_coefficient_functions[evt](
                                class.pathogens[class.hosts[h].pathogens[p]].sequence,
                                class.immunities[im].imprinted_sequence,
                                class.immunities[im].matured_sequence
                            )
                            for im in class.hosts[h].immunities for p in 1:length(class.hosts[h].pathogens)
                        ])
                end
            end
            if length(class.hosts[h].immunities) > 0
                class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] =
                    class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] *
                    sum([
                        class.parameters.immunity_types[
                            class.immunities[im].type
                        ].static_coefficient_functions[evt](
                            class.immunities[im].imprinted_sequence,
                            class.immunities[im].matured_sequence
                        )
                        for im in 1:length(class.hosts[h].immunities)
                    ])
            end
        end
    end
end

function hostWeights!(class::Class)
    hostWeightsPathogen!(class)
    hostWeightsImmunity!(class)
    hostWeightsHost!(class)
    hostWeightsReceive!(class)
end

function setRates!(class::Class)
    hostWeights!(class)
end

function newClass!(id::String, parameters::ClassParameters, population::Population)
    class = Class(
        id, parameters,
        Dict{Int64,Pathogen}(), Dict{Int64,Immunity}(), 0, 0,
        Vector{Host}(undef, 0),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
    )
    push!(population.classes, class)
end
