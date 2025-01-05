using StaticArrays

## Coefficient and weight calculation (from bottom to top):

# Bottom-level coefficients (Immunity and Pathogen)

function immunityStaticCoefficients(
    imprinted_pathogen::String, matured_pathogen::String, type::ImmunityType, class::Class)
    return [
        type.static_coefficient_functions[evt_id](
            imprinted_pathogen, matured_pathogen
        )
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

function immunitySpecificCoefficient(pathogen::Pathogen, immunity::Immunity, class::Class, coefficient::Int64)
    return class.parameters.immunity_types[immunity.type].specific_coefficient_functions[coefficient](
        class.pathogens[immunity.imprinted_pathogen].sequence,
        class.pathogens[immunity.matured_pathogen].sequence,
        pathogen.sequence
    )
end

function pathogenSequenceCoefficients(sequence::String, class::Class)
    return [
        class.parameters.pathogen_coefficient_functions[evt_id](sequence)
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

# Intrahost-level representation

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
                class.parameters.immunity_types[
                    class.immunities[im].type
                ].specific_coefficient_functions[INTRAHOST_FITNESS](
                    class.immunities[im].imprinted_sequence,
                    class.immunities[im].matured_sequence,
                    class.pathogens[host.pathogens[p]].sequence
                )
                for im in host.immunities
            ])
    end
    max_coef = argmax(fracs)
    fracs = zeros(Float64, length(host.pathogens))
    fracs[max_coef] = 1.0
    host.pathogen_fractions = fracs
end

# Intrahost-level weights

function pathogenWeights!(host::Host, class::Class)
    for evt in PATHOGEN_EVENTS
        for p in 1:length(host.pathogens)
            host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] =
                class.parameters.pathogen_coefficient_functions[evt](
                    class.pathogens[host.pathogens[p]].sequence
                )
            if evt == CLEARANCE
                # Clearance likelihoods are not proportional to intrahost fraction,
                # instead, we assume they are uniform (could be inversely proportional?)
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
                        class.parameters.immunity_types[
                            class.immunities[im].type
                        ].specific_coefficient_functions[evt](
                            class.immunities[im].imprinted_sequence,
                            class.immunities[im].matured_sequence,
                            class.pathogens[host.pathogens[p]].sequence
                        )
                        for im in host.immunities
                    ])
            end
        end
    end
end

function immunityWeights!(host::Host, class::Class)
    for evt in IMMUNITY_EVENTS
        for im in 1:length(host.immunities)
            host.immunity_weights[evt-IMMUNITY_EVENTS[1]+1, im] =
            # No immunity fraction here, just presence of immunity is enough
                class.parameters.immunity_types[
                    class.immunities[im].type
                ].static_coefficient_functions[evt](
                    class.immunities[im].imprinted_sequence,
                    class.immunities[im].matured_sequence
                )
        end
    end
end

# Intra-Class level

function hostWeightsPathogen!(class::Class)
    for evt in PATHOGEN_EVENTS
        for h in 1:length(class.hosts)
            class.host_weights[evt, h] = sum(
                @views class.hosts[h].pathogen_weights[evt+PATHOGEN_EVENTS[1]-1, :]
            )
        end
    end
end

function hostWeightsImmunity!(class::Class)
    for evt in IMMUNITY_EVENTS
        for h in 1:length(class.hosts)
            class.host_weights[evt, h] = sum(
                @views class.hosts[h].immunity_weights[evt+IMMUNITY_EVENTS[1]-1, :]
            )
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

# Intra-Population level

function classWeights!(population::Population)
    for evt in EVENTS
        for c in 1:length(population.classes)
            population.class_weights[evt, c] =
                sum(@views population.classes[c].host_weights[evt, :]) *
                population.classes[c].parameters.base_coefficients[evt]
        end
    end
end

function classWeightsReceive!(population::Population)
    for evt in CHOICE_MODIFIERS[begin:end-2]
        for c in 1:length(population.classes)
            population.class_weights_receive[evt-CHOICE_MODIFIERS[1]+1, c] =
                sum(@views population.classes[c].host_weights_receive[evt, :]) *
                population.classes[c].parameters.base_coefficients[evt]
        end
    end
end

# Intra-Model (Population) level

function populationWeights!(model::Model)
    for evt in EVENTS
        for p in 1:length(model.populations)
            model.population_weights[evt, p] = sum(@views model.populations[p].class_weights[evt, :])
        end
    end
end

function populationWeightsReceive!(model::Model)
    for evt in CHOICE_MODIFIERS[begin:end-2]
        for p in 1:length(model.populations)
            model.population_weights_receive[evt-CHOICE_MODIFIERS[1]+1, p] =
                sum(@views model.populations[p].class_weights_receive[evt, :])
        end
    end
end

# Intra-Model (Event) level

function rates!(model::Model)
    for evt in EVENTS
        model.event_rates[evt] = sum([
            model.population_weights[evt, p]
            #TODO: This loop has to be a switch with a case for every event type computing the rate the correct way
            for p in 1:length(model.populations)
        ])
    end
end

# Rate calculation:

function setRates!(host::Host, class::Class)
    pathogenFractions!(host, class)
    immunityWeights!(host, class)
    pathogenWeights!(host, class)
end

function setRates!(class::Class)
    hostWeightsPathogen!(class)
    hostWeightsImmunity!(class)
    hostWeightsHost!(class)
    hostWeightsReceive!(class)
end

function setRates!(population::Population)
    classWeights!(population)
    classWeightsReceive!(population)
end

function setRates!(model::Model)
    populationWeights!(model)
    populationWeightsReceive!(model)
    rates!(model)
end

# Weight change propagation:

function propagateWeightChanges!(change::SVector{NUM_COEFFICIENTS,Float64}, population::Population, model::Model)
    model.population_weights[:, model.population_dict[population.id]] .+= change[begin:NUM_EVENTS]
    model.population_weights_receive[:, model.population_dict[population.id]] .+= change[end-NUM_CHOICE_MODIFIERS+1:end-3]

    model.event_rates .+= change[begin:NUM_EVENTS]
end

function propagateWeightChanges!(change::SVector{NUM_COEFFICIENTS,Float64}, class::Class, population::Population, model::Model)
    population.class_weights[:, population.class_dict[class.id]] .+= change[begin:NUM_EVENTS]
    population.class_weights_receive[:, population.class_dict[class.id]] .+= change[end-NUM_CHOICE_MODIFIERS+1:end-2]

    propagateWeightChanges!(change, population, model)
end

function propagateWeightChanges!(change::SVector{NUM_COEFFICIENTS,Float64}, host_idx::Int64, class::Class, population::Population, model::Model)
    class.host_weights[:, host_idx] .+= class.parameters.base_coefficients[begin:NUM_EVENTS] .* change[begin:NUM_EVENTS] #class.host_weights[:,host_idx]
    class.host_weights_receive[:, host_idx] .+= class.parameters.base_coefficients[end-NUM_CHOICE_MODIFIERS+1:end-1] .* change[end-NUM_CHOICE_MODIFIERS+1:end-1] #class.host_weights_receive[begin:end-1,host_idx]

    propagateWeightChanges!(class.parameters.base_coefficients .* change, class, population, model)
end
