using StaticArrays

## Coefficient and weight calculation (from bottom to top):

# Bottom-level coefficients (Response and Pathogen)

function responseStaticCoefficients(
    imprinted_pathogen_genome::String, matured_pathogen_genome::String, type::ResponseType)
    return [
        type.static_coefficient_functions[evt_id](
            imprinted_pathogen_genome, matured_pathogen_genome
        )
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

function responseSpecificCoefficient(pathogen::Pathogen, response::Response, coefficient::Int64)
    return response.type.specific_coefficient_functions[coefficient](
        response.imprinted_pathogen.sequence,
        response.matured_pathogen.sequence,
        pathogen.sequence
    )
end

function pathogenSequenceCoefficients(sequence::String, type::PathogenType)
    return [
        type.coefficient_functions[evt_id](sequence)
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

# Intrahost-level representation

function pathogenFractions!(host::Host)
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
                    host.pathogens[p].coefficients[INTRAHOST_FITNESS] * sum([
                        re.type.specific_coefficient_functions[INTRAHOST_FITNESS](
                            re.imprinted_pathogen.sequence,
                            re.matured_pathogen.sequence,
                            host.pathogens[p].sequence
                        )
                        for re in host.responses
                    ])
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
    host.pathogen_fractions = fracs
end

# Intrahost-level weights

function pathogenWeights!(p::Int64, host::Host, evt::Int64)
    host.pathogen_weights[evt, p] =
        host.pathogens[p].type.coefficient_functions[evt](
            host.pathogens[p].sequence
        )
    if evt == CLEARANCE
        # Clearance likelihoods are not proportional to intrahost fraction,
        # instead, we assume they are uniform (could be inversely proportional?)
        host.pathogen_weights[evt, p] =
            host.pathogen_weights[evt, p] /
            length(host.pathogens)
    else
        # All other pathogen event likelihoods are proportional to intrahost fraction
        host.pathogen_weights[evt, p] =
            host.pathogen_weights[evt, p] *
            host.pathogen_fractions[p]
    end
    if length(host.responses) > 0
        host.pathogen_weights[evt, p] =
            host.pathogen_weights[evt, p] *
            sum([
                re.type.specific_coefficient_functions[evt](
                    re.imprinted_pathogen.sequence,
                    re.matured_pathogen.sequence,
                    host.pathogens[p].sequence
                )
                for re in host.responses
            ])
    end
end

function hostWeightsPathogen!(host_idx::Int64, class::Class, evt::Int64)
    class.host_weights[evt, host_idx] = 0.0
    for p in 1:length(class.hosts[host_idx].pathogens)
        pathogenWeights!(p, class.hosts[host_idx], evt)
        class.host_weights[evt, host_idx] += class.hosts[host_idx].pathogen_weights[evt, p]
    end
end

# function pathogenWeights!(host_idx::Int64, class::Class)
#     for evt in PATHOGEN_EVENTS
#         pathogenWeights!(host_idx, class, evt)
#     end
# end

function responseWeights!(re::Int64, host::Host, evt::Int64)
    host.response_weights[evt-RESPONSE_EVENTS[1]+1, re] =
    # No Response fraction here, just presence of Response is enough
        host.responses[re].type.static_coefficient_functions[evt](
            host.responses[re].imprinted_pathogen.sequence,
            host.responses[re].matured_pathogen.sequence
        )
end

function hostWeightsResponse!(host_idx::Int64, class::Class, evt::Int64)
    class.host_weights[evt, host_idx] = 0.0
    for re in 1:length(class.hosts[host_idx].responses)
        responseWeights!(re, class.hosts[host_idx], evt)
        class.host_weights[evt, host_idx] += class.hosts[host_idx].response_weights[evt-RESPONSE_EVENTS[1]+1, re]
    end
end

# function responseWeights!(host_idx::Int64, class::Class)
#     for evt in RESPONSE_EVENTS
#         responseWeights!(host_idx, class, evt)
#     end
# end

# Intra-Class level

function hostWeightsHost!(h::Int64, class::Class, evt::Int64)
    class.host_weights[evt, h] = 1.0
    if length(class.hosts[h].pathogens) > 0
        class.host_weights[evt, h] =
            class.host_weights[evt, h] *
            sum([
                class.hosts[h].pathogen_fractions[p] *
                class.hosts[h].pathogens[p].type.coefficient_functions[evt](
                    class.hosts[h].pathogens[p].sequence
                )
                for p in 1:length(class.hosts[h].pathogens)
            ])
        if length(class.hosts[h].responses) > 0
            class.host_weights[evt, h] =
                class.host_weights[evt, h] *
                sum([
                    class.hosts[h].pathogen_fractions[p] *
                    re.type.specific_coefficient_functions[evt](
                        re.imprinted_pathogen.sequence,
                        re.matured_pathogen.sequence,
                        class.hosts[h].pathogens[p].sequence
                    )
                    for re in class.hosts[h].responses for p in 1:length(class.hosts[h].pathogens)
                ])
        end
    end
    if length(class.hosts[h].responses) > 0
        class.host_weights[evt, h] =
            class.host_weights[evt, h] *
            sum([
                re.type.static_coefficient_functions[evt](
                    re.imprinted_pathogen.sequence,
                    re.matured_pathogen.sequence
                )
                for re in class.hosts[h].responses
            ])
    end
end

function hostWeightsReceive!(h::Int64, class::Class, evt::Int64)
    class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] = 1.0
    if length(class.hosts[h].pathogens) > 0
        class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] =
            class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] *
            sum([
                class.hosts[h].pathogen_fractions[p] *
                class.hosts[h].pathogens[p].type.coefficient_functions[evt](
                    class.hosts[h].pathogens[p].sequence
                )
                for p in 1:length(class.hosts[h].pathogens)
            ])
        if length(class.hosts[h].responses) > 0
            class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] =
                class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] *
                sum([
                    class.hosts[h].pathogen_fractions[p] *
                    re.type.specific_coefficient_functions[evt](
                        class.hosts[h].pathogens[p].sequence,
                        re.imprinted_pathogen.sequence,
                        re.matured_pathogen.sequence
                    )
                    for re in class.hosts[h].responses for p in 1:length(class.hosts[h].pathogens)
                ])
        end
    end
    if length(class.hosts[h].responses) > 0
        class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] =
            class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] *
            sum([
                re.type.static_coefficient_functions[evt](
                    re.imprinted_pathogen.sequence,
                    re.matured_pathogen.sequence
                )
                for re in class.hosts[h].responses
            ])
    end
end

function hostWeights!(host_idx::Int64, class::Class, population::Population, model::Model)
    pathogenFractions!(class.hosts[host_idx])
    prev = 0.0
    for weight in PATHOGEN_EVENTS
        prev = class.host_weights[weight, host_idx]
        hostWeightsPathogen!(host_idx, class, weight)
        if prev != class.host_weights[weight, host_idx] && class.parameters.base_coefficients[weight] != 0.0
            propagateWeightChanges!(
                class.parameters.base_coefficients[weight] * (class.host_weights[weight, host_idx] - prev),
                class, population, weight, model
            )
        end
    end
    for weight in RESPONSE_EVENTS
        prev = class.host_weights[weight, host_idx]
        hostWeightsResponse!(host_idx, class, weight)
        if prev != class.host_weights[weight, host_idx] && class.parameters.base_coefficients[weight] != 0.0
            propagateWeightChanges!(
                class.parameters.base_coefficients[weight] * (class.host_weights[weight, host_idx] - prev),
                class, population, weight, model
            )
        end
    end
    for weight in HOST_EVENTS
        prev = class.host_weights[weight, host_idx]
        hostWeightsHost!(host_idx, class, weight)
        if prev != class.host_weights[weight, host_idx] && class.parameters.base_coefficients[weight] != 0.0
            propagateWeightChanges!(
                class.parameters.base_coefficients[weight] * (class.host_weights[weight, host_idx] - prev),
                class, population, weight, model
            )
        end
    end
    for weight in CHOICE_MODIFIERS[begin:end-1]
        prev = class.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, host_idx]
        hostWeightsReceive!(host_idx, class, weight)
        if prev != class.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, host_idx] && class.parameters.base_coefficients[weight] != 0.0 && weight < NUM_COEFFICIENTS - 2
            propagateWeightReceiveChanges!(
                class.parameters.base_coefficients[weight] *
                (class.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, host_idx] - prev),
                class, population, weight, model
            )
        end
    end
end

# Intra-Population level: None needed

# Intra-Model (Event) level: None needed

# Rate level: None needed

# Weight change propagation:

function propagateWeightChanges!(change::SVector{NUM_COEFFICIENTS,Float64}, population::Population, model::Model)
    model.population_weights[:, model.population_dict[population.id]] .+= change[begin:NUM_EVENTS]
    model.population_weights_receive[:, model.population_dict[population.id]] .+= change[end-NUM_CHOICE_MODIFIERS+1:end-3]

    model.population_weights_receive_sums .+= change[end-NUM_CHOICE_MODIFIERS+1:end-3]
    model.event_rates .+= change[begin:NUM_EVENTS]
    model.event_rates_sum += sum(change[begin:NUM_EVENTS])
end

function propagateWeightChanges!(change::SVector{NUM_COEFFICIENTS,Float64}, class::Class, population::Population, model::Model)
    population.class_weights[:, population.class_dict[class.id]] .+= change[begin:NUM_EVENTS]
    population.class_weights_receive[:, population.class_dict[class.id]] .+= change[end-NUM_CHOICE_MODIFIERS+1:end-1]

    propagateWeightChanges!(change, population, model)
end

function propagateWeightChanges!(
    change::SVector{NUM_COEFFICIENTS,Float64}, host_idx::Int64, class::Class, population::Population, model::Model)
    class.host_weights[:, host_idx] .+= change[begin:NUM_EVENTS]
    class.host_weights_receive[:, host_idx] .+= change[end-NUM_CHOICE_MODIFIERS+1:end-1]

    propagateWeightChanges!(class.parameters.base_coefficients .* change, class, population, model)
end

function propagateWeightChanges!(change::Float64, population::Population, evt::Int64, model::Model)
    model.population_weights[evt, model.population_dict[population.id]] += change

    model.event_rates[evt] += change
    model.event_rates_sum += change
end

function propagateWeightChanges!(change::Float64, class::Class, population::Population, evt::Int64, model::Model)
    population.class_weights[evt, population.class_dict[class.id]] += change

    propagateWeightChanges!(change, population, evt, model)
end

function propagateWeightChanges!(change::Float64, host_idx::Int64, class::Class, population::Population, evt::Int64, model::Model)
    class.host_weights[evt, host_idx] += change

    propagateWeightChanges!(class.parameters.base_coefficients[evt] * change, class, population, evt, model)
end

function propagateWeightReceiveChanges!(change::Float64, population::Population, evt::Int64, model::Model)
    model.population_weights_receive[evt-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]] += change

    model.population_weights_receive_sums[evt-CHOICE_MODIFIERS[1]+1] += change
end

function propagateWeightReceiveChanges!(change::Float64, class::Class, population::Population, evt::Int64, model::Model)
    population.class_weights_receive[evt-CHOICE_MODIFIERS[1]+1, population.class_dict[class.id]] += change

    if evt < NUM_COEFFICIENTS - 3
        propagateWeightReceiveChanges!(change, population, evt, model)
    end
end

function propagateWeightReceiveChanges!(change::Float64, host_idx::Int64, class::Class, population::Population, evt::Int64, model::Model)
    class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, host_idx] += change

    if evt < NUM_COEFFICIENTS - 1
        propagateWeightReceiveChanges!(class.parameters.base_coefficients[evt] * change, class, population, evt, model)
    end
end
