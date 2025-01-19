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

function pathogenWeights!(p::Int64, host::Host, class::Class, evt::Int64)
    host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] =
        host.pathogens[p].type.coefficient_functions[evt](
            host.pathogens[p].sequence
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
    if length(host.responses) > 0
        host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] =
            host.pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, p] *
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

function pathogenWeights!(p::Int64, host::Host, class::Class)
    for evt in PATHOGEN_EVENTS
        pathogenWeights!(p, host, class, evt)
    end
end

function pathogenWeights!(host::Host, class::Class)
    for p in 1:length(host.pathogens)
        pathogenWeights!(p, host, class)
    end
end

function responseWeights!(re::Int64, host::Host, evt::Int64)
    host.response_weights[evt-RESPONSE_EVENTS[1]+1, re] =
    # No Response fraction here, just presence of Response is enough
        host.responses[re].type.static_coefficient_functions[evt](
            host.responses[re].imprinted_pathogen.sequence,
            host.responses[re].matured_pathogen.sequence
        )
end

function responseWeights!(re::Int64, host::Host)
    for evt in RESPONSE_EVENTS
        responseWeights!(re, host, evt)
    end
end

function responseWeights!(host::Host, class::Class)
    for re in 1:length(host.responses)
        responseWeights!(re, host)
    end
end

# Intra-Class level

function hostWeightsPathogen!(h::Int64, class::Class, evt::Int64)
    class.host_weights[evt, h] = sum(
        @views class.hosts[h].pathogen_weights[evt-PATHOGEN_EVENTS[1]+1, :]
    )
end

function hostWeightsPathogen!(h::Int64, class::Class)
    for evt in PATHOGEN_EVENTS
        hostWeightsPathogen!(h, class, evt)
    end
end

function hostWeightsPathogen!(class::Class)
    for h in 1:length(class.hosts)
        hostWeightsPathogen!(h, class)
    end
end

function hostWeightsResponse!(h::Int64, class::Class, evt::Int64)
    class.host_weights[evt, h] = sum(
        @views class.hosts[h].response_weights[evt-RESPONSE_EVENTS[1]+1, :]
    )
end

function hostWeightsResponse!(h::Int64, class::Class)
    for evt in RESPONSE_EVENTS
        hostWeightsResponse!(h, class, evt)
    end
end

function hostWeightsResponse!(class::Class)
    for h in 1:length(class.hosts)
        hostWeightsResponse!(h, class)
    end
end

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

function hostWeightsHost!(h::Int64, class::Class)
    for evt in HOST_EVENTS
        hostWeightsHost!(h, class, evt)
    end
end

function hostWeightsHost!(class::Class)
    for h in 1:length(class.hosts)
        hostWeightsHost!(h, class)
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

function hostWeightsReceive!(h::Int64, class::Class)
    for evt in CHOICE_MODIFIERS[begin:end-1] # excludes intrahost fitness
        hostWeightsReceive!(h, class, evt)
    end
end

function hostWeightsReceive!(class::Class)
    for h in 1:length(class.hosts)
        hostWeightsReceive!(h, class)
    end
end

function hostWeights!(host_idx::Int64, class::Class, population::Population, model::Model)
    pathogenFractions!(class.hosts[host_idx])
    pathogenWeights!(class.hosts[host_idx], class)
    responseWeights!(class.hosts[host_idx], class)
    prev = 0.0
    for weight in PATHOGEN_EVENTS
        prev = class.host_weights[weight]
        hostWeightsPathogen!(host_idx, class, weight)
        if prev != class.host_weights[weight]
            propagateWeightChanges!(
                class.parameters.base_coefficients[weight] * (class.host_weights[weight] - prev),
                class, population, weight, model
            )
        end
    end
    for weight in RESPONSE_EVENTS
        prev = class.host_weights[weight]
        hostWeightsResponse!(host_idx, class, weight)
        if prev != class.host_weights[weight]
            propagateWeightChanges!(
                class.parameters.base_coefficients[weight] * (class.host_weights[weight] - prev),
                class, population, weight, model
            )
        end
    end
    for weight in HOST_EVENTS
        prev = class.host_weights[weight]
        hostWeightsHost!(host_idx, class, weight)
        if prev != class.host_weights[weight]
            propagateWeightChanges!(
                class.parameters.base_coefficients[weight] * (class.host_weights[weight] - prev),
                class, population, weight, model
            )
        end
    end
    for weight in CHOICE_MODIFIERS[begin:end-1]
        prev = class.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1]
        hostWeightsReceive!(host_idx, class, weight)
        if prev != class.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1] && weight < NUM_COEFFICIENTS - 2
            propagateWeightReceiveChanges!(
                class.parameters.base_coefficients[weight] *
                (class.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1] - prev),
                class, population, weight, model
            )
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
    responseWeights!(host, class)
    pathogenWeights!(host, class)
end

function setRates!(class::Class)
    hostWeightsPathogen!(class)
    hostWeightsResponse!(class)
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

function propagateWeightChanges!(
    change::SVector{NUM_COEFFICIENTS,Float64}, host_idx::Int64, class::Class, population::Population, model::Model)
    class.host_weights[:, host_idx] .+= change[begin:NUM_EVENTS]
    class.host_weights_receive[:, host_idx] .+= change[end-NUM_CHOICE_MODIFIERS+1:end-1]

    propagateWeightChanges!(class.parameters.base_coefficients .* change, class, population, model)
end

function propagateWeightChanges!(change::Float64, population::Population, evt::Int64, model::Model)
    model.population_weights[evt, model.population_dict[population.id]] += change

    model.event_rates[evt] += change
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
end

function propagateWeightReceiveChanges!(change::Float64, class::Class, population::Population, evt::Int64, model::Model)
    population.class_weights_receive[evt-CHOICE_MODIFIERS[1]+1, population.class_dict[class.id]] += change

    if evt < NUM_COEFFICIENTS - 3
        propagateWeightReceiveChanges!(change, population, evt, model)
    end
end

function propagateWeightReceiveChanges!(change::Float64, host_idx::Int64, class::Class, population::Population, evt::Int64, model::Model)
    class.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, host_idx] += change

    if evt < NUM_COEFFICIENTS - 2
        propagateWeightReceiveChanges!(class.parameters.base_coefficients[evt] * change, class, population, evt, model)
    end
end
