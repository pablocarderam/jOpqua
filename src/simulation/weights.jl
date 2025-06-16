using StaticArrays
using Flexle

## Coefficient and weight calculation (from bottom to top):

# Bottom-level coefficients (Response and Pathogen)

function responseStaticSpecificCoefficients(
    host_sequence::String, imprinted_pathogen_sequence::String,
    matured_pathogen_sequence::String, type::ResponseType)
    return [
        type.static_specific_coefficient_functions[evt_id](
            host_sequence, imprinted_pathogen_sequence, matured_pathogen_sequence
        )
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

function responseStaticHostwideCoefficients(
    host_sequence::String, imprinted_pathogen_sequence::String,
    matured_pathogen_sequence::String, type::ResponseType)
    return [
        type.static_hostwide_coefficient_functions[evt_id](
            host_sequence, imprinted_pathogen_sequence, matured_pathogen_sequence
        )
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

function responseStaticSpecificCoefficient(response::Response, host::Host, coefficient::Int64)
    return response.type.static_specific_coefficient_functions[coefficient](
        host.sequence,
        isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
        isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
    )
end

function responseInteractionSpecificCoefficient(pathogen::Pathogen, response::Response, host::Host, coefficient::Int64)
    return response.type.interaction_specific_coefficient_functions[coefficient](
        host.sequence,
        isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
        isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
        pathogen.sequence
    )
end

function responseStaticHostwideCoefficient(response::Response, host::Host, coefficient::Int64)
    return response.type.static_hostwide_coefficient_functions[coefficient](
        host.sequence,
        isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
        isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
    )
end

function responseInteractionHostwideCoefficient(pathogen::Pathogen, response::Response, host::Host, coefficient::Int64)
    return response.type.interaction_hostwide_coefficient_functions[coefficient](
        host.sequence,
        isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
        isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
        pathogen.sequence
    )
end

function reactivityCoefficient(pathogen::Pathogen, response::Response, host::Host)
    return response.type.reactivityCoefficient(
        host.sequence,
        isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
        isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
        pathogen.sequence
    )
end

function nonsamplingValue(coef::Int64, host::Host, population::Population)
    return host.coefficients[coef] * hostwideNetCoefficient(host, population, coef) * population.parameters.base_coefficients[coef]
end

function nonsamplingValue(coef::Int64, pathogen::Pathogen, host::Host, population::Population)
    return pathogen.specific_coefficients[coef] * population.parameters.weightedInteractionPathogen(
               pathogen, host, coef
           ) * host.coefficients[coef] * hostwideNetCoefficient(host, population, coef) * population.parameters.base_coefficients[coef]
end

function nonsamplingValue(coef::Int64, response::Response, host::Host, population::Population)
    return response.specific_coefficients[coef] * population.parameters.weightedInteractionResponse(
               response, host, coef
           ) * host.coefficients[coef] * hostwideNetCoefficient(host, population, coef) * population.parameters.base_coefficients[coef]
end

function pathogenSequenceSpecificCoefficients(sequence::String, type::PathogenType)
    return [
        type.specific_coefficient_functions[evt_id](sequence)
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

function pathogenSequenceHostwideCoefficients(sequence::String, type::PathogenType)
    return [
        type.hostwide_coefficient_functions[evt_id](sequence)
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

function hostSequenceCoefficients(sequence::String, type::HostType)
    return MVector{NUM_COEFFICIENTS,Float64}([
        type.coefficient_functions[evt_id](sequence)
        for evt_id in 1:NUM_COEFFICIENTS
    ])
end

# Intrahost-level weights

function hostwideNetCoefficient(host::Host, population::Population, coef::Int64)
    net_hostwide_coefficient = 1.0
    if length(host.pathogens) > 0
        net_hostwide_coefficient = 0.0
        for p in 1:length(host.pathogens)
            net_hostwide_coefficient += host.pathogens[p].hostwide_coefficients[coef] * host.pathogen_fractions[p]
        end
    end
    for re in 1:length(host.responses)
        net_hostwide_coefficient = net_hostwide_coefficient * host.responses[re].hostwide_coefficients[coef]
    end

    return net_hostwide_coefficient * population.parameters.weightedInteractionHostwide(host, coef) * host.coefficients[coef]
end

function pathogenWeights!(p::Int64, host::Host, population::Population, evt::Int64)
    host.pathogen_weights[evt, p] = host.pathogens[p].specific_coefficients[evt]
    if evt == CLEARANCE || evt in PATHOGEN_NONSAMPLING_COEFFICIENTS
        # Clearance likelihoods are not proportional to intrahost fraction,
        # instead, we assume the rate of loss compounds
        # (could be inversely proportional to fraction?)
        # host.pathogen_weights[evt, p] = host.pathogen_weights[evt, p]
    elseif evt == RECOMBINANT_ESTABLISHMENT && length(host.pathogens) < 2
        # if nobody to recombine with, no recombination happens
        host.pathogen_weights[evt, p] = 0.0
        #TODO: regarding recombination and mutation establishment, each of these
        # depend on other parameters within Pathogen (mean_mutations_per_replication
        # and mean_recombination_crossovers); since in the current version we don't
        # consider intrahost population sizes or replication rates, then it's
        # pointless to incorporate these parameters in weight calculation here, but
        # eventually in a future version with better intrahost dynamics and popgen
        # this will be the case
    else
        # All other pathogen event likelihoods are proportional to intrahost fraction
        host.pathogen_weights[evt, p] =
            host.pathogen_weights[evt, p] *
            host.pathogen_fractions[p]
    end
    if length(host.responses) > 0
        host.pathogen_weights[evt, p] =
            host.pathogen_weights[evt, p] * population.parameters.weightedInteractionPathogen(
                host.pathogens[p], host, evt
            )
    end
end

function hostWeightsPathogen!(host_idx::Int64, population::Population, evt::Int64)
    population.host_weights[evt, host_idx] = 0.0
    for p in 1:length(population.hosts[host_idx].pathogens)
        pathogenWeights!(p, population.hosts[host_idx], population, evt)
        population.host_weights[evt, host_idx] += population.hosts[host_idx].pathogen_weights[evt, p]
        # we sum here because we have already weighted by fraction
        # (for everything except clearance, see above)
    end

    population.host_weights[evt, host_idx] = (
        population.host_weights[evt, host_idx] *
        hostwideNetCoefficient(population.hosts[host_idx], population, evt)
    )

    # population.host_weights_with_coefficient[evt, host_idx] = (
    #     population.host_weights[evt, host_idx] *
    #     population.parameters.base_coefficients[evt]
    # )
    setindex!(
        population.host_weights_with_coefficient_sampler[evt],
        population.host_weights[evt, host_idx] *
        population.parameters.base_coefficients[evt],
        host_idx
    )
end

function responseWeights!(re::Int64, host::Host, population::Population, evt::Int64)
    host.response_weights[evt-RESPONSE_EVENTS[1]+1, re] = host.responses[re].specific_coefficients[evt] * population.parameters.weightedInteractionResponse(
        host.responses[re], host, evt
    )
    # No Response fraction weighting here, just presence of Response is enough
end

function hostWeightsResponse!(host_idx::Int64, population::Population, evt::Int64)
    population.host_weights[evt, host_idx] = 0.0
    for re in 1:length(population.hosts[host_idx].responses)
        responseWeights!(re, population.hosts[host_idx], population, evt)
        population.host_weights[evt, host_idx] += population.hosts[host_idx].response_weights[evt-RESPONSE_EVENTS[1]+1, re]
        # We sum here because having more responses does imply more events, e.g. losses of response
    end

    population.host_weights[evt, host_idx] = (
        population.host_weights[evt, host_idx] *
        hostwideNetCoefficient(population.hosts[host_idx], population, evt)
    )

    # population.host_weights_with_coefficient[evt, host_idx] = (
    #     population.host_weights[evt, host_idx] *
    #     population.parameters.base_coefficients[evt]
    # )
    setindex!(
        population.host_weights_with_coefficient_sampler[evt],
        population.host_weights[evt, host_idx] *
        population.parameters.base_coefficients[evt],
        host_idx
    )
end

# Intra-Population level

function hostWeightsHost!(h::Int64, population::Population, evt::Int64)
    #TODO:
    # I think it might be okay here, but then we need to remove redundant user parameters that are
    # "host hostwide coefficients" in pathogens and responses
    population.host_weights[evt, h] = 1.0 * hostwideNetCoefficient(population.hosts[h], population, evt)
    # population.host_weights_with_coefficient[evt, h] = (
    #     population.host_weights[evt, h] *
    #     population.parameters.base_coefficients[evt]
    # )
    setindex!(
        population.host_weights_with_coefficient_sampler[evt],
        population.host_weights[evt, h] *
        population.parameters.base_coefficients[evt],
        h
    )
end

function hostWeightsReceive!(h::Int64, population::Population, evt::Int64)
    #TODO:
    # I think it might be okay here, but then we need to remove redundant user parameters that are
    # "host hostwide coefficients" in pathogens and responses
    population.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] = 1.0 * hostwideNetCoefficient(population.hosts[h], population, evt)
    # population.host_weights_receive_with_coefficient[evt-CHOICE_MODIFIERS[1]+1, h] = (
    #     population.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] *
    #     population.parameters.base_coefficients[evt]
    # )
    setindex!(
        population.host_weights_receive_with_coefficient_sampler[evt-CHOICE_MODIFIERS[1]+1],
        population.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] *
        population.parameters.base_coefficients[evt],
        h
    )
end

# #TODO:
# # I think it might be okay here, but then we need to remove redundant user parameters that are
# # "host hostwide coefficients" in pathogens and responses
# function hostWeightsNonsamplingPathogen!(host_idx::Int64, population::Population, evt::Int64)
#     for p in 1:length(population.hosts[host_idx].pathogens)
#         if length(host.responses) > 0
#             host.pathogen_weights[evt, p] =
#                 host.pathogen_weights[evt, p] * population.parameters.weightedInteractionPathogen(
#                     host.pathogens[p], host, evt
#                 )
#         end
#     end
#     population.hosts[host_idx].coefficients[evt] = hostwideNetCoefficient(population.hosts[host_idx], population, evt)
# end

# function hostWeightsNonsamplingResponse!(host_idx::Int64, population::Population, evt::Int64)
#     for re in 1:length(population.hosts[host_idx].responses)
#         responseWeights!(re, population.hosts[host_idx], population, evt)
#     end
#     population.hosts[host_idx].coefficients[evt] = hostwideNetCoefficient(population.hosts[host_idx], population, evt)
# end

# function hostWeightsNonsamplingHost!(host_idx::Int64, population::Population, evt::Int64)
#     population.hosts[host_idx].coefficients[evt] = hostwideNetCoefficient(population.hosts[host_idx], population, evt)
# end

function hostWeights!(host_idx::Int64, population::Population, model::Model)
    population.hosts[host_idx].pathogen_fractions = population.parameters.pathogenFractions(
        population.hosts[host_idx], population.parameters.weightedInteractionPathogen
    )
    prev = 0.0
    for weight in PATHOGEN_EVENTS
        prev = population.host_weights[weight, host_idx]
        hostWeightsPathogen!(host_idx, population, weight)
        if prev != population.host_weights[weight, host_idx] && population.parameters.base_coefficients[weight] != 0.0
            change = population.host_weights[weight, host_idx] - prev
            if weight == CONTACT
                population.contact_sum += change
                change = change * model.population_contact_weights_receive_sums[model.population_dict[population.id]]
            end
            propagateWeightChanges!(
                population.parameters.base_coefficients[weight] *
                change,
                population, weight, model
            )
        end
    end
    for weight in RESPONSE_EVENTS
        prev = population.host_weights[weight, host_idx]
        hostWeightsResponse!(host_idx, population, weight)
        if prev != population.host_weights[weight, host_idx] && population.parameters.base_coefficients[weight] != 0.0
            propagateWeightChanges!(
                population.parameters.base_coefficients[weight] *
                (population.host_weights[weight, host_idx] - prev),
                population, weight, model
            )
        end
    end
    for weight in HOST_EVENTS
        prev = population.host_weights[weight, host_idx]
        hostWeightsHost!(host_idx, population, weight)
        if prev != population.host_weights[weight, host_idx] && population.parameters.base_coefficients[weight] != 0.0
            change = population.host_weights[weight, host_idx] - prev
            if weight == TRANSITION
                population.transition_sum += change
                change = change * model.population_transition_weights_receive_sums[model.population_dict[population.id]]
            end
            propagateWeightChanges!(
                population.parameters.base_coefficients[weight] *
                change,
                population, weight, model
            )
        end
    end
    for weight in CHOICE_MODIFIERS[begin:end-1]
        prev = population.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, host_idx]
        hostWeightsReceive!(host_idx, population, weight)
        if (
            prev != population.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, host_idx] &&
            population.parameters.base_coefficients[weight] != 0.0 &&
            weight < INTRAHOST_FITNESS
        )

            propagateWeightReceiveChanges!(
                population.parameters.base_coefficients[weight] *
                (population.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, host_idx] - prev),
                population, weight, model
            )
        end
    end
    # for coef in PATHOGEN_NONSAMPLING_COEFFICIENTS
    #     hostWeightsNonsamplingPathogen!(host_idx, population, coef)
    # end
    # for coef in RESPONSE_NONSAMPLING_COEFFICIENTS
    #     hostWeightsNonsamplingResponse!(host_idx, population, coef)
    # end
    # for coef in HOST_NONSAMPLING_COEFFICIENTS
    #     hostWeightsNonsamplingHost!(host_idx, population, coef)
    # end
end

# Intra-Event/Weight level:

function calculateWeight!(population::Population, evt::Int64, model::Model)
    model.population_weights[evt, model.population_dict[population.id]] = sum([
        population.parameters.base_coefficients[evt] * population.host_weights[evt, h]
        for h in 1:length(population.hosts)
    ])
    population.host_weights_with_coefficient_sampler[evt].sum = model.population_weights[evt, model.population_dict[population.id]]
end

function calculateWeightReceive!(population::Population, weight::Int64, model::Model)
    model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]] = sum([
        population.parameters.base_coefficients[weight] * population.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, h]
        for h in 1:length(population.hosts)
    ])
    population.host_weights_receive_with_coefficient_sampler[weight-CHOICE_MODIFIERS[1]+1].sum = model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]]
end

function calculatePopulationContactWeightReceiveMatrix!(pop_idx_i::Int64, model::Model)
    model.population_contact_weights_receive_sums[pop_idx_i] = 0
    for pop_idx_j in values(model.population_dict)
        model.population_contact_weights_receive[pop_idx_j, pop_idx_i] = (
            model.populations[pop_idx_i].population_contact_coefficients[pop_idx_j] *
            model.population_weights_receive[RECEIVE_CONTACT-CHOICE_MODIFIERS[1]+1, pop_idx_j] /
            max(
                length(model.populations[pop_idx_j].hosts) *
                model.populations[pop_idx_j].parameters.constant_contact_density,
                1.0
            ))
        model.population_contact_weights_receive_sums[pop_idx_i] += model.population_contact_weights_receive[pop_idx_j, pop_idx_i]
    end
end

function calculatePopulationTransitionWeightReceiveMatrix!(pop_idx_i::Int64, model::Model)
    model.population_transition_weights_receive_sums[pop_idx_i] = 0
    for pop_idx_j in values(model.population_dict)
        model.population_transition_weights_receive[pop_idx_j, pop_idx_i] = (
            model.populations[pop_idx_i].population_transition_coefficients[pop_idx_j] *
            model.population_weights_receive[RECEIVE_TRANSITION-CHOICE_MODIFIERS[1]+1, pop_idx_j] /
            max(
                length(model.populations[pop_idx_j].hosts) *
                model.populations[pop_idx_j].parameters.constant_transition_density,
                1.0
            ))
        model.population_transition_weights_receive_sums[pop_idx_i] += model.population_transition_weights_receive[pop_idx_j, pop_idx_i]
    end
end

function updatePopulationContactWeightReceiveMatrix!(pop_idx_1::Int64, pop_idx_2::Int64, change::Float64, model::Model)
    model.population_contact_weights_receive[pop_idx_2, pop_idx_1] += change
    if becomesZero(model.population_contact_weights_receive_sums[pop_idx_1], change)
        prev = model.population_contact_weights_receive_sums[pop_idx_1]
        calculatePopulationContactWeightReceiveMatrix!(pop_idx_1, model)
        change = model.population_contact_weights_receive_sums[pop_idx_1] - prev
    end
    model.population_contact_weights_receive_sums[pop_idx_1] += change
    propagateWeightChanges!(
        model.populations[pop_idx_1].parameters.base_coefficients[CONTACT] *
        model.populations[pop_idx_1].contact_sum * change,
        model.populations[pop_idx_1], CONTACT, model
    )
end

function updatePopulationTransitionWeightReceiveMatrix!(pop_idx_1::Int64, pop_idx_2::Int64, change::Float64, model::Model)
    model.population_transition_weights_receive[pop_idx_2, pop_idx_1] += change
    if becomesZero(model.population_transition_weights_receive_sums[pop_idx_1], change)
        prev = model.population_transition_weights_receive_sums[pop_idx_1]
        calculatePopulationTransitionWeightReceiveMatrix!(pop_idx_1, model)
        change = model.population_transition_weights_receive_sums[pop_idx_1] - prev
    end
    model.population_transition_weights_receive_sums[pop_idx_1] += change
    propagateWeightChanges!(
        model.populations[pop_idx_1].parameters.base_coefficients[TRANSITION] *
        model.populations[pop_idx_1].transition_sum * change,
        model.populations[pop_idx_1], TRANSITION, model
    )
end

# Intra-model level:

function calculateRate!(evt::Int64, model::Model)
    model.event_rates[evt] = sum([
        model.population_weights[evt, pop_idx] for pop_idx in values(model.population_dict)
    ])
end

function calculateWeightReceive!(weight::Int64, model::Model)
    model.population_weights_receive_sums[weight-CHOICE_MODIFIERS[1]+1] = sum([
        model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, pop_idx] for pop_idx in values(model.population_dict)
    ])
end

# Weight change propagation:

function propagateWeightChanges!(change::Float64, evt::Int64, model::Model)
    if becomesZero(model.event_rates[evt], change)
        calculateRate!(evt, model)
    else
        model.event_rates[evt] += change
    end
    model.event_rates_sum = sum(model.event_rates) # this single line is all the "super-model level" calculation code haha
    # model.event_rates[evt] += change
    # model.event_rates_sum += change
end

function propagateWeightChanges!(change::Float64, population::Population, evt::Int64, model::Model)
    if change != 0.0
        if becomesZero(model.population_weights[evt, model.population_dict[population.id]], change)
            prev = model.population_weights[evt, model.population_dict[population.id]]
            calculateWeight!(population, evt, model)
            change = model.population_weights[evt, model.population_dict[population.id]] - prev
        else
            model.population_weights[evt, model.population_dict[population.id]] += change
        end
        propagateWeightChanges!(change, evt, model)
    end
end

function propagateWeightChanges!(change::Float64, host_idx::Int64, population::Population, evt::Int64, model::Model)
    if change != 0
        if length(population.hosts) > 0
            if evt == CONTACT
                population.contact_sum += change
                change = change * model.population_contact_weights_receive_sums[model.population_dict[population.id]]
            elseif evt == TRANSITION
                population.transition_sum += change
                change = change * model.population_transition_weights_receive_sums[model.population_dict[population.id]]
            end
        else
            if evt == CONTACT
                population.contact_sum = 0.0
            elseif evt == TRANSITION
                population.transition_sum = 0.0
            end
            change = -model.population_weights[evt, model.population_dict[population.id]]
        end

        population.host_weights[evt, host_idx] += change

        if population.parameters.base_coefficients[evt] != 0
            # population.host_weights_with_coefficient[evt, host_idx] += (
            #     change *
            #     population.parameters.base_coefficients[evt]
            # )
            setindex!(
                population.host_weights_with_coefficient_sampler[evt],
                population.host_weights[evt, host_idx] *
                population.parameters.base_coefficients[evt],
                host_idx
            )

            propagateWeightChanges!(population.parameters.base_coefficients[evt] * change, population, evt, model)
        end
    end
end

function propagateWeightReceiveChanges!(change::Float64, population::Population, evt::Int64, model::Model)
    if change != 0.0
        if becomesZero(model.population_weights_receive[evt-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]], change)
            prev = model.population_weights_receive[evt-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]]
            calculateWeightReceive!(population, evt, model)
            change = model.population_weights_receive[evt-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]] - prev
        else
            model.population_weights_receive[evt-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]] += change
        end
        if becomesZero(model.population_weights_receive_sums[evt-CHOICE_MODIFIERS[1]+1], change)
            prev = model.population_weights_receive_sums[evt-CHOICE_MODIFIERS[1]+1]
            calculateWeightReceive!(evt, model)
            change = model.population_weights_receive_sums[evt-CHOICE_MODIFIERS[1]+1] - prev
        else
            model.population_weights_receive_sums[evt-CHOICE_MODIFIERS[1]+1] += change
        end

        if evt == RECEIVE_CONTACT
            for p in 1:length(model.populations)
                change_p = (
                    change * model.populations[p].population_contact_coefficients[model.population_dict[population.id]] /
                    max(length(population.hosts) * population.parameters.constant_contact_density, 1.0)
                )
                # Contact rates assume scaling area if constant_contact_density is true:
                # large populations are equally
                # dense as small ones, so contact is constant (divide by total hosts).
                # If you don't want this to happen, modify each population's
                # receive contact coefficient accordingly.
                model.population_contact_weights_receive[model.population_dict[population.id], p] += change_p
                if approxZero(model.population_contact_weights_receive_sums[p] + change_p, t=ERROR_TOLERANCE)
                    change_p = -model.population_contact_weights_receive_sums[p]
                end
                model.population_contact_weights_receive_sums[p] += change_p
                propagateWeightChanges!(
                    change_p * model.populations[p].parameters.base_coefficients[CONTACT] * model.populations[p].contact_sum,
                    model.populations[p], CONTACT, model
                )
            end
        elseif evt == RECEIVE_TRANSITION
            for p in 1:length(model.populations)
                change_p = (
                    change * model.populations[p].population_transition_coefficients[model.population_dict[population.id]] /
                    max(length(population.hosts) * population.parameters.constant_transition_density, 1.0)
                )
                # Transition receive weights are assumed to be independent of population
                # (constant_transition_density=false), but can be modified if desired.
                model.population_transition_weights_receive[model.population_dict[population.id], p] += change_p
                model.population_transition_weights_receive_sums[p] += change_p
                propagateWeightChanges!(
                    change_p * model.populations[p].parameters.base_coefficients[TRANSITION] * model.populations[p].transition_sum,
                    model.populations[p], TRANSITION, model
                )
            end
        end
    end
end

function propagateWeightReceiveChanges!(change::Float64, host_idx::Int64, population::Population, evt::Int64, model::Model)
    if change != 0.0
        population.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, host_idx] += change
        # population.host_weights_receive_with_coefficient[evt-CHOICE_MODIFIERS[1]+1, host_idx] += (
        #     change *
        #     population.parameters.base_coefficients[evt]
        # )
        setindex!(
            population.host_weights_receive_with_coefficient_sampler[evt-CHOICE_MODIFIERS[1]+1],
            population.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, host_idx] *
            population.parameters.base_coefficients[evt],
            host_idx
        )

        if evt < INTRAHOST_FITNESS && population.parameters.base_coefficients[evt] != 0.0
            propagateWeightReceiveChanges!(population.parameters.base_coefficients[evt] * change, population, evt, model)
        end
    end
end

function propagateWeightsOnAddHost!(host_num::Int64, population::Population, model::Model)
    for p in 1:length(model.populations)
        model.population_contact_weights_receive[model.population_dict[population.id], p] -= (
            model.population_contact_weights_receive[model.population_dict[population.id], p] /
            max(host_num * population.parameters.constant_contact_density, 1)
        )
        model.population_contact_weights_receive_sums[p] -= (
            model.population_contact_weights_receive_sums[p] /
            max(host_num * population.parameters.constant_contact_density, 1)
        )
        propagateWeightChanges!(
            -model.population_weights[CONTACT, p] /
            max(host_num * population.parameters.constant_contact_density, 1),
            model.populations[p], CONTACT, model
        )
        model.population_transition_weights_receive[model.population_dict[population.id], p] -= (
            model.population_transition_weights_receive[model.population_dict[population.id], p] /
            max(host_num * population.parameters.constant_transition_density, 1)
        )
        model.population_transition_weights_receive_sums[p] -= (
            model.population_transition_weights_receive_sums[p] /
            max(host_num * population.parameters.constant_transition_density, 1)
        )
        propagateWeightChanges!(
            -model.population_weights[TRANSITION, p] /
            max(host_num * population.parameters.constant_transition_density, 1),
            model.populations[p], TRANSITION, model
        )
    end

    for coef in EVENTS
        if START_COEFFICIENTS[coef] != 0.0
            propagateWeightChanges!(
                START_COEFFICIENTS[coef], host_num, population, coef, model
            )
        end
    end
    for coef in CHOICE_MODIFIERS[begin:end-1]
        if START_COEFFICIENTS[coef] != 0.0
            propagateWeightReceiveChanges!(
                START_COEFFICIENTS[coef], host_num, population, coef, model
            )
        end
    end
end
