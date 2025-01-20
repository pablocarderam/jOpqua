function randChoose(rand_n::Float64, rates, rates_sum::Float64; regenerate_rand::Bool=false)
        #TODO: rates is not type defined to accommodate @views, maybe this is bad
    norm_rand_n = rand_n * rates_sum
    cum_sum = 0.0
    idx = 0
    chosen_rate = 0
    for r in rates
        idx += 1
        cum_sum += r
        chosen_rate = r
        if cum_sum > norm_rand_n
            break
        end
    end

    if regenerate_rand
        return idx, (norm_rand_n - cum_sum + chosen_rate) / chosen_rate
    else
        return idx, 0.0
    end
end

function choosePathogen(host_idx::Int64, class_idx::Int64, population_idx::Int64, weight::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.populations[population_idx].classes[class_idx].hosts[host_idx].pathogen_weights[weight,:]),
        model.populations[population_idx].classes[class_idx].host_weights[weight, host_idx],
            # assumes class's host weights are the sum of all pathogen weights, should be true
        regenerate_rand=true
    )
end

function chooseResponse(host_idx::Int64, class_idx::Int64, population_idx::Int64, weight::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.populations[population_idx].classes[class_idx].hosts[host_idx].response_weights[weight-RESPONSE_EVENTS[1]+1,:]),
        model.populations[population_idx].classes[class_idx].host_weights[weight, host_idx],
            # assumes class's host weights are the sum of all response weights, should be true
        regenerate_rand=true
    )
end

function chooseHost(class_idx::Int64, population_idx::Int64, weight::Int64, model::Model, rand_n::Float64)
    if weight > NUM_EVENTS
        println(weight)
        return randChoose(
            rand_n,
            @views(model.populations[population_idx].classes[class_idx].host_weights_receive[weight-CHOICE_MODIFIERS[1]+1,:]) .*
                model.populations[population_idx].classes[class_idx].parameters.base_coefficients[weight],
            model.populations[population_idx].class_weights_receive[weight-CHOICE_MODIFIERS[1]+1, class_idx],
            regenerate_rand=true
        )
    else
        return randChoose(
            rand_n,
            @views(model.populations[population_idx].classes[class_idx].host_weights[weight,:]) .*
                model.populations[population_idx].classes[class_idx].parameters.base_coefficients[weight],
            model.populations[population_idx].class_weights[weight, class_idx],
            regenerate_rand=true
        )
    end
end

function chooseClass(population_idx::Int64, weight::Int64, model::Model, rand_n::Float64)
    if weight > NUM_EVENTS
        return randChoose(
            rand_n,
            @views(model.populations[population_idx].class_weights_receive[weight-CHOICE_MODIFIERS[1]+1,:]),
            model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, population_idx],
            regenerate_rand=true
        )
    else
        return randChoose(
            rand_n,
            @views(model.populations[population_idx].class_weights[weight,:]),
            model.population_weights[weight, population_idx],
            regenerate_rand=true
        )
    end
end

function choosePopulation(weight::Int64, model::Model, rand_n::Float64)
    if weight > NUM_EVENTS
        return randChoose(
                rand_n,
                @views(model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1,:]),
                model.population_weights_receive_sums[weight-CHOICE_MODIFIERS[1]+1],
                regenerate_rand=true
            )
    else
        return randChoose(
            rand_n,
            @views(model.population_weights[weight,:]),
            model.event_rates[weight],
            regenerate_rand=true
        )
    end
end

function chooseEvent(model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        model.event_rates,
        model.event_rates_sum,
        regenerate_rand=true
    )
end
