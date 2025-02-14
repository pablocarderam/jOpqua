function randChoose(rand_n::Float64, rates::AbstractVector{Float64}, rates_sum::Float64; regenerate_rand::Bool=false)
    #TODO: rates is not type-defined to accommodate @views, maybe this is bad
    if rates_sum > 0.0
        norm_rand_n = rand_n * rates_sum
        cum_sum = 0.0
        idx = 0
        chosen_rate = 0.0
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
    else
        return NaN, rand_n
    end
end

function rejectionSample(rates::AbstractVector{Float64}, r_max::Float64)
    indices = eachindex(rates)
    while true
        i = rand(indices)
        if rates[i] > rand() * r_max
            return i
        end
    end
end

function choosePathogen(host_idx::Int64, population_idx::Int64, weight::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.populations[population_idx].hosts[host_idx].pathogen_weights[weight, :]),
        model.populations[population_idx].host_weights[weight, host_idx],
        # assumes population's host weights are the sum of all pathogen weights, should be true
        regenerate_rand=true
    )
end

function chooseResponse(host_idx::Int64, population_idx::Int64, weight::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.populations[population_idx].hosts[host_idx].response_weights[weight-RESPONSE_EVENTS[1]+1, :]),
        model.populations[population_idx].host_weights[weight, host_idx],
        # assumes population's host weights are the sum of all response weights, should be true
        regenerate_rand=true
    )
end

function chooseHost(population_idx::Int64, weight::Int64, model::Model, rand_n::Float64)
    if weight > NUM_EVENTS
        return randChoose(
            rand_n,
            @views(model.populations[population_idx].host_weights_receive_with_coefficient[weight-CHOICE_MODIFIERS[1]+1, :]),
            model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, population_idx],
            #TODO: revise, maybe this is another field?
            regenerate_rand=true
        )
    else
        return randChoose(
            rand_n,
            @views(model.populations[population_idx].host_weights_with_coefficient[weight, :]),
            model.population_weights[weight, population_idx],
            #TODO: revise, maybe this is another field?
            regenerate_rand=true
        )
    end
end

function choosePopulation(weight::Int64, model::Model, rand_n::Float64)
    if weight > NUM_EVENTS
        return randChoose(
            rand_n,
            @views(model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, :]),
            model.population_weights_receive_sums[weight-CHOICE_MODIFIERS[1]+1],
            regenerate_rand=true
        )
    else
        return randChoose(
            rand_n,
            @views(model.population_weights[weight, :]),
            model.event_rates[weight],
            regenerate_rand=true
        )
    end
end

function choosePopulationReceiveContact(emitter_idx::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.population_contact_weights_receive[:, emitter_idx]),
        model.population_contact_weights_receive_sums[emitter_idx],
        regenerate_rand=true
    )
end

function choosePopulationReceiveTransition(emitter_idx::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.population_transition_weights_receive[:, emitter_idx]),
        model.population_transition_weights_receive_sums[emitter_idx],
        regenerate_rand=true
    )
end

function chooseEvent(model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        model.event_rates,
        model.event_rates_sum,
        regenerate_rand=true
    )
end

function choosePathogen(weight::Int64, model::Model, rand_n::Float64)
    pop_idx, rand_n = choosePopulation(weight, model, rand_n)
    host_idx, rand_n = chooseHost(pop_idx, weight, model, rand_n)
    pathogen_idx, rand_n = choosePathogen(host_idx, pop_idx, weight, model, rand_n)

    return pathogen_idx, host_idx, pop_idx, rand_n
end

function chooseResponse(weight::Int64, model::Model, rand_n::Float64)
    pop_idx, rand_n = choosePopulation(weight, model, rand_n)
    host_idx, rand_n = chooseHost(pop_idx, weight, model, rand_n)
    response_idx, rand_n = chooseResponse(host_idx, pop_idx, weight, model, rand_n)

    return response_idx, host_idx, pop_idx, rand_n
end

function chooseHost(weight::Int64, model::Model, rand_n::Float64)
    pop_idx, rand_n = choosePopulation(weight, model, rand_n)
    host_idx, rand_n = chooseHost(pop_idx, weight, model, rand_n)

    return host_idx, pop_idx, rand_n
end

function newFlexLevel(i::Int64, w::Float64; next::Union{FlexLevel, Nothing}=nothing)
    return FlexLevel(logBounds(w), w, w, [i], next)
end

function addToFlexLevel!(i::Int64, w::Float64, level::FlexLevel)
    push!(level.indices, i)
    level.sum += w
    if w > level.max
        level.max = w
    end
end

function verifyFlexlevSampler(sampler::FlexlevSampler)
    errors = 0
    
    # ]all weights in sampler.weights are represented in some level?
    for i in eachindex(sampler.weights)
        if !inSampler(i, sampler)
            errors += 1
            @printf "Error %i: index %i (weight %f) not in any level\n" errors i sampler.weights[i]
        end
    end

    # all weights in levels are represented in sampler.weights, and if so, in the correct level?
    curr = sampler.levels
    num_weights = length(sampler.weights)
    while !isnothing(curr)
        for i in curr.indices
            if !(1 <= i <= num_weights)
                errors += 1
                @printf "Error %i: index %i present in level (%f, %f) but not in weights\n" errors i curr.bounds[1] curr.bounds[2]
            elseif !(curr.bounds[1] <= sampler.weights[i] < curr.bounds[2])
                errors += 1
                @printf "Error %i: index %i (weight %f) present in level (%f, %f)\n" errors i sampler.weights[i] curr.bounds[1] curr.bounds[2]
            end
        end
        curr = curr.next
    end

    # all level sums/maxes correct?
    curr = sampler.levels
    while !isnothing(curr)
        expected_sum = sum(sampler.weights[i] for i in curr.indices)
        if expected_sum != curr.sum
            errors += 1
            @printf "Error %i: level (%f, %f) sum incorrect (expected %f, got %f)\n" errors curr.bounds[1] curr.bounds[2] expected_sum curr.sum
        end
        expected_max = maximum(sampler.weights[i] for i in curr.indices)
        if expected_max != curr.max
            errors += 1
            @printf "Error %i: level (%f, %f) max incorrect (expected %f, got %f)\n" errors curr.bounds[1] curr.bounds[2] expected_max curr.max
        end
        curr = curr.next
    end

    @printf"Final error count: %i" errors
    return errors
end

function inSampler(i::Int64, sampler::FlexlevSampler)
    curr = sampler.levels
    while !isnothing(curr)
        if i in curr.indices
            return true
        end
        curr = curr.next
    end
    return false
end

function printFlexlevSampler(sampler::FlexlevSampler)
    @printf "\nFlexlevSampler (%f, %f)\n" sampler.min_level sampler.max_level
    show(sampler.weights)
    println()
    curr = sampler.levels
    while !isnothing(curr)
        @printf "\nLevel (%f, %f)\n" curr.bounds[1] curr.bounds[2]
        @printf "sum=%f\n" curr.sum
        @printf "max=%f\n" curr.max
        println(curr.indices)
        curr = curr.next
    end
end
