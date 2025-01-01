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

function rates!(model::Model)
    for evt in EVENTS
        model.event_rates[evt] = sum([
            model.population_weights[evt, p]
            #TODO: This loop has to be a switch with a case for every event type computing the rate the correct way
            for p in 1:length(model.populations)
        ])
    end
end

function setRates!(model::Model)
    populationWeights!(model)
    populationWeightsReceive!(model)
    rates!(model)
end

function newPopulation(id::String)
    model = Model(
        id,
        Vector{Population}(undef, 0),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 3, 0),
        Vector{Float64}(undef, 0),
    )
    return model
end
