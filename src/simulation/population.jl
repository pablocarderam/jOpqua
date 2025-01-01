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

function setRates!(population::Population)
    classWeights!(population)
    classWeightsReceive!(population)
end

function newPopulation!(id::String, model::Model)
    population = Population(
        id,
        Vector{Class}(undef, 0),
        Dict{String,Int64}(),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 2, 0),
    )
    push!(model.populations, population)
end
