using Random

function simulate!(
        model::Model, time_vector::Vector{Float64};
        host_samples_population::Dict{String, Int64}=Dict{String, Int64}())

    # Interventions
    interventions = sort(model.interventions, by=intervention -> intervention.time)
    intervention_tracker = 0 # keeps track of what the next intervention should be

    # History
    his_tracker = 1
    compartment_vars = Dict{String, Matrix{Int64}}()
    for p in model.populations
        compartment_vars[p.id] = Matrix{Int64}(
            undef, NUM_COMPARTMENTS, length(time_vector)
        )
        # tracks uninfected naive, infected naive, uninfected immune, infected immune,
        # and dead for this population at each time point
    end
    host_samples = Dict{String, Matrix{StaticHost}}()
    for p in model.populations
        if haskey(host_samples_population, p.id)
            host_samples[p.id] = Matrix{StaticHost}(
                undef, host_samples_population[p.id], length(time_vector)
            )
        end
    end

    # Gillespie
    time = 0.0
    dt = 0.0
    rand_n = 0.0
    evt_idx = 0
    evt_count = 0
    while time < time_vector[end]
        if (model.event_rates_sum > 0.0)
            dt = randexp() / model.event_rates_sum
            if (intervention_tracker < length(interventions) &&
                time + dt >= interventions[intervention_tracker].time)
                # if there are any interventions left and if it is time to make one,
                time = interventions[intervention_tracker].time
                interventions[intervention_tracker].intervention(model)
                intervention_tracker += 1
            else
                time += dt

                if time > time_vector[end]
                    time = time_vector[end]
                end

                rand_n = rand()
                evt_idx, rand_n = randChoose(
                    rand_n, model.event_rates,
                    model.event_rates_sum, regenerate_rand=true
                )
                # println((time, model.event_rates, evt_idx))
                EVENT_FUNCTIONS[evt_idx](model, rand_n)

                # alternative sampling method:
                # evt_idx = sample(1:length(evt_funcs), Weights(model.event_rates))
                # EVENT_FUNCTIONS[evt_idx](model, rand_n)

                evt_count += 1
            end
        else
            if (intervention_tracker < length(interventions) &&
                time >= interventions[intervention_tracker].time)
                # if there are any interventions left
                time = interventions[intervention_tracker].time
                interventions[intervention_tracker].intervention(model)
                intervention_tracker += 1
            else
                time = time_vector[end]
            end
        end
        while his_tracker <= length(time_vector) && time >= time_vector[his_tracker]
            for p in model.populations
                compartment_vars[p.id][:,his_tracker] = p.compartment_vars
                if haskey(host_samples_population, p.id)
                    for i in 1:min(host_samples_population[p.id], length(p.hosts))
                        host_samples[p.id][i,his_tracker] = staticHost(p.hosts[i])
                    end
                end
            end
            his_tracker += 1
        end
    end

    println(evt_count)

    return model, compartment_vars, host_samples
end
