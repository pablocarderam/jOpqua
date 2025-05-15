using Random

function simulate!(
    model::Model, time_vector::Vector{Float64};
    population_host_samples::Dict{String,Int64}=Dict{String,Int64}(),
    interventions::Vector{Intervention}=Vector{Intervention}(undef, 0))

    # Interventions
    interventions = sort(interventions, by=intervention -> intervention.time)
    intervention_tracker = 0 # keeps track of what the next intervention should be

    # History
    his_tracker = 1
    compartment_vars = Dict{String,Matrix{Int64}}()
    for p in model.populations
        compartment_vars[p.id] = Matrix{Int64}(
            undef, NUM_COMPARTMENTS, length(time_vector)
        )
        # tracks uninfected naive, infected naive, uninfected immune, infected immune,
        # and dead for this population at each time point
    end
    host_samples_idxs = Dict{String,Vector{Int64}}() # this contains indexes of hosts to be sampled
    max_sample_idx = Dict{String,Int64}() # this is used to track max of each vector in host_samples_idxs
    host_samples = Dict{String,Matrix{StaticHost}}() # this stores history
    for p in model.populations
        if haskey(population_host_samples, p.id)
            if isa(population_host_samples[p.id], Number)
                host_samples_idxs[p.id] = sort(rand(1:length(p.hosts), population_host_samples[p.id]))
            else
                host_samples_idxs[p.id] = population_host_samples[p.id]
            end
            max_sample_idx[p.id] = maximum(host_samples_idxs[p.id])
            host_samples[p.id] = Matrix{StaticHost}(
                undef, length(host_samples_idxs[p.id]), length(time_vector)
            )
        end
    end

    # Gillespie
    dt = 0.0
    rand_n = 0.0
    evt_idx = 0
    evt_count = 0
    while model.time < time_vector[end]
        if (model.event_rates_sum > 0.0)
            dt = randexp() / model.event_rates_sum
            if (intervention_tracker < length(interventions) &&
                model.time + dt >= interventions[intervention_tracker].time)
                # if there are any interventions left and if it is time to make one,
                model.time = interventions[intervention_tracker].time
                interventions[intervention_tracker].intervention(model)
                intervention_tracker += 1
            else
                model.time += dt

                if model.time > time_vector[end]
                    model.time = time_vector[end]
                end

                rand_n = rand()
                evt_idx, rand_n = randChoose(
                    rand_n, model.event_rates,
                    model.event_rates_sum, regenerate_rand=true
                )
                # println((model.time, model.event_rates, evt_idx))
                EVENT_FUNCTIONS[evt_idx](model, rand_n)

                # alternative sampling method:
                # evt_idx = sample(1:length(evt_funcs), Weights(model.event_rates))
                # EVENT_FUNCTIONS[evt_idx](model, rand_n)

                evt_count += 1
            end
        else
            if (intervention_tracker < length(interventions) &&
                model.time >= interventions[intervention_tracker].time)
                # if there are any interventions left
                model.time = interventions[intervention_tracker].time
                interventions[intervention_tracker].intervention(model)
                intervention_tracker += 1
            else
                model.time = time_vector[end]
            end
        end
        while his_tracker <= length(time_vector) && model.time >= time_vector[his_tracker]
            for p in model.populations
                compartment_vars[p.id][:, his_tracker] = p.compartment_vars
                if haskey(host_samples_idxs, p.id)
                    if max_sample_idx[p.id] > length(p.hosts)
                        host_samples_idxs[p.id] = sort(
                            rand(1:length(p.hosts), length(host_samples_idxs[p.id]))
                        )
                    end
                    for (i, host_idx) in enumerate(host_samples_idxs[p.id])
                        host_samples[p.id][i, his_tracker] = staticHost(p.hosts[host_idx])
                    end
                end
            end
            his_tracker += 1
            # println((model.time, evt_idx))
            # println((model.event_rates, model.populations[1].compartment_vars, model.population_weights_receive, sum(model.populations[1].host_weights[CONTACT,:]), sum(model.populations[1].host_weights_with_coefficient[CONTACT,:]), model.population_contact_weights_receive_sums, model.populations[1].contact_sum, 1.05*model.population_contact_weights_receive_sums*model.populations[1].contact_sum))
        end
    end

    println(evt_count)

    return model, Output(model, time_vector, compartment_vars, host_samples)
end
