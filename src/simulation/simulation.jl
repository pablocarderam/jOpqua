using Random

function simulate!(model::Model, time_vector::Vector{Float64})
    # Gillespie
    time = 0.0
    dt = 0.0
    rand_n = 0.0
    # t_len = length(time_vector)
    evt_idx = 0
    evt_count = 0
    while time < time_vector[end]
        if (model.event_rates_sum > 0.0)
            dt = randexp() / model.event_rates_sum
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
            event_functions[evt_idx](model, rand_n)

            # alternative sampling method:
            # evt_idx = sample(1:length(evt_funcs), Weights(model.event_rates))
            # event_functions[evt_idx](model, rand_n)

            evt_count += 1
        else
            time = time_vector[end]
        end
    end

    println(evt_count)

    # return model
end
