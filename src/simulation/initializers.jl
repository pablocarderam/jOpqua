using StaticArrays

function newPathogen!(sequence::String, population::Population, type::PathogenType)
    population.pathogens[sequence] = Pathogen(
        sequence, pathogenSequenceCoefficients(sequence, type),
        type.mean_effective_inoculum * type.inoculumCoefficient(sequence),
        type.mean_mutations_per_replication * type.mutationCoefficient(sequence),
        type.mean_recombination_crossovers * type.recombinationCoefficient(sequence),
        type
    )

    return population.pathogens[sequence]
end

function newResponse!(
    imprinted_pathogen::Pathogen, matured_pathogen::Pathogen, parent::Tuple{String,String,String},
    population::Population, type::ResponseType)
    population.responses[(imprinted_pathogen.sequence, matured_pathogen.sequence, type.id)] = Response(
        parent, imprinted_pathogen, matured_pathogen,
        responseStaticCoefficients(imprinted_pathogen.sequence, matured_pathogen.sequence, type),
        type
    )

    return population.responses[(imprinted_pathogen.sequence, matured_pathogen.sequence, type.id)]
end

function newHost!(population::Population, model::Model)
    push!(population.hosts, Host(
        length(population.hosts) + 1,
        Vector{Pathogen}(undef, 0), Vector{Response}(undef, 0),
        Vector{Float64}(undef, 0),
        Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, 0),
        Matrix{Float64}(undef, NUM_RESPONSE_EVENTS, 0),
    ))
    population.host_weights = catCol(population.host_weights, zeros(Float64, NUM_EVENTS))
    population.host_weights_receive = catCol(population.host_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS - 1))

    propagateWeightChanges!(
        - model.population_weights[CONTACT, model.population_dict[population.id]] /
        (population.total_hosts + 1),
        population, CONTACT, model
    )
    #TODO: update inter-pop contacts
    population.total_hosts += 1

    for coef in EVENTS
        if START_COEFFICIENTS[coef] != 0
            propagateWeightChanges!(
                START_COEFFICIENTS[coef], length(population.hosts), population, coef, model
            )
        end
    end
    for coef in CHOICE_MODIFIERS[begin:end-1]
        if START_COEFFICIENTS[coef] != 0
            propagateWeightReceiveChanges!(
                START_COEFFICIENTS[coef], length(population.hosts), population, coef, model
            )
        end
    end

    # propagateWeightChanges!(
    #     START_COEFFICIENTS
    # )

    return population.hosts[end]
end

function newPopulation!(id::String, parameters::PopulationType, model::Model)
    push!(model.populations, Population(
        id, parameters,
        Dict{String,Pathogen}(), Dict{Tuple{String,String,String},Response}(),
        Vector{Host}(undef, 0),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
        0, 0.0, 0.0
    ))
    model.population_dict[id] = length(model.populations)
    model.population_weights = catCol(model.population_weights, zeros(Float64, NUM_EVENTS))
    model.population_weights_receive = catCol(model.population_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS - 1))

    return model.populations[end]
end

function newModel()
    return Model(
        Vector{Population}(undef, 0),
        Dict{String,Int64}(),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
        zeros(SVector{NUM_CHOICE_MODIFIERS - 1,Float64}),
        zeros(SVector{NUM_EVENTS,Float64}),
        0.0
    )
end
