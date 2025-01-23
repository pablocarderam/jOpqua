using StaticArrays

function newPathogen!(sequence::String, class::Class, type::PathogenType)
    class.pathogens[sequence] = Pathogen(
        sequence, pathogenSequenceCoefficients(sequence, type),
        type.mean_effective_inoculum * type.inoculumCoefficient(sequence),
        type.mean_mutations_per_replication * type.mutationCoefficient(sequence),
        type.mean_recombination_crossovers * type.recombinationCoefficient(sequence),
        type
    )

    return class.pathogens[sequence]
end

function newResponse!(
    imprinted_pathogen::Pathogen, matured_pathogen::Pathogen, parent::Tuple{String,String,String},
    class::Class, type::ResponseType)
    class.responses[(imprinted_pathogen.sequence, matured_pathogen.sequence, type.id)] = Response(
        parent, imprinted_pathogen, matured_pathogen,
        responseStaticCoefficients(imprinted_pathogen.sequence, matured_pathogen.sequence, type),
        type
    )

    return class.responses[(imprinted_pathogen.sequence, matured_pathogen.sequence, type.id)]
end

function newHost!(class::Class, population::Population, model::Model)
    push!(class.hosts, Host(
        length(class.hosts) + 1,
        Vector{Pathogen}(undef, 0), Vector{Response}(undef, 0),
        Vector{Float64}(undef, 0),
        Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, 0),
        Matrix{Float64}(undef, NUM_RESPONSE_EVENTS, 0),
    ))
    class.host_weights = catCol(class.host_weights, zeros(Float64, NUM_EVENTS))
    class.host_weights_receive = catCol(class.host_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS - 1))

    model.population_weights[INTRA_POPULATION_CONTACT, model.populations_dict[population.id]] = (
        model.population_weights[INTRA_POPULATION_CONTACT, model.populations_dict[population.id]] *
        population.total_hosts / (population.total_hosts + 1)
    )
    #TODO: update inter-pop contacts
    population.total_hosts += 1

    propagateWeightChanges!(
        START_COEFFICIENTS
        # class change, inter-population contact and migration numbers per class or
        # population are fractions that sum to one, so no need to account for in here
        , length(class.hosts), class, population, model
    )

    return class.hosts[end]
end

function newClass!(id::String, parameters::ClassParameters, population::Population)
    push!(population.classes, Class(
        id, parameters,
        Dict{String,Pathogen}(), Dict{Tuple{String,String,String},Response}(),
        Vector{Host}(undef, 0),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
    ))
    population.class_dict[id] = length(population.classes)
    population.class_weights = catCol(population.class_weights, zeros(Float64, NUM_EVENTS))
    population.class_weights_receive = catCol(population.class_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS - 1))

    return population.classes[end]
end

function newPopulation!(id::String, model::Model)
    push!(model.populations, Population(
        id,
        Vector{Class}(undef, 0),
        Dict{String,Int64}(),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
        0, 0.0
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
