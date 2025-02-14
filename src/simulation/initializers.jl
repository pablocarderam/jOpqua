using StaticArrays

function newPathogen!(
        sequence::String, population::Population, type::PathogenType;
        parents::MVector{2, Union{Pathogen, Nothing}}=MVector{2, Union{Pathogen, Nothing}}([nothing, nothing]))
    population.pathogens[sequence] = Pathogen(
        parents, sequence, pathogenSequenceCoefficients(sequence, type),
        type.mean_effective_inoculum * type.inoculumCoefficient(sequence) * population.parameters.inoculum_coefficient,
        type.mean_mutations_per_replication * type.mutationCoefficient(sequence) * population.parameters.mutation_coefficient,
        type.mean_recombination_crossovers * type.recombinationCoefficient(sequence) * population.parameters.recombination_coefficient,
        type
    )

    return population.pathogens[sequence]
end

function newResponse!(
    imprinted_pathogen::Pathogen, matured_pathogen::Pathogen,
    population::Population, type::ResponseType;
    parents::MVector{2, Union{Response, Nothing}}=MVector{2, Union{Response, Nothing}}([nothing, nothing]))
    population.responses[(imprinted_pathogen.sequence, matured_pathogen.sequence, type.id)] = Response(
        parents, imprinted_pathogen, matured_pathogen,
        responseStaticCoefficients(imprinted_pathogen.sequence, matured_pathogen.sequence, type),
        type
    )

    return population.responses[(imprinted_pathogen.sequence, matured_pathogen.sequence, type.id)]
end

function newHost!(population::Population, model::Model)
    addHostToPopulation!(
        Host(
            length(population.hosts) + 1,
            Vector{Pathogen}(undef, 0), Vector{Response}(undef, 0),
            Vector{Float64}(undef, 0),
            Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, 0),
            Matrix{Float64}(undef, NUM_RESPONSE_EVENTS, 0),
        ),
        population, model
    )

    return population.hosts[end]
end

function newPopulation!(id::String, parameters::PopulationType, model::Model)
    push!(model.populations, Population(
        id, parameters,
        Dict{String,Pathogen}(), Dict{Tuple{String,String,String},Response}(),
        Vector{Host}(undef, 0),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
        0.0, 0.0,
        zeros(Float64, length(model.populations)),
        zeros(Float64, length(model.populations)),
    ))
    model.population_dict[id] = length(model.populations)
    model.population_weights = catCol(model.population_weights, zeros(Float64, NUM_EVENTS))
    model.population_weights_receive = catCol(model.population_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS - 1))

    for pop in model.populations
        push!(pop.population_contact_coefficients, 0.0)
        push!(pop.population_transition_coefficients, 0.0)
    end

    if length(model.populations) > 1 # if more than just this population we just added
        model.population_contact_weights_receive = catCol(
            model.population_contact_weights_receive, zeros(Float64, length(model.populations))
        )
        model.population_transition_weights_receive = catCol(
            model.population_transition_weights_receive, zeros(Float64, length(model.populations))
        )
    else
        model.population_contact_weights_receive = Matrix{Float64}(undef, 1, 1)
        model.population_transition_weights_receive = Matrix{Float64}(undef, 1, 1)
    end

    push!(model.population_contact_weights_receive_sums, 0.0)
    push!(model.population_transition_weights_receive_sums, 0.0)

    # By default, self contact is set to 1.0
    setPopulationContactCoefficient!(length(model.populations), length(model.populations), 1.0, model)

    return model.populations[end]
end

function newModel()
    return Model(
        Vector{Population}(undef, 0),
        Dict{String,Int64}(),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
        zeros(SVector{NUM_CHOICE_MODIFIERS - 1,Float64}),
        Matrix{Float64}(undef, 0, 0),
        Matrix{Float64}(undef, 0, 0),
        Vector{Float64}(undef, 0),
        Vector{Float64}(undef, 0),
        Vector{Intervention}(undef, 0),
        zeros(SVector{NUM_EVENTS,Float64}),
        0.0
    )
end
