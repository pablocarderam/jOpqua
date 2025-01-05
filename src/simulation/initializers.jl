using StaticArrays

function newPathogen!(sequence::String, class::Class)
    # class.pathogen_count += 1
    push!(class.pathogens,Pathogen(
        length(class.pathogens) + 1, sequence, pathogenSequenceCoefficients(sequence, class)
    ))

    return class.pathogens[end]
end

function newImmunity!(imprinted_pathogen::Pathogen, matured_pathogen::Pathogen, class::Class, type::ImmunityType)
    # class.immunity_count += 1
    push!(class.immunities, Immunity(
        length(class.immunities) + 1, imprinted_pathogen.id, matured_pathogen.id,
        immunityStaticCoefficients(imprinted_pathogen.sequence, matured_pathogen.sequence, type, class),
        type.id
    ))

    return class.immunities[end]
end

function newHost!(class::Class,population::Population,model::Model)
    push!(class.hosts, Host(
        length(class.hosts) + 1,
        Vector{Int64}(undef, 0), Vector{Int64}(undef, 0),
        Vector{Float64}(undef, 0),
        Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, 0),
        Matrix{Float64}(undef, NUM_IMMUNITY_EVENTS, 0),
    ))
    class.host_weights = catCol(
        class.host_weights,
        [
            0.0,0.0,0.0,0.0,0.0,
            0.0,1.0,1.0,1.0,1.0,
            # class change, inter-population contact and migration numbers per class or
            # population are fractions that sum to one, so no need to account for in here
        ] .* class.parameters.base_coefficients[begin:NUM_EVENTS]
    )
    class.host_weights_receive = catCol(
        class.host_weights_receive,
        [0.0,0.0,0.0,1.0] .* class.parameters.base_coefficients[end-NUM_CHOICE_MODIFIERS+1:end-1]
    )

    propagateWeightChanges!(SVector{NUM_COEFFICIENTS,Float64}([class.host_weights[:,length(class.hosts)]...,class.host_weights_receive[:,length(class.hosts)]...,0.0]), length(class.hosts), class, population, model)

    return class.hosts[end]
end

function newClass!(id::String, parameters::ClassParameters, population::Population)
    push!(population.classes, Class(
        id, parameters,
        Vector{Pathogen}(undef, 0), Vector{Immunity}(undef, 0),
        Vector{Host}(undef, 0),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
    ))
    population.class_dict[id] = length(population.classes)
    population.class_weights = catCol(population.class_weights, zeros(Float64, NUM_EVENTS))
    population.class_weights_receive = catCol(population.class_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS - 2))

    return population.classes[end]
end

function newPopulation!(id::String, model::Model)
    push!(model.populations, Population(
        id,
        Vector{Class}(undef, 0),
        Dict{String,Int64}(),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 2, 0),
    ))
    model.population_dict[id] = length(model.populations)
    model.population_weights = catCol(model.population_weights, zeros(Float64, NUM_EVENTS))
    model.population_weights_receive = catCol(model.population_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS - 3))

    return model.populations[end]
end

function newModel()
    return Model(
        Vector{Population}(undef, 0),
        Dict{String,Int64}(),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 3, 0),
        zeros(SVector{jOpqua.NUM_EVENTS,Float64}),
    )
end
