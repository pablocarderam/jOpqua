using StaticArrays
using FunctionWrappers
import FunctionWrappers: FunctionWrapper

# Model parameter setup initializers

function newPathogenType(
        id::String;

        template::PathogenType=DEFAULT_PATHOGEN_TYPE,

        num_loci::Union{Nothing,Int64}=nothing,
        possible_alleles::Union{Nothing,String}=nothing,
        mean_effective_inoculum::Union{Nothing,Float64}=nothing,
        mean_mutations_per_replication::Union{Nothing,Float64}=nothing,
        mean_recombination_crossovers::Union{Nothing,Float64}=nothing,
        vertical_transmission::Union{Nothing,Float64}=nothing,

        inoculumCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing, # takes seq argument, returns Float64
        mutationCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing, # takes seq argument, returns Float64
        recombinationCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing, # takes seq argument, returns Float64
        verticalTransmission::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing, # takes seq argument, returns Float64

        # Each element takes seq argument, returns Float64
        mutantEstablishmentCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        clearanceCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        responseAcquisitionCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        recombinantEstablishmentCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        contactCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        responseLossCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        birthCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        deathCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        transitionCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        receiveTransitionCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        receiveContactCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        intrahostFitnessCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing,
        )

    isnothing(num_loci) ? num_loci=template.num_loci : num_loci=num_loci
    isnothing(possible_alleles) ? possible_alleles=template.possible_alleles : possible_alleles=possible_alleles
    isnothing(mean_effective_inoculum) ? mean_effective_inoculum=template.mean_effective_inoculum : mean_effective_inoculum=mean_effective_inoculum
    isnothing(mean_mutations_per_replication) ? mean_mutations_per_replication=template.mean_mutations_per_replication : mean_mutations_per_replication=mean_mutations_per_replication
    isnothing(mean_recombination_crossovers) ? mean_recombination_crossovers=template.mean_recombination_crossovers : mean_recombination_crossovers=mean_recombination_crossovers
    isnothing(vertical_transmission) ? vertical_transmission=template.vertical_transmission : vertical_transmission=vertical_transmission
    isnothing(inoculumCoefficient) ? inoculumCoefficient=template.inoculumCoefficient : inoculumCoefficient=inoculumCoefficient
    isnothing(mutationCoefficient) ? mutationCoefficient=template.mutationCoefficient : mutationCoefficient=mutationCoefficient
    isnothing(recombinationCoefficient) ? recombinationCoefficient=template.recombinationCoefficient : recombinationCoefficient=recombinationCoefficient
    isnothing(verticalTransmission) ? verticalTransmission=template.verticalTransmission : verticalTransmission=verticalTransmission

    isnothing(mutantEstablishmentCoefficient) ? mutantEstablishmentCoefficient=template.coefficient_functions[MUTANT_ESTABLISHMENT] : mutantEstablishmentCoefficient=mutantEstablishmentCoefficient
    isnothing(clearanceCoefficient) ? clearanceCoefficient=template.coefficient_functions[CLEARANCE] : clearanceCoefficient=clearanceCoefficient
    isnothing(responseAcquisitionCoefficient) ? responseAcquisitionCoefficient=template.coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionCoefficient=responseAcquisitionCoefficient
    isnothing(recombinantEstablishmentCoefficient) ? recombinantEstablishmentCoefficient=template.coefficient_functions[RECOMBINANT_ESTABLISHMENT] : recombinantEstablishmentCoefficient=recombinantEstablishmentCoefficient
    isnothing(contactCoefficient) ? contactCoefficient=template.coefficient_functions[CONTACT] : contactCoefficient=contactCoefficient
    isnothing(responseLossCoefficient) ? responseLossCoefficient=template.coefficient_functions[RESPONSE_LOSS] : responseLossCoefficient=responseLossCoefficient
    isnothing(birthCoefficient) ? birthCoefficient=template.coefficient_functions[BIRTH] : birthCoefficient=birthCoefficient
    isnothing(deathCoefficient) ? deathCoefficient=template.coefficient_functions[DEATH] : deathCoefficient=deathCoefficient
    isnothing(transitionCoefficient) ? transitionCoefficient=template.coefficient_functions[TRANSITION] : transitionCoefficient=transitionCoefficient
    isnothing(receiveTransitionCoefficient) ? receiveTransitionCoefficient=template.coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionCoefficient=receiveTransitionCoefficient
    isnothing(receiveContactCoefficient) ? receiveContactCoefficient=template.coefficient_functions[RECEIVE_CONTACT] : receiveContactCoefficient=receiveContactCoefficient
    isnothing(intrahostFitnessCoefficient) ? intrahostFitnessCoefficient=template.coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessCoefficient=intrahostFitnessCoefficient

    return PathogenType(
        id,
        num_loci,
        possible_alleles,
        mean_effective_inoculum,
        mean_mutations_per_replication,
        mean_recombination_crossovers,
        vertical_transmission,
        inoculumCoefficient, mutationCoefficient, recombinationCoefficient, verticalTransmission,
        SA[ # order defined in COEFFICIENTS
            mutantEstablishmentCoefficient, clearanceCoefficient, responseAcquisitionCoefficient,
            recombinantEstablishmentCoefficient, contactCoefficient, responseLossCoefficient,
            birthCoefficient, deathCoefficient, transitionCoefficient,
            receiveTransitionCoefficient, receiveContactCoefficient, intrahostFitnessCoefficient
        ],
    )
end

function newResponseType(
        id::String;

        template::ResponseType=DEFAULT_RESPONSE_TYPE,

        inherit_response::Union{Nothing,Float64}=nothing,

        infectionCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        # takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient
        reactivityCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        # takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient

        # Each takes host, imprinted, matured sequences and returns Float64 coefficient
        mutantEstablishmentStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        clearanceStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        responseAcquisitionStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        recombinantEstablishmentStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        contactStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        responseLossStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        birthStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        deathStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        transitionStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        receiveTransitionStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        receiveContactStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,
        intrahostFitnessStaticCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String}}}=nothing,

        # Each takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient
        mutantEstablishmentSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        clearanceSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        responseAcquisitionSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        recombinantEstablishmentSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        contactSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        responseLossSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        birthSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        deathSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        transitionSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        receiveTransitionSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        receiveContactSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        intrahostFitnessSpecificCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String,String,String,String}}}=nothing,
        )

    isnothing(inherit_response) ? inherit_response=template.inherit_response : inherit_response=inherit_response
    isnothing(infectionCoefficient) ? infectionCoefficient=template.infectionCoefficient : infectionCoefficient=infectionCoefficient
    isnothing(reactivityCoefficient) ? reactivityCoefficient=template.reactivityCoefficient : reactivityCoefficient=reactivityCoefficient

    isnothing(mutantEstablishmentStaticCoefficient) ? mutantEstablishmentStaticCoefficient=template.static_coefficient_functions[MUTANT_ESTABLISHMENT] : mutantEstablishmentStaticCoefficient=mutantEstablishmentStaticCoefficient
    isnothing(clearanceStaticCoefficient) ? clearanceStaticCoefficient=template.static_coefficient_functions[CLEARANCE] : clearanceStaticCoefficient=clearanceStaticCoefficient
    isnothing(responseAcquisitionStaticCoefficient) ? responseAcquisitionStaticCoefficient=template.static_coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionStaticCoefficient=responseAcquisitionStaticCoefficient
    isnothing(recombinantEstablishmentStaticCoefficient) ? recombinantEstablishmentStaticCoefficient=template.static_coefficient_functions[RECOMBINANT_ESTABLISHMENT] : recombinantEstablishmentStaticCoefficient=recombinantEstablishmentStaticCoefficient
    isnothing(contactStaticCoefficient) ? contactStaticCoefficient=template.static_coefficient_functions[CONTACT] : contactStaticCoefficient=contactStaticCoefficient
    isnothing(responseLossStaticCoefficient) ? responseLossStaticCoefficient=template.static_coefficient_functions[RESPONSE_LOSS] : responseLossStaticCoefficient=responseLossStaticCoefficient
    isnothing(birthStaticCoefficient) ? birthStaticCoefficient=template.static_coefficient_functions[BIRTH] : birthStaticCoefficient=birthStaticCoefficient
    isnothing(deathStaticCoefficient) ? deathStaticCoefficient=template.static_coefficient_functions[DEATH] : deathStaticCoefficient=deathStaticCoefficient
    isnothing(transitionStaticCoefficient) ? transitionStaticCoefficient=template.static_coefficient_functions[TRANSITION] : transitionStaticCoefficient=transitionStaticCoefficient
    isnothing(receiveTransitionStaticCoefficient) ? receiveTransitionStaticCoefficient=template.static_coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionStaticCoefficient=receiveTransitionStaticCoefficient
    isnothing(receiveContactStaticCoefficient) ? receiveContactStaticCoefficient=template.static_coefficient_functions[RECEIVE_CONTACT] : receiveContactStaticCoefficient=receiveContactStaticCoefficient
    isnothing(intrahostFitnessStaticCoefficient) ? intrahostFitnessStaticCoefficient=template.static_coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessStaticCoefficient=intrahostFitnessStaticCoefficient
    isnothing(mutantEstablishmentSpecificCoefficient) ? mutantEstablishmentSpecificCoefficient=template.specific_coefficient_functions[MUTANT_ESTABLISHMENT] : mutantEstablishmentSpecificCoefficient=mutantEstablishmentSpecificCoefficient
    isnothing(clearanceSpecificCoefficient) ? clearanceSpecificCoefficient=template.specific_coefficient_functions[CLEARANCE] : clearanceSpecificCoefficient=clearanceSpecificCoefficient
    isnothing(responseAcquisitionSpecificCoefficient) ? responseAcquisitionSpecificCoefficient=template.specific_coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionSpecificCoefficient=responseAcquisitionSpecificCoefficient
    isnothing(recombinantEstablishmentSpecificCoefficient) ? recombinantEstablishmentSpecificCoefficient=template.specific_coefficient_functions[RECOMBINANT_ESTABLISHMENT] : recombinantEstablishmentSpecificCoefficient=recombinantEstablishmentSpecificCoefficient
    isnothing(contactSpecificCoefficient) ? contactSpecificCoefficient=template.specific_coefficient_functions[CONTACT] : contactSpecificCoefficient=contactSpecificCoefficient
    isnothing(responseLossSpecificCoefficient) ? responseLossSpecificCoefficient=template.specific_coefficient_functions[RESPONSE_LOSS] : responseLossSpecificCoefficient=responseLossSpecificCoefficient
    isnothing(birthSpecificCoefficient) ? birthSpecificCoefficient=template.specific_coefficient_functions[BIRTH] : birthSpecificCoefficient=birthSpecificCoefficient
    isnothing(deathSpecificCoefficient) ? deathSpecificCoefficient=template.specific_coefficient_functions[DEATH] : deathSpecificCoefficient=deathSpecificCoefficient
    isnothing(transitionSpecificCoefficient) ? transitionSpecificCoefficient=template.specific_coefficient_functions[TRANSITION] : transitionSpecificCoefficient=transitionSpecificCoefficient
    isnothing(receiveTransitionSpecificCoefficient) ? receiveTransitionSpecificCoefficient=template.specific_coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionSpecificCoefficient=receiveTransitionSpecificCoefficient
    isnothing(receiveContactSpecificCoefficient) ? receiveContactSpecificCoefficient=template.specific_coefficient_functions[RECEIVE_CONTACT] : receiveContactSpecificCoefficient=receiveContactSpecificCoefficient
    isnothing(intrahostFitnessSpecificCoefficient) ? intrahostFitnessSpecificCoefficient=template.specific_coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessSpecificCoefficient=intrahostFitnessSpecificCoefficient

    return ResponseType(
        id,
        inherit_response,
        infectionCoefficient,
        reactivityCoefficient,
        SA[ # order defined in COEFFICIENTS
            mutantEstablishmentStaticCoefficient, clearanceStaticCoefficient, responseAcquisitionStaticCoefficient,
            recombinantEstablishmentStaticCoefficient, contactStaticCoefficient, responseLossStaticCoefficient,
            birthStaticCoefficient, deathStaticCoefficient, transitionStaticCoefficient,
            receiveTransitionStaticCoefficient, receiveContactStaticCoefficient, intrahostFitnessStaticCoefficient
        ],
        SA[ # order defined in COEFFICIENTS
            mutantEstablishmentSpecificCoefficient, clearanceSpecificCoefficient, responseAcquisitionSpecificCoefficient,
            recombinantEstablishmentSpecificCoefficient, contactSpecificCoefficient, responseLossSpecificCoefficient,
            birthSpecificCoefficient, deathSpecificCoefficient, transitionSpecificCoefficient,
            receiveTransitionSpecificCoefficient, receiveContactSpecificCoefficient, intrahostFitnessSpecificCoefficient
        ],
    )
end

function newPopulationType(
        id::String;

        template::PopulationType=DEFAULT_POPULATION_TYPE,

        constant_contact_density::Union{Nothing,Bool}=nothing,
        constant_transition_density::Union{Nothing,Bool}=nothing,

        inoculum_coefficient::Union{Nothing,Float64}=nothing,
        mutation_coefficient::Union{Nothing,Float64}=nothing,
        recombination_coefficient::Union{Nothing,Float64}=nothing,

        host_num_loci::Union{Nothing,Int64}=nothing,
        host_possible_alleles::Union{Nothing,String}=nothing,
        host_mean_mutations_per_replication::Union{Nothing,Float64}=nothing,
        host_sexual_reproduction::Union{Nothing,Bool}=nothing,
        host_mean_recombination_crossovers::Union{Nothing,Float64}=nothing,

        hostSexualCompatibility::Union{Nothing,FunctionWrapper{Bool,Tuple{String,String}}}=nothing,
        hostMutationCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing, # takes seq argument, returns Float64
        hostRecombinationCoefficient::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing, # takes seq argument, returns Float64

        # Rate coefficients:
        mutant_establishment_coefficient::Union{Nothing,Float64}=nothing,
        clearance_coefficient::Union{Nothing,Float64}=nothing,
        response_acquisition_coefficient::Union{Nothing,Float64}=nothing,
        recombinant_establishment_coefficient::Union{Nothing,Float64}=nothing,
        contact_coefficient::Union{Nothing,Float64}=nothing,
        response_loss_coefficient::Union{Nothing,Float64}=nothing,
        birth_coefficient::Union{Nothing,Float64}=nothing,
        death_coefficient::Union{Nothing,Float64}=nothing,
        transition_coefficient::Union{Nothing,Float64}=nothing,

        # Receive coefficients
        receive_transition_coefficient::Union{Nothing,Float64}=nothing,
        receive_contact_coefficient::Union{Nothing,Float64}=nothing,
        intrahost_fitness_coefficient::Union{Nothing,Float64}=nothing,

        pathogenFractions::Union{Nothing, FunctionWrapper{Vector{Float64},Tuple{
                    Host,
                    FunctionWrapper{Float64,Tuple{Pathogen,Host,Int64}}
                    }}}=nothing,
        # Takes Host entity and Population's weightedResponse function,
        # returns vector with fractional representation of each pathogen present
        weightedResponse::Union{Nothing,FunctionWrapper{Float64,Tuple{Pathogen,Host,Int64}}}=nothing,
        # Takes Pathogen entity, Host entity, and event number;
        # returns aggregated response coefficient against that Pathogen for that event
        infectionProbability::Union{Nothing,FunctionWrapper{Float64,Tuple{Pathogen,Host}}}=nothing,
        # Takes Pathogen and Host entities,
        # returns probability that a contact results in successful infection given the Responses in Host

        developResponses::Union{Nothing,FunctionWrapper{Vector{Response},Tuple{Pathogen,Host,Vector{Response}}}}=nothing,
        # takes in Pathogen, Host, list of Responses as arguments, returns Response entities to be added
        # (this handles how many and which responses to choose when adding a response to a host)
        )

    isnothing(constant_contact_density) ? constant_contact_density=template.constant_contact_density : constant_contact_density=constant_contact_density
    isnothing(constant_transition_density) ? constant_transition_density=template.constant_transition_density : constant_transition_density=constant_transition_density
    isnothing(inoculum_coefficient) ? inoculum_coefficient=template.inoculum_coefficient : inoculum_coefficient=inoculum_coefficient
    isnothing(mutation_coefficient) ? mutation_coefficient=template.mutation_coefficient : mutation_coefficient=mutation_coefficient
    isnothing(recombination_coefficient) ? recombination_coefficient=template.recombination_coefficient : recombination_coefficient=recombination_coefficient

    isnothing(host_num_loci) ? host_num_loci=template.host_num_loci : host_num_loci=host_num_loci
    isnothing(host_possible_alleles) ? host_possible_alleles=template.host_possible_alleles : host_possible_alleles=host_possible_alleles
    isnothing(host_mean_mutations_per_replication) ? host_mean_mutations_per_replication=template.host_mean_mutations_per_replication : host_mean_mutations_per_replication=host_mean_mutations_per_replication
    isnothing(host_sexual_reproduction) ? host_sexual_reproduction=template.host_sexual_reproduction : host_sexual_reproduction=host_sexual_reproduction
    isnothing(host_mean_recombination_crossovers) ? host_mean_recombination_crossovers=template.host_mean_recombination_crossovers : host_mean_recombination_crossovers=host_mean_recombination_crossovers

    isnothing(hostSexualCompatibility) ? hostSexualCompatibility=template.hostSexualCompatibility : hostSexualCompatibility=hostSexualCompatibility
    isnothing(hostMutationCoefficient) ? hostMutationCoefficient=template.hostMutationCoefficient : hostMutationCoefficient=hostMutationCoefficient
    isnothing(hostRecombinationCoefficient) ? hostRecombinationCoefficient=template.hostRecombinationCoefficient : hostRecombinationCoefficient=hostRecombinationCoefficient

    isnothing(mutant_establishment_coefficient) ? mutant_establishment_coefficient=template.base_coefficients[MUTANT_ESTABLISHMENT] : mutant_establishment_coefficient=mutant_establishment_coefficient
    isnothing(clearance_coefficient) ? clearance_coefficient=template.base_coefficients[CLEARANCE] : clearance_coefficient=clearance_coefficient
    isnothing(response_acquisition_coefficient) ? response_acquisition_coefficient=template.base_coefficients[RESPONSE_ACQUISITION] : response_acquisition_coefficient=response_acquisition_coefficient
    isnothing(recombinant_establishment_coefficient) ? recombinant_establishment_coefficient=template.base_coefficients[RECOMBINANT_ESTABLISHMENT] : recombinant_establishment_coefficient=recombinant_establishment_coefficient
    isnothing(contact_coefficient) ? contact_coefficient=template.base_coefficients[CONTACT] : contact_coefficient=contact_coefficient
    isnothing(response_loss_coefficient) ? response_loss_coefficient=template.base_coefficients[RESPONSE_LOSS] : response_loss_coefficient=response_loss_coefficient
    isnothing(birth_coefficient) ? birth_coefficient=template.base_coefficients[BIRTH] : birth_coefficient=birth_coefficient
    isnothing(death_coefficient) ? death_coefficient=template.base_coefficients[DEATH] : death_coefficient=death_coefficient
    isnothing(transition_coefficient) ? transition_coefficient=template.base_coefficients[TRANSITION] : transition_coefficient=transition_coefficient
    isnothing(receive_transition_coefficient) ? receive_transition_coefficient=template.base_coefficients[RECEIVE_TRANSITION] : receive_transition_coefficient=receive_transition_coefficient
    isnothing(receive_contact_coefficient) ? receive_contact_coefficient=template.base_coefficients[RECEIVE_CONTACT] : receive_contact_coefficient=receive_contact_coefficient
    isnothing(intrahost_fitness_coefficient) ? intrahost_fitness_coefficient=template.base_coefficients[INTRAHOST_FITNESS] : intrahost_fitness_coefficient=intrahost_fitness_coefficient

    isnothing(pathogenFractions) ? pathogenFractions=template.pathogenFractions : pathogenFractions=pathogenFractions
    isnothing(weightedResponse) ? weightedResponse=template.weightedResponse : weightedResponse=weightedResponse
    isnothing(infectionProbability) ? infectionProbability=template.infectionProbability : infectionProbability=infectionProbability
    isnothing(developResponses) ? developResponses=template.developResponses : developResponses=developResponses

    return PopulationType(
        id,
        constant_contact_density, constant_transition_density,
        inoculum_coefficient,
        mutation_coefficient,
        recombination_coefficient,
        host_num_loci,
        host_possible_alleles,
        host_mean_mutations_per_replication,
        host_sexual_reproduction,
        host_mean_recombination_crossovers,
        hostSexualCompatibility,
        hostMutationCoefficient, # takes seq argument, returns Float64
        hostRecombinationCoefficient, # takes seq argument, returns Float64
        SA[ # order defined in COEFFICIENTS
            mutant_establishment_coefficient, clearance_coefficient, response_acquisition_coefficient,
            recombinant_establishment_coefficient, contact_coefficient, response_loss_coefficient,
            birth_coefficient, death_coefficient, transition_coefficient,
            receive_transition_coefficient, receive_contact_coefficient, intrahost_fitness_coefficient
        ],
        pathogenFractions,
        weightedResponse,
        infectionProbability,
        developResponses,
    )
end

# Model entity initializers

function newPathogen!(
    sequence::String, population::Population, type::PathogenType;
    parents::MVector{2,Union{Pathogen,Nothing}}=MVector{2,Union{Pathogen,Nothing}}([nothing, nothing]))
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
    imprinted_pathogen::Union{Pathogen,Nothing}, matured_pathogen::Union{Pathogen,Nothing}, host_sequence::String,
    population::Population, type::ResponseType;
    parents::MVector{2,Union{Response,Nothing}}=MVector{2,Union{Response,Nothing}}([nothing, nothing]))
    population.responses[(host_sequence, imprinted_pathogen.sequence, matured_pathogen.sequence, type.id)] = Response(
        parents, host_sequence, imprinted_pathogen, matured_pathogen,
        responseStaticCoefficients(host_sequence, imprinted_pathogen.sequence, matured_pathogen.sequence, type),
        type
    )

    return population.responses[(host_sequence, imprinted_pathogen.sequence, matured_pathogen.sequence, type.id)]
end

function newHost!(sequence::String, population::Population, model::Model)
    addHostToPopulation!(
        Host(
            length(population.hosts) + 1,
            sequence,
            population.parameters.host_mean_mutations_per_replication * population.type.hostMutationCoefficient(
                sequence
            ),
            population.parameters.host_mean_recombination_crossovers * population.type.hostRecombinationCoefficient(
                sequence
            ),
            Vector{Pathogen}(undef, 0), Vector{Response}(undef, 0),
            Vector{Float64}(undef, 0),
            Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, 0),
            Matrix{Float64}(undef, NUM_RESPONSE_EVENTS, 0),
        ),
        population, model
    )

    return population.hosts[end]
end

function staticHost(host::Host)
    return StaticHost(host.id, host.sequence, copy(host.pathogens), copy(host.responses))
end

function newPopulation!(id::String, parameters::PopulationType, model::Model)
    push!(model.populations, Population(
        id, parameters,
        Dict{String,Pathogen}(), Dict{Tuple{String,String,String,String},Response}(),
        Vector{Host}(undef, 0),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS - 1, 0),
        0.0, 0.0,
        zeros(Float64, length(model.populations)),
        zeros(Float64, length(model.populations)),
        zeros(Int64, NUM_COMPARTMENTS),
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
        zeros(SVector{NUM_EVENTS,Float64}),
        0.0
    )
end
