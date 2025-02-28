using StaticArrays
using FunctionWrappers
import FunctionWrappers: FunctionWrapper

# Lower level (Pathogens, Responses)
# Parameters structs (entity types):

struct PathogenType
    id::String

    num_loci::Int64
    possible_alleles::String

    mean_effective_inoculum::Float64
    mean_mutations_per_replication::Float64
    mean_recombination_crossovers::Float64

    inoculumCoefficient::FunctionWrapper{Float64,Tuple{String}} # takes seq argument, returns Float64
    mutationCoefficient::FunctionWrapper{Float64,Tuple{String}} # takes seq argument, returns Float64
    recombinationCoefficient::FunctionWrapper{Float64,Tuple{String}} # takes seq argument, returns Float64

    verticalTransmissionCoefficient::FunctionWrapper{Float64,Tuple{String}} # takes seq argument, returns Float64
    # acquireResponseUponClearanceCoefficient::FunctionWrapper{Float64,Tuple{String}} # takes seq argument, returns Float64

    coefficient_functions::SVector{NUM_COEFFICIENTS,FunctionWrapper{Float64,Tuple{String}}}
    # Each element takes seq argument, returns Float64
end

struct ResponseType
    id::String
    inherit_response::Float64
    inoculumCoefficient::FunctionWrapper{Float64,Tuple{String,String,String,String}}
    # takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient
    verticalTransmissionCoefficient::FunctionWrapper{Float64,Tuple{String,String,String,String}}
    # takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient
    infectionCoefficient::FunctionWrapper{Float64,Tuple{String,String,String,String}}
    # takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient
    reactivityCoefficient::FunctionWrapper{Float64,Tuple{String,String,String,String}}
    # takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient
    static_coefficient_functions::SVector{NUM_COEFFICIENTS,FunctionWrapper{Float64,Tuple{String,String,String}}}
    # each takes host, imprinted, matured sequences and returns Float64 coefficient
    specific_coefficient_functions::SVector{NUM_COEFFICIENTS,FunctionWrapper{Float64,Tuple{String,String,String,String}}}
    # each takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient
end

# Model entities:

struct Pathogen
    parents::MVector{2,Union{Pathogen,Nothing}} # parent pathogen objects, if any
    birth_time::Float64
    # This is only useful for response lineage tracing, but not the simulation?
    sequence::String
    coefficients::SVector{NUM_COEFFICIENTS,Float64}
    mean_effective_inoculum::Float64
    mean_mutations_per_replication::Float64
    mean_recombination_crossovers::Float64
    vertical_transmission_coefficient::Float64
    type::PathogenType
end

struct Response
    parents::MVector{2,Union{Response,Nothing}} # parent response objects, if any
    birth_time::Float64
    # This is only useful for response lineage tracing, but not the simulation?
    host_sequence::String
    imprinted_pathogen::Union{Pathogen,Nothing} # This will track the Pathogen imprinted in the naive response
    matured_pathogen::Union{Pathogen,Nothing}
    coefficients::SVector{NUM_COEFFICIENTS,Float64} # static coefficients
    type::ResponseType
end

struct StaticHost
    id::Int64
    sequence::String
    pathogens::Vector{Pathogen} # size MAX_PATHOGENS
    responses::Vector{Response} # size MAX_RESPONSES
end

mutable struct Host
    id::Int64

    parents::MVector{2,Union{Host,Nothing}} # parent response objects, if any
    birth_time::Float64

    sequence::String
    mean_mutations_per_replication::Float64
    mean_recombination_crossovers::Float64

    pathogens::Vector{Pathogen} # size MAX_PATHOGENS
    responses::Vector{Response} # size MAX_RESPONSES

    pathogen_fractions::Vector{Float64} # size MAX_PATHOGENS

    pathogen_weights::Matrix{Float64} # size NUM_PATHOGEN_EVENTS x MAX_PATHOGENS
    response_weights::Matrix{Float64} # size NUM_RESPONSE_EVENTS x MAX_RESPONSES

    # We could handle receiving rates at this level, but the only one relevant
    # to entities within hosts is recombination, and that uses the same rates
    # already calculated above since it's symmetrical. We therefore don't define
    # any additional matrices here and we handle all receiving rates for hosts
    # and larger at the Population level.
end

# Higher-level (Populations, Model)
# Parameters structs (entity types):

struct PopulationType
    id::String

    constant_contact_density::Bool
    constant_transition_density::Bool

    inoculum_coefficient::Float64
    mutation_coefficient::Float64
    recombination_coefficient::Float64
    vertical_transmission_coefficient::Float64

    host_num_loci::Int64
    host_possible_alleles::String
    host_mean_mutations_per_replication::Float64
    host_sexual_reproduction::Bool
    host_mean_recombination_crossovers::Float64

    hostSexualCompatibility::FunctionWrapper{Bool,Tuple{String,String}}
    hostMutationCoefficient::FunctionWrapper{Float64,Tuple{String}} # takes seq argument, returns Float64
    hostRecombinationCoefficient::FunctionWrapper{Float64,Tuple{String}} # takes seq argument, returns Float64

    base_coefficients::SVector{NUM_COEFFICIENTS - 1,Float64}

    pathogenFractions::FunctionWrapper{Vector{Float64},Tuple{Host,FunctionWrapper{Float64,Tuple{Pathogen,Host,Int64}}}}
    # Takes Host entity and Population's weightedResponse function,
    # returns vector with fractional representation of each pathogen present
    weightedResponse::FunctionWrapper{Float64,Tuple{Pathogen,Host,Int64}}
    # Takes Pathogen entity, Host entity, and event number;
    # returns aggregated response coefficient against that Pathogen for that event
    infectionProbability::FunctionWrapper{Float64,Tuple{Pathogen,Host}}
    # Takes Pathogen and Host entities,
    # returns probability that a contact results in successful infection given the Responses in Host

    developResponses::FunctionWrapper{Vector{Response},Tuple{Pathogen,Host,Dict{Tuple{String,String,String,String},Response},Dict{String,ResponseType},Float64}}
    # takes in Pathogen, Host, population's Dict of Responses, population type's dictionary of ResponseTypes, and birth time as arguments, returns Response entities to be added
    # (this handles how many and which responses to choose when adding a response to a host)
    # The list of Responses is the Population level vector of Responses,
    # and is used to return a reference to an existing Response struct rather than
    # creating a new instance of the struct

    response_types::Dict{String,ResponseType}
end

mutable struct Population
    id::String
    parameters::PopulationType

    pathogens::Dict{String,Pathogen}
    responses::Dict{Tuple{String,String,String,String},Response}
    # keys are tuples of host genome, imprinted genome, matured genome, type ID
    hosts::Vector{Host} # size MAX_HOSTS

    host_weights::Matrix{Float64} # size NUM_EVENTS x MAX_HOSTS
    host_weights_receive::Matrix{Float64}
    # size NUM_CHOICE_MODIFIERS-1 x MAX_HOSTS; -1 excludes intrahost fitness

    host_weights_with_coefficient::Matrix{Float64} # size NUM_EVENTS x MAX_HOSTS
    host_weights_receive_with_coefficient::Matrix{Float64}
    # size NUM_CHOICE_MODIFIERS-1 x MAX_HOSTS; -1 excludes intrahost fitness

    contact_sum::Float64
    transition_sum::Float64

    population_contact_coefficients::Vector{Float64} # size POPULATIONS
    population_transition_coefficients::Vector{Float64} # size POPULATIONS

    compartment_vars::MVector{NUM_COMPARTMENTS,Int64}
    # uninfected naive, infected naive, uninfected immune, infected immune, dead
end

mutable struct Model
    populations::Vector{Population} # size POPULATIONS
    population_dict::Dict{String,Int64}
    # Holds Population ids => Population indexes, used for transition to get matching Population

    population_weights::Matrix{Float64} # size NUM_EVENTS x POPULATIONS
    population_weights_receive::Matrix{Float64}
    # size NUM_CHOICE_MODIFIERS-1 x POPULATIONS; -1 excludes intrahost fitness

    population_weights_receive_sums::MVector{NUM_CHOICE_MODIFIERS - 1,Float64}
    # size NUM_CHOICE_MODIFIERS-1; -1 excludes intrahost fitness

    population_contact_weights_receive::Matrix{Float64}
    # size POPULATIONS x POPULATIONS; receiver is on row, emitter is on column
    population_transition_weights_receive::Matrix{Float64}
    # size POPULATIONS x POPULATIONS; receiver is on row, emitter is on column

    population_contact_weights_receive_sums::Vector{Float64}
    # size POPULATIONS
    population_transition_weights_receive_sums::Vector{Float64}
    # size POPULATIONS

    event_rates::MVector{NUM_EVENTS,Float64}
    event_rates_sum::Float64

    time::Float64
end

struct Intervention
    time::Float64
    intervention::FunctionWrapper{Nothing,Tuple{Model}} # Takes Model as argument
end

struct Output
    model::Model
    time::Vector{Float64}
    compartment_vars::Dict{String,Matrix{Int64}}
    host_samples::Dict{String,Matrix{StaticHost}}
    # Matrix is of size host sample size x time vector length
end
