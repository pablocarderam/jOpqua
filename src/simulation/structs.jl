using StaticArrays

# Parameters structs (entity types):
struct PathogenType
    id::String

    num_loci::Int64
    possible_alleles::String

    mean_effective_inoculum::Float64
    mean_mutations_per_replication::Float64
    mean_recombination_crossovers::Int64

    vertical_transmission::Float64

    coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # Each element takes seq argument, returns Float64
    inoculumCoefficient::Function # takes seq argument, returns Float64
    mutationCoefficient::Function # takes seq argument, returns Float64
    recombinationCoefficient::Function # takes seq argument, returns Float64

    verticalTransmission::Function # takes seq argument, returns Float64
end

struct ResponseType
    id::String
    static_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # each takes imprinted, matured sequences and returns Float64 coefficient
    specific_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # each takes imprinted, matured, and infecting sequences and returns Float64 coefficient
    inherit_response::Float64
    infectionCoefficient::Function # takes imprinted, matured, and infecting sequences and returns Float64 coefficient
    reactivityCoefficient::Function # takes imprinted, matured, and infecting sequences and returns Float64 coefficient

end

struct PopulationType
    id::String

    base_coefficients::SVector{NUM_COEFFICIENTS,Float64}

    constant_contact_density::Bool
    constant_transition_density::Bool

    pathogenFractions::Function
    # Takes Host and Population entities, returns vector with fractional representation of each pathogen present
    weightedResponse::Function
    # Takes Pathogen entity, Host entity, and event number;
    # returns aggregated response coefficient against that Pathogen for that event
    infectionProbability::Function
    # Takes Pathogen and Host entities,
    # returns probability that a contact results in successful infection given the Responses in Host

    developResponses::Function
    # takes in Pathogen, Host, Population as arguments, returns Response objects to be added
    # (this handles how many and which responses to choose when adding a response to a host)
    #TODO: maybe doesn't need Population? Probably does to ensure responses don't already exist

    inoculum_coefficient::Float64
    mutation_coefficient::Float64
    recombination_coefficient::Float64
end

# Model entities:
struct Pathogen
    parents::MVector{2, Union{Pathogen, Nothing}} # parent pathogen objects, if any
    # This is only useful for response lineage tracing, but not the simulation?
    sequence::String
    coefficients::SVector{NUM_COEFFICIENTS,Float64}
    mean_effective_inoculum::Float64
    mean_mutations_per_replication::Float64
    mean_recombination_crossovers::Float64
    type::PathogenType
end

struct Response
    parents::MVector{2, Union{Response, Nothing}} # parent response objects, if any
    # This is only useful for response lineage tracing, but not the simulation?
    imprinted_pathogen::Pathogen # This will track the Pathogen imprinted in the naive response
    matured_pathogen::Pathogen
    coefficients::SVector{NUM_COEFFICIENTS,Float64} # static coefficients
    type::ResponseType
end

mutable struct Host
    id::Int64

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

mutable struct Population
    id::String
    parameters::PopulationType

    pathogens::Dict{String,Pathogen}
    responses::Dict{Tuple{String,String,String},Response}
    # keys are tuples of imprinted genome, matured genome, type ID
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
end

struct Intervention
    time::Float64
    intervention::Function
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

    interventions::Vector{Intervention}

    event_rates::MVector{NUM_EVENTS,Float64}
    event_rates_sum::Float64
end

mutable struct FlexLevel
    bounds::Tuple{Float64, Float64}
    sum::Float64
    max::Float64
    indices::Vector{Int64}
    next::Union{FlexLevel, Nothing}
end

mutable struct FlexlevSampler
    min_level::Float64
    max_level::Float64
    levels::FlexLevel
    weights::AbstractVector{Float64}
end