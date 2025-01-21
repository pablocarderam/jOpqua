using StaticArrays

struct PathogenType
    id::String

    num_loci::Int64
    possible_alleles::String
    mean_recombination_crossovers::Int64
    mean_effective_inoculum::Float64

    coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # Each element takes seq argument, returns Float64
end

struct ResponseType
    id::String
    static_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # each takes imprinted, matured sequences and returns Float64 coefficient
    specific_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # each takes imprinted, matured, and infecting sequences and returns Float64 coefficient
end

struct ClassParameters
    id::String

    base_coefficients::SVector{NUM_COEFFICIENTS,Float64}

    class_change_fractions::Dict{String,Float64} # size CLASSES, must sum to 1

    inter_population_contact_fractions::Dict{String,Float64} # size POPULATIONS, must sum to 1
    migration_fractions::Dict{String,Float64} # size POPULATIONS, must sum to 1

    pathogen_types::Dict{String,PathogenType}
    response_types::Dict{String,ResponseType}
    developResponses::Function # takes in Pathogen, Host, Class as arguments, returns Response objects to be added
    # (this handles how many and which responses to choose when adding a response to a host)
    #TODO: maybe doesn't need Class? Probably does to ensure responses don't already exist
end

struct Pathogen
    sequence::String
    coefficients::SVector{NUM_COEFFICIENTS,Float64}
    type::PathogenType
end

struct Response
    parent::Tuple{String,String,String} # imprinted genome, matured genome, type ID
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
    # and larger at the Class level.
end

mutable struct Class
    id::String
    parameters::ClassParameters

    pathogens::Dict{String,Pathogen}
    responses::Dict{Tuple{String,String,String},Response}
    # keys are tuples of imprinted genome, matured genome, type ID
    hosts::Vector{Host} # size MAX_HOSTS

    host_weights::Matrix{Float64} # size NUM_EVENTS x MAX_HOSTS
    host_weights_receive::Matrix{Float64}
    # size NUM_CHOICE_MODIFIERS-1 x MAX_HOSTS; -1 excludes intrahost fitness
end

mutable struct Population
    id::String

    classes::Vector{Class} # size CLASSES
    class_dict::Dict{String,Int64}
    # Holds Class ids => Class indexes, used for migration to get matching class

    class_weights::Matrix{Float64} # size NUM_EVENTS x CLASSES
    class_weights_receive::Matrix{Float64}
    # size NUM_CHOICE_MODIFIERS-2 x CLASSES; -2 excludes intrahost fitness, host receive contact rates
end

mutable struct Model
    populations::Vector{Population} # size POPULATIONS
    population_dict::Dict{String,Int64}
    # Holds Population ids => Population indexes, used for migration to get matching Population

    population_weights::Matrix{Float64} # size NUM_EVENTS x POPULATIONS
    population_weights_receive::Matrix{Float64}
    # size NUM_CHOICE_MODIFIERS-3 x POPULATIONS; -3 excludes intrahost fitness, host receive contact rates, and class change

    population_weights_receive_sums::MVector{NUM_CHOICE_MODIFIERS-3,Float64}
    # size NUM_CHOICE_MODIFIERS-3; -3 excludes intrahost fitness, host receive contact rates, and class change

    event_rates::MVector{NUM_EVENTS,Float64}
    event_rates_sum::Float64
end
