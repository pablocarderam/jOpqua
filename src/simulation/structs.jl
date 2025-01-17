using StaticArrays

struct ImmunityType
    id::String
    static_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # each takes imprinted, matured sequences and returns Float64 coefficient
    specific_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # each takes imprinted, matured, and infecting sequences and returns Float64 coefficient
    immunodominance::Function # takes imprinted, matured, and infecting sequences and returns Float64 coefficient
end

struct ClassParameters
    id::String

    base_coefficients::SVector{NUM_COEFFICIENTS,Float64}

    class_change_fractions::Dict{String,Float64} # size CLASSES, must sum to 1

    inter_population_contact_fractions::Dict{String,Float64} # size POPULATIONS, must sum to 1
    migration_fractions::Dict{String,Float64} # size POPULATIONS, must sum to 1

    pathogen_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # Each element takes seq argument, returns Float64

    immunity_types::Dict{String,ImmunityType}
    acquireImmunities::Function # takes in Pathogen, Host, Class as arguments, returns Immunity objects to be added
    # (this handles how many and which immunities to choose when immunizing a host)

end

struct Pathogen
    id::Int64
    sequence::String
    coefficients::SVector{NUM_COEFFICIENTS,Float64}
end

struct Immunity
    id::Int64
    imprinted_pathogen::Pathogen
    matured_pathogen::Pathogen
    coefficients::SVector{NUM_COEFFICIENTS,Float64} # static coefficients
    type::ImmunityType
end

mutable struct Host
    id::Int64

    pathogens::Vector{Pathogen} # size MAX_PATHOGENS
    immunities::Vector{Immunity} # size MAX_IMMUNITIES

    pathogen_fractions::Vector{Float64} # size MAX_PATHOGENS

    pathogen_weights::Matrix{Float64} # size NUM_PATHOGEN_EVENTS x MAX_PATHOGENS
    immunity_weights::Matrix{Float64} # size NUM_IMMUNITY_EVENTS x MAX_IMMUNITIES

    # We could handle receiving rates at this level, but the only one relevant
    # to entities within hosts is recombination, and that uses the same rates
    # already calculated above since it's symmetrical. We therefore don't define
    # any additional matrices here and we handle all receiving rates for hosts
    # and larger at the Class level.
end

mutable struct Class
    id::String
    parameters::ClassParameters

    pathogens::Vector{Pathogen}
    immunities::Vector{Immunity}
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

# struct Model
#     id::Int64
#     populations::MVector{POPULATIONS,Population}
#     pathogen_rates::MVector{POPULATIONS,Float64}
#     immunity_rates::MVector{POPULATIONS,Float64}
#     host_rates::MVector{POPULATIONS,Float64}
# end

mutable struct Model
    populations::Vector{Population} # size POPULATIONS
    population_dict::Dict{String,Int64}
    # Holds Population ids => Population indexes, used for migration to get matching Population

    population_weights::Matrix{Float64} # size NUM_EVENTS x POPULATIONS
    population_weights_receive::Matrix{Float64}
    # size NUM_CHOICE_MODIFIERS-3 x POPULATIONS; -3 excludes intrahost fitness, host receive contact rates, and class change

    event_rates::MVector{NUM_EVENTS,Float64}
end
