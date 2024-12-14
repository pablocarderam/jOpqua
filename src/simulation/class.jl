mutable struct Class
    id::String
    parameters::ClassParameters
    pathogens::Vector{Pathogen}
    immunities::Vector{Immunity}
    hosts::Vector{Host} # size MAX_HOSTS

    pathogen_rates::Matrix{Float64} # size NUM_EVENTS x MAX_HOSTS
    immunity_rates::Matrix{Float64} # size NUM_EVENTS x MAX_HOSTS
    host_rates::Matrix{Float64} # size NUM_EVENTS x MAX_HOSTS
end

function newClass!(population::Population, parameters::ClassParameters, id::String)
    class = Class(
        id, parameters,
        Vector{Pathogen}(undef, 0), Vector{Immunity}(undef, 0), Vector{Host}(undef, 0),
        Vector{Float64}(undef, 0), Vector{Float64}(undef, 0), Vector{Float64}(undef, 0)
    )
    # pathogenRates!(class)
    # immunityRates!(class)
    # hostRates!(class)
    push!(population.classes, class)
end

#TODO: at this level we have pathogens and immunities change rates for host events;
# we also have actual rate calculation from dimensionless weight coefficients times class rates
