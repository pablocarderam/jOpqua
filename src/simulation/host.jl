struct Host
    id::Int64
    pathogens::Vector{MAX_PATHOGENS,Pathogen}
    pathogen_rates::Vector{MAX_PATHOGENS,Float64}
    immunities::Vector{MAX_IMMUNITIES,Immunity}
    immunity_rates::Vector{MAX_IMMUNITIES,Float64}
    host_rates::Vector{NUM_HOST_EVENTS,Float64}
end
