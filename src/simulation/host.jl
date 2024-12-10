struct Host
    id::Int64
    pathogens::Vector{MAX_PATHOGENS,Pathogen}
    pathogen_rates::Vector{MAX_PATHOGENS,Float64}
    immunities::Vector{MAX_IMMUNITIES,Immunity}
    immunity_rates::Vector{MAX_IMMUNITIES,Float64}
    event_rates::Vector{NUM_HOST_EVENTS+2,Float64} # plus pathogens and immunities
end
