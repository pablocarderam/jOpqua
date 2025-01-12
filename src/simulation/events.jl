using StaticArrays

function infect!(pathogen::Pathogen, host::Host, class::Class, population::Population, model::Model)
    push!(host.pathogens, pathogen)
    push!(host.pathogen_fractions, 0.0)
    host.pathogen_weights = catCol(host.pathogen_weights, zeros(Float64,NUM_PATHOGEN_EVENTS))

    hostWeights!(host.id, class, population, model)
end

function immunize!(immunity::Immunity, host::Host, class::Class, population::Population, model::Model)
    push!(host.immunities, immunity)
    host.immunity_weights = catCol(host.immunity_weights, zeros(Float64,NUM_IMMUNITY_EVENTS))

    hostWeights!(host.id, class, population, model)
end

function clear!(pathogen_idx::Int64, host::Host, class::Class, population::Population, model::Model)
    deleteat!(host.pathogens, pathogen_idx)
    deleteat!(host.pathogen_fractions, pathogen_idx)
    host.pathogen_weights = host.pathogen_weights[:, begin:end .!= pathogen_idx]

    hostWeights!(host.id, class, population, model)
end

function deimmunize!(immunity_idx::Int64, host::Host, class::Class, population::Population, model::Model)
    deleteat!(host.immunities, immunity_idx)
    host.immunity_weights = host.immunity_weights[:, begin:end .!= immunity_idx]

    hostWeights!(host.id, class, population, model)
end
