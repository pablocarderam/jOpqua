using StaticArrays

# Basic actions

# function mutantPathogen(pathogen::Pathogen, host::Host, class::Class)

#     mut = newPathogen!(sequence, class)

#     return mut
# end

function addPathogenToHost!(pathogen::Pathogen, host::Host, class::Class, population::Population, model::Model)
    push!(host.pathogens, pathogen)
    push!(host.pathogen_fractions, 0.0)
    host.pathogen_weights = catCol(host.pathogen_weights, zeros(Float64, NUM_PATHOGEN_EVENTS))

    hostWeights!(host.id, class, population, model)
end

function addResponseToHost!(response::Response, host::Host, class::Class, population::Population, model::Model)
    push!(host.responses, response)
    host.response_weights = catCol(host.response_weights, zeros(Float64, NUM_RESPONSE_EVENTS))

    hostWeights!(host.id, class, population, model)
end

function removePathogenFromHost!(pathogen_idx::Int64, host::Host, class::Class, population::Population, model::Model)
    deleteat!(host.pathogens, pathogen_idx)
    deleteat!(host.pathogen_fractions, pathogen_idx)
    host.pathogen_weights = host.pathogen_weights[:, begin:end.!=pathogen_idx]

    hostWeights!(host.id, class, population, model)
end

function removeResponseFromHost!(response_idx::Int64, host::Host, class::Class, population::Population, model::Model)
    deleteat!(host.responses, response_idx)
    host.response_weights = host.response_weights[:, begin:end.!=response_idx]

    hostWeights!(host.id, class, population, model)
end

# Model events
