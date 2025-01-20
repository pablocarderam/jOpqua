using StaticArrays
using Random

# Basic actions

function mutantPathogen!(pathogen::Pathogen, host::Host, class::Class)
    locus = rand(1:pathogen.type.num_loci)
    seq = pathogen.sequence[1:locus-1] * rand(pathogen.type.possible_alleles) * pathogen.sequence[locus+1:end]

    return newPathogen!(seq, class)
end

function addPathogenToHost!(pathogen::Pathogen, host_idx::Int64, class::Class, population::Population, model::Model)
    push!(class.hosts[host_idx].pathogens, pathogen)
    push!(class.hosts[host_idx].pathogen_fractions, 0.0)
    class.hosts[host_idx].pathogen_weights = catCol(
        class.hosts[host_idx].pathogen_weights, zeros(Float64, NUM_PATHOGEN_EVENTS)
    )

    hostWeights!(host_idx, class, population, model)
end

function addResponseToHost!(response::Response, host_idx::Int64, class::Class, population::Population, model::Model)
    push!(class.hosts[host_idx].responses, response)
    class.hosts[host_idx].response_weights = catCol(
        class.hosts[host_idx].response_weights, zeros(Float64, NUM_RESPONSE_EVENTS)
    )

    hostWeights!(host_idx, class, population, model)
end

function removePathogenFromHost!(pathogen_idx::Int64, host_idx::Int64, class::Class, population::Population, model::Model)
    deleteat!(class.hosts[host_idx].pathogens, pathogen_idx)
    deleteat!(class.hosts[host_idx].pathogen_fractions, pathogen_idx)
    class.hosts[host_idx].pathogen_weights = class.hosts[host_idx].pathogen_weights[:, begin:end.!=pathogen_idx]

    hostWeights!(host_idx, class, population, model)
end

function removeResponseFromHost!(response_idx::Int64, host_idx::Int64, class::Class, population::Population, model::Model)
    deleteat!(class.hosts[host_idx].responses, response_idx)
    class.hosts[host_idx].response_weights = class.hosts[host_idx].response_weights[:, begin:end.!=response_idx]

    hostWeights!(host_idx, class, population, model)
end

# Model events

# function establishMutant!(pathogen::Pathogen, host::Host, class::Class, population::Population, model::Model)
#     mut = mutantPathogen!(pathogen, host, class)
#     addPathogenToHost!(mut, host, class, population, model)
# end
