using StaticArrays
using StatsBase
using PoissonRandom
using Random

# General actions

function mutantPathogen!(pathogen::Pathogen, population::Population)
    loci = rand(
        1:pathogen.type.num_loci, zeroTruncatedPoisson(pathogen.mean_mutations_per_replication)
    )
    seq = ""
    for locus in loci
        seq = seq * pathogen.sequence[length(seq)+1:locus-1] * rand(pathogen.type.possible_alleles)
    end
    seq = seq * pathogen.sequence[length(seq)+1:end]

    if haskey(population.pathogens, seq)
        return population.pathogens[seq]
    else
        return newPathogen!(seq, population, pathogen.type)
    end
end

function recombinantPathogens!(pathogen_1::Pathogen, pathogen_2::Pathogen, population::Population)
    children = MVector(pathogen_1.sequence, pathogen_2.sequence)

    if pathogen_1.mean_recombination_crossovers > 0 && pathogen_2.mean_recombination_crossovers > 0
        num_evts = zeroTruncatedPoisson(mean(
            pathogen_1.mean_recombination_crossovers,
            pathogen_2.mean_recombination_crossovers
        ))
        loci = rand(1:pathogen_1.type.num_loci, num_evts)
        for l in loci
            temp_child_1 = children[1]
            children[1] = children[1][1:l-1] + children[2][l:end]
            children[2] = children[1][1:l-1] + temp_child_1[l:end]
        end
    else
        num_evts = 0
    end

    children = MVector{2}([split(seq, CHROMOSOME_SEPARATOR) for seq in children])
    parent = rand(0:1, length(children[1]))

    children = SVector([
        join([children[parent[i]+1][i] for i in 1:length(children[1])], CHROMOSOME_SEPARATOR),
        join([children[(parent[i]!=true)+1][i] for i in 1:length(children[2])], CHROMOSOME_SEPARATOR)
    ])

    # for seq in children
    #     if !haskey(class.pathogens, seq)
    #         newPathogen!(seq, class, pathogen_1.type)
    #     end
    # end

    # return SA[class.pathogens[children[1]], class.pathogens[children[2]]]

    if !haskey(population.pathogens, children[1])
        newPathogen!(children[1], population, pathogen_1.type)
    end
    return population.pathogens[children[1]]
end

function addPathogenToHost!(pathogen::Pathogen, host_idx::Int64, population::Population, model::Model)
    push!(population.hosts[host_idx].pathogens, pathogen)
    push!(population.hosts[host_idx].pathogen_fractions, 0.0)
    population.hosts[host_idx].pathogen_weights = catCol(
        population.hosts[host_idx].pathogen_weights, zeros(Float64, NUM_PATHOGEN_EVENTS)
    )

    hostWeights!(host_idx, population, model)
end

function addResponseToHost!(response::Response, host_idx::Int64, population::Population, model::Model)
    push!(population.hosts[host_idx].responses, response)
    population.hosts[host_idx].response_weights = catCol(
        population.hosts[host_idx].response_weights, zeros(Float64, NUM_RESPONSE_EVENTS)
    )

    hostWeights!(host_idx, population, model)
end

function removePathogenFromHost!(pathogen_idx::Int64, host_idx::Int64, population::Population, model::Model)
    deleteat!(population.hosts[host_idx].pathogens, pathogen_idx)
    deleteat!(population.hosts[host_idx].pathogen_fractions, pathogen_idx)
    population.hosts[host_idx].pathogen_weights = population.hosts[host_idx].pathogen_weights[:, begin:end.!=pathogen_idx]

    hostWeights!(host_idx, population, model)
end

function removeResponseFromHost!(response_idx::Int64, host_idx::Int64, population::Population, model::Model)
    deleteat!(population.hosts[host_idx].responses, response_idx)
    population.hosts[host_idx].response_weights = population.hosts[host_idx].response_weights[:, begin:end.!=response_idx]

    hostWeights!(host_idx, population, model)
end

function attemptInfection!(pathogen::Pathogen,
        host_idx::Int64, pop_idx::Int64, model::Model)
    if !(
        pathogen in
        model.populations[pop_idx].hosts[host_idx].pathogens
        ) &&
        rand() < infectionProbability(
            pathogen,
            model.populations[pop_idx].hosts[host_idx]
        )

        addPathogenToHost!(
            pathogen, host_idx, model.populations[pop_idx], model
        )
    end
end


function hostContact!(model::Model, rand_n::Float64)
    host_idx_1, pop_idx_1, rand_n = chooseHost(CONTACT, model, rand_n)
    pop_idx_2 = pop_idx_1
    # if !same_population
    #     pop_idx_2, rand_n = choosePopulation(RECEIVE_CONTACT, model, rand_n)
    #     #TODO: this sampling needs to be done according to contact rates from pop1
    # end
    # class_idx_2, rand_n = chooseClass(pop_idx_2, RECEIVE_CONTACT, model, rand_n)
    host_idx_2, rand_n = chooseHost(pop_idx_2, RECEIVE_CONTACT, model, rand_n)

    if host_idx_1 != host_idx_2 || pop_idx_1 != pop_idx_2
        inocula = MVector{length(model.populations[pop_idx_1].hosts[host_idx_1].pathogens)}([
            pois_rand(
                model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx].mean_effective_inoculum *
                model.populations[pop_idx_1].hosts[host_idx_1].pathogen_fractions[p_idx]
            )
            for p_idx in 1:length(model.populations[pop_idx_1].hosts[host_idx_1].pathogens)
        ])
        if !any(inocula .!= 0)
            idx, rand_n = randChoose(
                rand_n,
                model.populations[pop_idx_1].hosts[host_idx_1].pathogen_fractions,
                1.0, regenerate_rand=true
            )
            inocula[idx] += 1
        end

        for p_idx in eachindex(inocula)
            if inocula[p_idx] > 0
                mut_prob = min(
                    1.0, 1.0-exp(
                        -model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx].mean_mutations_per_replication
                    )
                )
                # probability from Poisson PMF with k=0
                rec_prob = 0.0
                if length(model.populations[pop_idx_1].hosts[host_idx_1].pathogens) > 1
                    rec_prob = min(
                        1.0, 1.0-exp(
                            -model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx].mean_recombination_crossovers
                        )
                    )
                    # probability from Poisson PMF with k=0
                end
                num_mut = binomial(inocula[p_idx], mut_prob * (1.0 - rec_prob))
                num_rec = binomial(inocula[p_idx], rec_prob * (1.0 - mut_prob))
                num_mut_rec = binomial(inocula[p_idx], mut_prob * rec_prob)
                if inocula[p_idx] - num_mut - num_rec - num_mut_rec > 0
                    attemptInfection!(
                        model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx],
                        host_idx_2, pop_idx_2, model
                    )
                end
                for _ in 1:num_mut
                    attemptInfection!(
                        mutantPathogen!(
                            model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx],
                            model.populations[pop_idx_1]
                        ), host_idx_2, pop_idx_2, model
                    )
                end
                for _ in 1:num_rec
                    p_idx_2 = choosePathogen(host_idx_1, pop_idx_1, RECOMBINANT_ESTABLISHMENT, model, rand_n)

                    if p_idx != p_idx_2 && (
                        model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx].type ==
                        model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx_2].type)

                        attemptInfection!(
                            recombinantPathogens!(
                                model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx],
                                model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx_2],
                                model.populations[pop_idx_1]
                            ), host_idx_2, pop_idx_2, model
                        )
                    end
                end
                for _ in 1:num_mut_rec
                    p_idx_2 = choosePathogen(host_idx_1, pop_idx_1, RECOMBINANT_ESTABLISHMENT, model, rand_n)

                    if p_idx != p_idx_2
                        recombinant = recombinantPathogens!(
                            model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx],
                            model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx_2],
                            model.populations[pop_idx_1]
                        )

                        attemptInfection!(
                            mutantPathogen!(recombinant, model.populations[pop_idx_1]),
                            host_idx_2, pop_idx_2, model
                        )
                    else
                        attemptInfection!(
                            mutantPathogen!(
                                model.populations[pop_idx_1].hosts[host_idx_1].pathogens[p_idx],
                                model.populations[pop_idx_1]
                            ),
                            host_idx_2, pop_idx_2, model
                        )
                    end
                end
            end
        end
    end
end

# Model events

function establishMutant!(model::Model, rand_n::Float64)
    pathogen_idx, host_idx, pop_idx = choosePathogen(MUTANT_ESTABLISHMENT, model, rand_n)

    mut = mutantPathogen!(
        model.populations[pop_idx].hosts[host_idx].pathogens[pathogen_idx],
        model.populations[pop_idx]
    )

    attemptInfection!(mut, host_idx, pop_idx, model)
end

function clearPathogen!(model::Model, rand_n::Float64)
    pathogen_idx, host_idx, pop_idx = choosePathogen(CLEARANCE, model, rand_n)

    removePathogenFromHost!(
        pathogen_idx, host_idx,
        model.populations[pop_idx], model
    )
end

function acquireResponse!(model::Model, rand_n::Float64)
    pathogen_idx, host_idx, pop_idx = choosePathogen(CLEARANCE, model, rand_n)

    responses = model.populations[pop_idx].parameters.developResponses(
        model.populations[pop_idx].hosts[host_idx].pathogens[pathogen_idx],
        model.populations[pop_idx].hosts[host_idx],
        model.populations[pop_idx],
    )

    for response in responses
        if !(response in model.populations[pop_idx].hosts[host_idx].responses)
            addResponseToHost!(
                response, host_idx,
                model.populations[pop_idx], model
            )
        end
    end
end

function establishRecombinant!(model::Model, rand_n::Float64)
    pathogen_idx_1, host_idx, pop_idx, rand_n = choosePathogen(RECOMBINANT_ESTABLISHMENT, model, rand_n)
    pathogen_idx_2 = choosePathogen(host_idx, pop_idx, RECOMBINANT_ESTABLISHMENT, model, rand_n)

    if pathogen_idx_1 != pathogen_idx_2 && pathogen_idx_1.type == pathogen_idx_2.type
        recombinant = recombinantPathogens!(
            model.populations[pop_idx].hosts[host_idx].pathogens[pathogen_idx_1],
            model.populations[pop_idx].hosts[host_idx].pathogens[pathogen_idx_2],
            model.populations[pop_idx]
        )

        attemptInfection!(recombinant, host_idx, pop_idx, model)
    end
end

event_functions = SA[
    establishMutant!, clearPathogen!, acquireResponse!,
    establishRecombinant!, hostContact!,
]
