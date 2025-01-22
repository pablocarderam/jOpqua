using StaticArrays
using StatsBase
using PoissonRandom
using Random

# Basic actions

function mutantPathogen!(pathogen::Pathogen, class::Class)
    loci = rand(
        1:pathogen.type.num_loci, zeroTruncatedPoisson(
            class.parameters.base_coefficients[MUTATION_PER_REPLICATION] *
            pathogen.coefficients[MUTATION_PER_REPLICATION]
        )
    )
    seq = ""
    for locus in loci
        seq = seq * pathogen.sequence[length(seq)+1:locus-1] * rand(pathogen.type.possible_alleles)
    end
    seq = seq * pathogen.sequence[length(seq)+1:end]

    if haskey(class.pathogens, seq)
        return class.pathogens[seq]
    else
        return newPathogen!(seq, class, pathogen.type)
    end
end

function recombinantPathogens!(pathogen_1::Pathogen, pathogen_2::Pathogen, class::Class)
    children = MVector(pathogen_1.sequence, pathogen_2.sequence)

    if pathogen_1.type.mean_recombination_crossovers > 0 && pathogen_2.type.mean_recombination_crossovers > 0
        num_evts = zeroTruncatedPoisson(mean(
            pathogen_1.type.mean_recombination_crossovers,
            pathogen_2.type.mean_recombination_crossovers
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

    if !haskey(class.pathogens, children[1])
        newPathogen!(children[1], class, pathogen_1.type)
    end
    return class.pathogens[children[1]]
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

function establishMutant!(model::Model, rand_n::Float64)
    pathogen_idx, host_idx, class_idx, pop_idx = choosePathogen(MUTANT_ESTABLISHMENT, model, rand_n)

    mut = mutantPathogen!(
        model.populations[pop_idx].classes[class_idx].hosts[host_idx].pathogens[pathogen_idx],
        model.populations[pop_idx].classes[class_idx]
    )

    if !(mut in model.populations[pop_idx].classes[class_idx].hosts[host_idx].pathogens)
        addPathogenToHost!(
            mut, host_idx,
            model.populations[pop_idx].classes[class_idx],
            model.populations[pop_idx], model
        )
    end
end

function clearPathogen!(model::Model, rand_n::Float64)
    pathogen_idx, host_idx, class_idx, pop_idx = choosePathogen(CLEARANCE, model, rand_n)

    removePathogenFromHost!(
        pathogen_idx, host_idx,
        model.populations[pop_idx].classes[class_idx],
        model.populations[pop_idx], model
    )
end

function acquireResponse!(model::Model, rand_n::Float64)
    pathogen_idx, host_idx, class_idx, pop_idx = choosePathogen(CLEARANCE, model, rand_n)

    responses = model.populations[pop_idx].classes[class_idx].parameters.developResponses(
        model.populations[pop_idx].classes[class_idx].hosts[host_idx].pathogens[pathogen_idx],
        model.populations[pop_idx].classes[class_idx].hosts[host_idx],
        model.populations[pop_idx].classes[class_idx],
    )

    for response in responses
        if !(response in model.populations[pop_idx].classes[class_idx].hosts[host_idx].responses)
            addResponseToHost!(
                response, host_idx,
                model.populations[pop_idx].classes[class_idx],
                model.populations[pop_idx], model
            )
        end
    end
end

function establishRecombinant!(model::Model, rand_n::Float64)
    pathogen_idx_1, host_idx, class_idx, pop_idx, rand_n = choosePathogen(RECOMBINANT_ESTABLISHMENT, model, rand_n)
    pathogen_idx_2 = choosePathogen(host_idx, class_idx, pop_idx, RECOMBINANT_ESTABLISHMENT, model, rand_n)

    if pathogen_idx_1 != pathogen_idx_2 && pathogen_idx_1.type == pathogen_idx_2.type
        recombinant = recombinantPathogens!(
            model.populations[pop_idx].classes[class_idx].hosts[host_idx].pathogens[pathogen_idx_1],
            model.populations[pop_idx].classes[class_idx].hosts[host_idx].pathogens[pathogen_idx_2],
            model.populations[pop_idx].classes[class_idx]
        )

        if !(recombinant in model.populations[pop_idx].classes[class_idx].hosts[host_idx].pathogens)
            addPathogenToHost!(
                recombinant, host_idx,
                model.populations[pop_idx].classes[class_idx],
                model.populations[pop_idx], model
            )
        end
    end
end

function intraPopulationContact!(model::Model, rand_n::Float64)
    host_idx_1, class_idx_1, pop_idx, rand_n = chooseHost(INTRA_POPULATION_CONTACT, model, rand_n)
    class_idx_2, rand_n = chooseClass(pop_idx, RECEIVE_CONTACT, model, rand_n)
    host_idx_2, rand_n = chooseHost(class_idx_2, pop_idx, RECEIVE_CONTACT, model, rand_n)

    if host_idx_1 != host_idx_2 || class_idx_1 != class_idx_2
        inocula = MVector{length(model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens)}([
            pois_rand(
                model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx].type.mean_effective_inoculum *
                model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogen_fractions[p_idx]
            )
            for p_idx in 1:length(model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens)
        ])
        if !any(inocula .!= 0)
            idx, rand_n = randChoose(
                rand_n,
                model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogen_fractions,
                1.0, regenerate_rand=true
            )
            inocula[idx] += 1
        end

        for p_idx in eachindex(inocula)
            if inocula[p_idx] > 0
                mut_prob = min(
                    1.0,
                    model.populations[pop_idx].classes[class_idx_1].parameters.base_coefficients[MUTATION_PER_REPLICATION] *
                    model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx].coefficients[MUTATION_PER_REPLICATION]
                )
                rec_prob = 0.0
                if length(model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens) > 1
                    rec_prob = min(
                        1.0,
                        model.populations[pop_idx].classes[class_idx_1].parameters.base_coefficients[RECOMBINATION_PER_REPLICATION] *
                        model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx].coefficients[RECOMBINATION_PER_REPLICATION]
                    )
                end
                num_mut = binomial(inocula[p_idx], mut_prob * (1.0 - rec_prob))
                num_rec = binomial(inocula[p_idx], rec_prob * (1.0 - mut_prob))
                num_mut_rec = binomial(inocula[p_idx], mut_prob * rec_prob)
                if inocula[p_idx] - num_mut - num_rec - num_mut_rec > 0
                    if !(
                        model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx] in
                        model.populations[pop_idx].classes[class_idx_2].hosts[host_idx_2].pathogens)

                        if rand() < infectionProbability(
                            model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx],
                            model.populations[pop_idx].classes[class_idx_2].hosts[host_idx_2])

                            addPathogenToHost!(
                                model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx],
                                host_idx_2, model.populations[pop_idx].classes[class_idx_2],
                                model.populations[pop_idx], model
                            )
                        end
                    end
                end
                for _ in 1:num_mut
                    mut = mutantPathogen!(
                        model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx],
                        model.populations[pop_idx].classes[class_idx_1]
                    )

                    if (!(
                        mut in model.populations[pop_idx].classes[class_idx_2].hosts[host_idx_2].pathogens
                    )) && rand() < infectionProbability(
                        mut, model.populations[pop_idx].classes[class_idx_2].hosts[host_idx_2]
                    )

                        addPathogenToHost!(
                            mut, host_idx_2,
                            model.populations[pop_idx].classes[class_idx_2],
                            model.populations[pop_idx], model
                        )
                    end
                end
                for _ in 1:num_rec
                    p_idx_2 = choosePathogen(host_idx_1, class_idx_1, pop_idx, RECOMBINANT_ESTABLISHMENT, model, rand_n)

                    if p_idx != p_idx_2 && (
                        model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx].type ==
                        model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx_2].type)

                        recombinant = recombinantPathogens!(
                            model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx],
                            model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx_2],
                            model.populations[pop_idx].classes[class_idx_1]
                        )

                        if (!(
                            recombinant in model.populations[pop_idx].classes[class_idx_2].hosts[host_idx_2].pathogens
                        )) && rand() < infectionProbability(
                            recombinant, model.populations[pop_idx].classes[class_idx_2].hosts[host_idx_2]
                        )
                            addPathogenToHost!(
                                recombinant, host_idx_2,
                                model.populations[pop_idx].classes[class_idx_2],
                                model.populations[pop_idx], model
                            )
                        end
                    end
                end
                for _ in 1:num_mut_rec
                    p_idx_2 = choosePathogen(host_idx_1, class_idx_1, pop_idx, RECOMBINANT_ESTABLISHMENT, model, rand_n)

                    if p_idx != p_idx_2
                        recombinant = recombinantPathogens!(
                            model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx],
                            model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx_2],
                            model.populations[pop_idx].classes[class_idx_1]
                        )

                        recombinant_mutant = mutantPathogen!(recombinant, model.populations[pop_idx].classes[class_idx_1])

                        if (!(
                            recombinant_mutant in model.populations[pop_idx].classes[class_idx_2].hosts[host_idx_2].pathogens
                        )) && rand() < infectionProbability( # we use > because < implies immunity
                            recombinant_mutant, model.populations[pop_idx].classes[class_idx_2].hosts[host_idx_2]
                        )
                            addPathogenToHost!(
                                recombinant_mutant, host_idx_2,
                                model.populations[pop_idx].classes[class_idx_2],
                                model.populations[pop_idx], model
                            )
                        end
                    else
                        mut = mutantPathogen!(
                            model.populations[pop_idx].classes[class_idx_1].hosts[host_idx_1].pathogens[p_idx],
                            model.populations[pop_idx].classes[class_idx_1]
                        )

                        if !(mut in model.populations[pop_idx].classes[class_idx_2].hosts[host_idx_2].pathogens)
                            addPathogenToHost!(
                                mut, host_idx_2,
                                model.populations[pop_idx].classes[class_idx_2],
                                model.populations[pop_idx], model
                            )
                        end
                    end
                end
            end
        end
    end
end

event_functions = SA[
    establishMutant!, clearPathogen!, acquireResponse!,
    establishRecombinant!, intraPopulationContact!
]
