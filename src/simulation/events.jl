using StaticArrays
using StatsBase
using PoissonRandom
using Random

# General actions

function mutantSequence!(
    sequence::String, num_loci::Int64, possible_alleles::String,
    mean_mutations_per_replication::Float64)

    chromosomes = split(sequence, CHROMOSOME_SEPARATOR)
    for (i, chromosome_pair) in enumerate(chromosomes)
        homologous_chromosomes = split(chromosome_pair, HOMOLOGOUS_CHROMOSOME_SEPARATOR)
        for (j, chromosome) in enumerate(homologous_chromosomes)
            loci = rand(
                1:length(chromosome), zeroTruncatedPoisson(mean_mutations_per_replication * length(chromosome) / num_loci)
            )
            seq = ""
            for locus in loci
                seq = seq * chromosome[length(seq)+1:locus-1] * rand(possible_alleles)
            end
            homologous_chromosomes[j] = seq * sequence[length(seq)+1:end]
        end
        chromosomes[i] = join(homologous_chromosomes, HOMOLOGOUS_CHROMOSOME_SEPARATOR)
    end

    return join(chromosomes, CHROMOSOME_SEPARATOR)
end

function mutantPathogen!(pathogen::Pathogen, population::Population, birth_time::Float64)
    seq = mutantSequence!(
        pathogen.sequence, pathogen.type.num_loci, pathogen.type.possible_alleles,
        pathogen.mean_mutations_per_replication
    )

    if haskey(population.pathogens, seq)
        return population.pathogens[seq]
    else
        return newPathogen!(
            seq, population, pathogen.type,
            parents=MVector{2,Union{Pathogen,Nothing}}([pathogen, nothing]),
            birth_time=birth_time
        )
    end
end

function recombinantSequences(
    seq_1::String, seq_2::String, num_loci::Int64,
    mean_recombination_crossovers_1::Float64, mean_recombination_crossovers_2::Float64)

    children = MVector{2,String}(seq_1, seq_2)

    if mean_recombination_crossovers_1 > 0.0 && mean_recombination_crossovers_2 > 0.0
        num_evts = zeroTruncatedPoisson(mean(
            mean_recombination_crossovers_1,
            mean_recombination_crossovers_2
        ))
        loci = rand(1:num_loci, num_evts)
        for l in loci
            temp_child_1 = children[1]
            children[1] = children[1][1:l-1] + children[2][l:end]
            children[2] = children[1][1:l-1] + temp_child_1[l:end]
        end
    else
        num_evts = 0
    end

    children = MVector{2,String}([split(seq, CHROMOSOME_SEPARATOR) for seq in children])
    parent = rand(0:1, length(children[1]))

    children = SVector{2,String}([
        join([children[parent[i]+1][i] for i in 1:length(children[1])], CHROMOSOME_SEPARATOR),
        join([children[(parent[i]!=true)+1][i] for i in 1:length(children[2])], CHROMOSOME_SEPARATOR)
    ])

    return children
end

function recombinantPathogens!(pathogen_1::Pathogen, pathogen_2::Pathogen, population::Population, birth_time::Float64)
    children = recombinantSequences(
        pathogen_1.sequence, pathogen_2.sequence, pathogen_1.type.num_loci,
        pathogen_1.mean_recombination_crossovers, pathogen_2.mean_recombination_crossovers
    )

    # for seq in children
    #     if !haskey(population.pathogens, seq)
    #         newPathogen!(
    #             seq, population, pathogen_1.type,
    #             parents=MVector{2, Union{Pathogen, Nothing}}([pathogen_1, pathogen_2]),
    #             birth_time = birth_time
    #         )
    #     end
    # end

    # return SA[population.pathogens[children[1]], population.pathogens[children[2]]]

    if !haskey(population.pathogens, children[1])
        newPathogen!(
            children[1], population, pathogen_1.type,
            parents=MVector{2,Union{Pathogen,Nothing}}([pathogen_1, pathogen_2]),
            birth_time=birth_time
        )
    end

    return population.pathogens[children[1]]
end

function generateGamete(
    sequence::String, num_loci::Int64,
    possible_alleles::String, mean_mutations_per_replication::Float64)
    gamete = join(
        [
            rand(split(homologous_chromosome_pair, HOMOLOGOUS_CHROMOSOME_SEPARATOR))
            for homologous_chromosome_pair in split(sequence, CHROMOSOME_SEPARATOR)
        ],
        CHROMOSOME_SEPARATOR
    )
    if mean_mutations_per_replication > 0.0
        gamete = mutantSequence!(
            gamete, num_loci, possible_alleles,
            mean_mutations_per_replication
        )
    end

    return gamete
end

function generateZygote(gamete_1::String, gamete_2::String)
    chromosomes_1 = split(gamete_1, CHROMOSOME_SEPARATOR)
    chromosomes_2 = split(gamete_2, CHROMOSOME_SEPARATOR)
    zygote = ""
    for chromosome_idx in eachindex(chromosomes_1)
        zygote = zygote * chromosomes_1[chromosome_idx] * HOMOLOGOUS_CHROMOSOME_SEPARATOR * chromosomes_2[chromosome_idx] * CHROMOSOME_SEPARATOR
    end

    return zygote[1:end-length(CHROMOSOME_SEPARATOR)]
end

function addPathogenToHost!(pathogen::Pathogen, host_idx::Int64, population::Population, model::Model)
    if length(population.hosts[host_idx].pathogens) == 0
        if length(population.hosts[host_idx].responses) == 0
            population.compartment_vars[INFECTED_NAIVE] += 1
            population.compartment_vars[UNINFECTED_NAIVE] -= 1
        else
            population.compartment_vars[INFECTED_IMMUNE] += 1
            population.compartment_vars[UNINFECTED_IMMUNE] -= 1
        end
    end

    push!(population.hosts[host_idx].pathogens, pathogen)
    push!(population.hosts[host_idx].pathogen_fractions, 0.0)
    population.hosts[host_idx].pathogen_weights = catCol(
        population.hosts[host_idx].pathogen_weights, zeros(Float64, NUM_PATHOGEN_EVENTS)
    )

    hostWeights!(host_idx, population, model)
end

function addResponseToHost!(response::Response, host_idx::Int64, population::Population, model::Model)
    if length(population.hosts[host_idx].responses) == 0
        if length(population.hosts[host_idx].pathogens) == 0
            population.compartment_vars[UNINFECTED_IMMUNE] += 1
            population.compartment_vars[UNINFECTED_NAIVE] -= 1
        else
            population.compartment_vars[INFECTED_IMMUNE] += 1
            population.compartment_vars[INFECTED_NAIVE] -= 1
        end
    end

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

    if length(population.hosts[host_idx].pathogens) == 0
        if length(population.hosts[host_idx].responses) == 0
            population.compartment_vars[INFECTED_NAIVE] -= 1
            population.compartment_vars[UNINFECTED_NAIVE] += 1
        else
            population.compartment_vars[INFECTED_IMMUNE] -= 1
            population.compartment_vars[UNINFECTED_IMMUNE] += 1
        end
    end

    hostWeights!(host_idx, population, model)
end

function removeResponseFromHost!(response_idx::Int64, host_idx::Int64, population::Population, model::Model)
    deleteat!(population.hosts[host_idx].responses, response_idx)
    population.hosts[host_idx].response_weights = population.hosts[host_idx].response_weights[:, begin:end.!=response_idx]

    if length(population.hosts[host_idx].responses) == 0
        if length(population.hosts[host_idx].pathogens) == 0
            population.compartment_vars[UNINFECTED_IMMUNE] -= 1
            population.compartment_vars[UNINFECTED_NAIVE] += 1
        else
            population.compartment_vars[INFECTED_IMMUNE] -= 1
            population.compartment_vars[INFECTED_NAIVE] += 1
        end
    end

    hostWeights!(host_idx, population, model)
end

function attemptInfection!(pathogen::Pathogen, host_idx::Int64, pop_idx::Int64, model::Model)
    if !(
        pathogen in
        model.populations[pop_idx].hosts[host_idx].pathogens
    ) &&
       rand() < model.populations[pop_idx].parameters.infectionProbability(
        pathogen,
        model.populations[pop_idx].hosts[host_idx]
    )

        addPathogenToHost!(
            pathogen, host_idx, model.populations[pop_idx], model
        )
    end
end

function addHostToPopulation!(new_host::Host, population::Population, model::Model)
    if length(new_host.pathogens) == 0
        if length(new_host.responses) == 0
            population.compartment_vars[UNINFECTED_NAIVE] += 1
        else
            population.compartment_vars[UNINFECTED_IMMUNE] += 1
        end
    else
        if length(new_host.responses) == 0
            population.compartment_vars[INFECTED_NAIVE] += 1
        else
            population.compartment_vars[INFECTED_IMMUNE] += 1
        end
    end

    push!(population.hosts, new_host)
    population.host_weights = catCol(population.host_weights, zeros(Float64, NUM_EVENTS))
    population.host_weights_receive = catCol(population.host_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS))
    population.host_weights_with_coefficient = catCol(population.host_weights_with_coefficient, zeros(Float64, NUM_EVENTS))
    population.host_weights_receive_with_coefficient = catCol(population.host_weights_receive_with_coefficient, zeros(Float64, NUM_CHOICE_MODIFIERS))

    # TODO: make compatible with `propagateWeightsOnAddHost!`
    for p in 1:length(model.populations)
        model.population_contact_weights_receive[model.population_dict[population.id], p] -= (
            model.population_contact_weights_receive[model.population_dict[population.id], p] /
            max(length(population.hosts) * population.parameters.constant_contact_density, 1)
        )
        model.population_contact_weights_receive_sums[p] -= (
            model.population_contact_weights_receive_sums[p] /
            max(length(population.hosts) * population.parameters.constant_contact_density, 1)
        )
        propagateWeightChanges!(
            -model.population_weights[CONTACT, p] /
            max(length(population.hosts) * population.parameters.constant_contact_density, 1),
            model.populations[p], CONTACT, model
        )
        model.population_transition_weights_receive[model.population_dict[population.id], p] -= (
            model.population_transition_weights_receive[model.population_dict[population.id], p] /
            max(length(population.hosts) * population.parameters.constant_transition_density, 1)
        )
        model.population_transition_weights_receive_sums[p] -= (
            model.population_transition_weights_receive_sums[p] /
            max(length(population.hosts) * population.parameters.constant_transition_density, 1)
        )
        propagateWeightChanges!(
            -model.population_weights[TRANSITION, p] /
            max(length(population.hosts) * population.parameters.constant_transition_density, 1),
            model.populations[p], TRANSITION, model
        )
    end

    for coef in EVENTS
        if START_COEFFICIENTS[coef] != 0.0
            propagateWeightChanges!(
                START_COEFFICIENTS[coef], length(population.hosts), population, coef, model
            )
        end
    end
    for coef in CHOICE_MODIFIERS[begin:end-1]
        if START_COEFFICIENTS[coef] != 0.0
            propagateWeightReceiveChanges!(
                START_COEFFICIENTS[coef], length(population.hosts), population, coef, model
            )
        end
    end
end

function addHostsToPopulation!(num_hosts::Int64, host_sequence::String, population::Population, model::Model)
    population.compartment_vars[UNINFECTED_NAIVE] += num_hosts

    # update matrices
    population.host_weights = hcat(population.host_weights, zeros(Float64, NUM_EVENTS, num_hosts))
    population.host_weights_receive = hcat(population.host_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS, num_hosts))
    population.host_weights_with_coefficient = hcat(population.host_weights_with_coefficient, zeros(Float64, NUM_EVENTS, num_hosts))
    population.host_weights_receive_with_coefficient = hcat(population.host_weights_receive_with_coefficient, zeros(Float64, NUM_CHOICE_MODIFIERS, num_hosts))

    num_starting_hosts = length(population.hosts)
    for i in 1:num_hosts
        push!(population.hosts, Host(
            length(population.hosts) + 1,
            MVector{2,Union{Host,Nothing}}([nothing, nothing]),
            model.time,
            host_sequence,
            population.parameters.host_mean_mutations_per_replication * population.parameters.hostMutationCoefficient(
                host_sequence
            ),
            population.parameters.host_mean_recombination_crossovers * population.parameters.hostRecombinationCoefficient(
                host_sequence
            ),
            Vector{Pathogen}(undef, 0), Vector{Response}(undef, 0),
            Vector{Float64}(undef, 0),
            Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, 0),
            Matrix{Float64}(undef, NUM_RESPONSE_EVENTS, 0),
            MVector{NUM_NON_SAMPLING_COEFFICIENTS,Float64}(
                START_COEFFICIENTS[NON_SAMPLING_COEFFICIENTS[begin]:NON_SAMPLING_COEFFICIENTS[end]]
            ),
        ))

        propagateWeightsOnAddHost!(num_starting_hosts + i, population, model)
    end
end

function removeHostFromPopulation!(host_idx::Int64, population::Population, model::Model)
    if length(population.hosts[host_idx].pathogens) == 0
        if length(population.hosts[host_idx].responses) == 0
            population.compartment_vars[UNINFECTED_NAIVE] -= 1
        else
            population.compartment_vars[UNINFECTED_IMMUNE] -= 1
        end
    else
        if length(population.hosts[host_idx].responses) == 0
            population.compartment_vars[INFECTED_NAIVE] -= 1
        else
            population.compartment_vars[INFECTED_IMMUNE] -= 1
        end
    end

    for coef in EVENTS
        if population.host_weights[host_idx, coef] != 0.0
            change = -population.host_weights[host_idx, coef]
            if coef == CONTACT
                population.contact_sum += change
                change = change * model.population_contact_weights_receive_sums[model.population_dict[population.id]]
            elseif coef == TRANSITION
                population.transition_sum += change
                change = change * model.population_transition_weights_receive_sums[model.population_dict[population.id]]
            end
            propagateWeightChanges!(
                change * population.parameters.base_coefficients[coef], population, coef, model
            )
        end
    end
    for coef in CHOICE_MODIFIERS[begin:end-1]
        if population.host_weights[host_idx, coef] != 0.0
            propagateWeightReceiveChanges!(
                -population.host_weights[host_idx, coef] * population.parameters.base_coefficients[coef], population, coef, model
            )
        end
    end

    for p in 1:length(model.populations)
        model.population_contact_weights_receive[model.population_dict[population.id], p] += (
            model.population_contact_weights_receive[model.population_dict[population.id], p] /
            max(length(population.hosts) * population.parameters.constant_contact_density, 1)
        )
        model.population_contact_weights_receive_sums[p] += (
            model.population_contact_weights_receive_sums[p] /
            max(length(population.hosts) * population.parameters.constant_contact_density, 1)
        )
        propagateWeightChanges!(
            model.population_weights[CONTACT, p] /
            max(length(population.hosts) * population.parameters.constant_contact_density, 1),
            model.populations[p], CONTACT, model
        )
        model.population_transition_weights_receive[model.population_dict[population.id], p] += (
            model.population_transition_weights_receive[model.population_dict[population.id], p] /
            max(length(population.hosts) * population.parameters.constant_transition_density, 1)
        )
        model.population_transition_weights_receive_sums[p] += (
            model.population_transition_weights_receive_sums[p] /
            max(length(population.hosts) * population.parameters.constant_transition_density, 1)
        )
        propagateWeightChanges!(
            model.population_weights[TRANSITION, p] /
            max(length(population.hosts) * population.parameters.constant_transition_density, 1),
            model.populations[p], TRANSITION, model
        )
    end

    population.host_weights = population.host_weights[1:end.!=host_idx, 1:end.!=host_idx]
    population.host_weights_receive = population.host_weights_receive[1:end.!=host_idx, 1:end.!=host_idx]
    population.host_weights_with_coefficient = population.host_weights_with_coefficient[1:end.!=host_idx, 1:end.!=host_idx]
    population.host_weights_receive_with_coefficient = population.host_weights_receive_with_coefficient[1:end.!=host_idx, 1:end.!=host_idx]

    deleteat!(population.hosts, host_idx)
end

function setPopulationContactCoefficient!(pop_idx_1::Int64, pop_idx_2::Int64, coefficient::Float64, model::Model)
    model.populations[pop_idx_1].population_contact_coefficients[pop_idx_2] = coefficient
    change = (
        model.populations[pop_idx_1].population_contact_coefficients[pop_idx_2] *
        model.population_weights_receive[RECEIVE_CONTACT-CHOICE_MODIFIERS[1]+1, pop_idx_2] /
        max(
            length(model.populations[pop_idx_2].hosts) *
            model.populations[pop_idx_2].parameters.constant_contact_density,
            1.0
        )) - model.population_contact_weights_receive[pop_idx_2, pop_idx_1]
    updatePopulationContactWeightReceiveMatrix!(pop_idx_1, pop_idx_2, change, model)
end

function setPopulationTransitionCoefficient!(pop_idx_1::Int64, pop_idx_2::Int64, coefficient::Float64, model::Model)
    model.populations[pop_idx_1].population_transition_coefficients[pop_idx_2] = coefficient
    change = (
        model.populations[pop_idx_1].population_transition_coefficients[pop_idx_2] *
        model.population_weights_receive[RECEIVE_TRANSITION-CHOICE_MODIFIERS[1]+1, pop_idx_2] /
        max(
            length(model.populations[pop_idx_2].hosts) *
            model.populations[pop_idx_2].parameters.constant_transition_density,
            1.0
        )) - model.population_transition_weights_receive[pop_idx_2, pop_idx_1]
    updatePopulationTransitionWeightReceiveMatrix!(pop_idx_1, pop_idx_2, change, model)
end

# Model events

function establishMutant!(model::Model, rand_n::Float64)
    pathogen_idx, host_idx, pop_idx = choosePathogen(MUTANT_ESTABLISHMENT, model, rand_n)

    mut = mutantPathogen!(
        model.populations[pop_idx].hosts[host_idx].pathogens[pathogen_idx],
        model.populations[pop_idx], model.time
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
    pathogen_idx, host_idx, pop_idx = choosePathogen(RESPONSE_ACQUISITION, model, rand_n)

    responses = model.populations[pop_idx].parameters.developResponses(
        model.populations[pop_idx].hosts[host_idx].pathogens[pathogen_idx],
        model.populations[pop_idx].hosts[host_idx],
        model.populations[pop_idx].responses,
        model.populations[pop_idx].parameters.response_types,
        model.time
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
            model.populations[pop_idx],
            model.time
        )

        attemptInfection!(recombinant, host_idx, pop_idx, model)
    end
end

function hostContact!(model::Model, rand_n::Float64)
    host_idx_1, pop_idx_1, rand_n = chooseHost(CONTACT, model, rand_n)
    pop_idx_2, rand_n = choosePopulationReceiveContact(pop_idx_1, model, rand_n)
    host_idx_2, rand_n = chooseHost(pop_idx_2, RECEIVE_CONTACT, model, rand_n)

    if host_idx_1 != host_idx_2 || pop_idx_1 != pop_idx_2
        host1 = model.populations[pop_idx_1].hosts[host_idx_1]
        # inocula = MVector{length(model.populations[pop_idx_1].hosts[host_idx_1].pathogens),Int64}([
        inocula = Vector{Int64}([
            pois_rand(
                host1.pathogens[p_idx].mean_effective_inoculum *
                host1.pathogen_fractions[p_idx]
            )
            for p_idx in 1:length(host1.pathogens)
        ])
        if all(inocula .== 0.0) # had originally written !any(inocula .!= 0.0), don't know if that's actually faster
            idx, rand_n = randChoose(rand_n, host1.pathogen_fractions, 1.0)
            inocula[idx] += 1
        end

        for p_idx in eachindex(inocula)
            if inocula[p_idx] > 0
                mut_prob = min(
                    1.0, 1.0 - exp(
                        -host1.pathogens[p_idx].mean_mutations_per_replication
                    )
                )
                # probability from Poisson PMF with k=0
                rec_prob = 0.0
                if length(host1.pathogens) > 1
                    rec_prob = min(
                        1.0, 1.0 - exp(
                            -host1.pathogens[p_idx].mean_recombination_crossovers
                        )
                    )
                    # probability from Poisson PMF with k=0
                end
                num_mut = binomial(inocula[p_idx], mut_prob * (1.0 - rec_prob))
                num_rec = binomial(inocula[p_idx], rec_prob * (1.0 - mut_prob))
                num_mut_rec = binomial(inocula[p_idx], mut_prob * rec_prob)
                if inocula[p_idx] - num_mut - num_rec - num_mut_rec > 0
                    attemptInfection!(
                        host1.pathogens[p_idx],
                        host_idx_2, pop_idx_2, model
                    )
                end
                for _ in 1:num_mut
                    attemptInfection!(
                        mutantPathogen!(
                            host1.pathogens[p_idx],
                            model.populations[pop_idx_1],
                            model.time
                        ), host_idx_2, pop_idx_2, model
                    )
                end
                for _ in 1:num_rec
                    p_idx_2, rand_n = choosePathogen(host_idx_1, pop_idx_1, RECOMBINANT_ESTABLISHMENT, model, rand())

                    if p_idx != p_idx_2 && (
                        host1.pathogens[p_idx].type ==
                        host1.pathogens[p_idx_2].type)

                        attemptInfection!(
                            recombinantPathogens!(
                                host1.pathogens[p_idx],
                                host1.pathogens[p_idx_2],
                                model.populations[pop_idx_1],
                                model.time
                            ), host_idx_2, pop_idx_2, model
                        )
                    end
                end
                for _ in 1:num_mut_rec
                    p_idx_2, rand_n = choosePathogen(host_idx_1, pop_idx_1, RECOMBINANT_ESTABLISHMENT, model, rand())

                    if p_idx != p_idx_2
                        recombinant = recombinantPathogens!(
                            host1.pathogens[p_idx],
                            host1.pathogens[p_idx_2],
                            model.populations[pop_idx_1],
                            model.time
                        )

                        attemptInfection!(
                            mutantPathogen!(recombinant, model.populations[pop_idx_1], model.time),
                            host_idx_2, pop_idx_2, model
                        )
                    else
                        attemptInfection!(
                            mutantPathogen!(
                                host1.pathogens[p_idx],
                                model.populations[pop_idx_1],
                                model.time
                            ),
                            host_idx_2, pop_idx_2, model
                        )
                    end
                end
            end
        end
    end
end

function loseResponse!(model::Model, rand_n::Float64)
    response_idx, host_idx, pop_idx = chooseResponse(RESPONSE_LOSS, model, rand_n)

    removeResponseFromHost!(
        response_idx, host_idx, model.populations[pop_idx], model
    )
end

function birth!(model::Model, rand_n::Float64)
    host_idx, pop_idx, rand_n = chooseHost(BIRTH, model, rand_n)
    parents = MVector{2,Union{Nothing,Host}}[model.populations[pop_idx].hosts[host_idx], nothing]
    if model.populations[pop_idx].host_sexual_reproduction
        host_idx_2, pop_idx_2, rand_n = chooseHost(BIRTH, model, rand_n)
        parents[2] = model.populations[pop_idx_2].hosts[host_idx_2]
        if !model.populations[pop_idx].sexualCompatibility(parents[1], parents[2])
            parents = MVector{2,Union{Nothing,Host}}[nothing, nothing]
        end
    end

    if !isnothing(parent[1])
        child_sequence = ""

        parent_gametes = MVector{2,Union{String,Nothing}}[
            generateGamete(
                parent[1].sequence,
                model.populations[pop_idx].host_num_loci,
                model.populations[pop_idx].host_possible_alleles,
                parents[1].host_mean_mutations_per_replication
            ), nothing
        ]

        if model.populations[pop_idx].diploid_host
            parent_gametes[1] = recombinantSequences(
                parent_gametes[1],
                generateGamete(
                    parent[1].sequence,
                    model.populations[pop_idx].host_num_loci,
                    model.populations[pop_idx].host_possible_alleles,
                    parents[1].host_mean_mutations_per_replication
                ),
                model.populations[pop_idx].host_num_loci,
                parents[1].mean_recombination_crossovers,
                parents[1].mean_recombination_crossovers
            )
            if !isnothing(parent[2])
                parent_gametes[2] = recombinantSequences(
                    parent_gametes[2],
                    generateGamete(
                        parent[2].sequence,
                        model.populations[pop_idx_2].host_num_loci,
                        model.populations[pop_idx_2].host_possible_alleles,
                        parents[2].host_mean_mutations_per_replication
                    ),
                    model.populations[pop_idx_2].host_num_loci,
                    parents[2].mean_recombination_crossovers,
                    parents[2].mean_recombination_crossovers
                )
            end
            child_sequence = generateZygote(parent_gametes[1], parent_gametes[2])
        else
            if !isnothing(parent[2]) && parents[1].host_mean_recombination_crossovers > 0.0 && parents[2].host_mean_recombination_crossovers > 0.0
                child_sequence = recombinantSequences(
                    parent_gametes[1],
                    generateGamete(
                        parent[2].sequence,
                        model.populations[pop_idx].host_num_loci,
                        model.populations[pop_idx].host_possible_alleles,
                        parents[2].host_mean_mutations_per_replication
                    ),
                    model.populations[pop_idx].host_num_loci,
                    parents[1].mean_recombination_crossovers,
                    parents[2].mean_recombination_crossovers
                )
            else
                child_sequence = parent_gametes[1]
            end
        end

        jOpqua.newHost!(child_sequence, model.populations[pop_idx], model, parents=parents, birth_time=model.time)

        for parent in parents
            if !isnothing(parent)
                for response in parent.responses
                    if response.type.inherit_response > 0.0 && rand() < response.type.inherit_response
                        addResponseToHost!(
                            response, length(model.populations[pop_idx].hosts),
                            model.populations[pop_idx], model
                        )
                    end
                end

                for pathogen_idx in 1:length(parent.pathogens)
                    if (parent.pathogens[pathogen_idx].type.vertical_transmission > 0.0 &&
                        rand() < parent.pathogens[pathogen_idx].vertical_transmission_coefficient *
                                 parent.pathogens[pathogen_idx].type.verticalTransmissionCoefficient(
                                     parent.pathogens[pathogen_idx].sequence
                                 ) *
                                 parent.pathogen_fractions[pathogen_idx])
                        attemptInfection!(
                            parent.pathogens[pathogen_idx],
                            length(model.populations[pop_idx].hosts), pop_idx, model
                        )
                    end
                end
            end
        end
    end
end

function death!(model::Model, rand_n::Float64)
    host_idx, pop_idx, rand_n = chooseHost(DEATH, model, rand_n)

    model.populations[pop_idx].compartment_vars[DEAD] += 1
    # done here because removeHostFromPopulation doesn't advance deaths,
    # but does remove from other vars

    removeHostFromPopulation!(host_idx, model.populations[pop_idx], model)
end

function transition!(model::Model, rand_n::Float64)
    host_idx, pop_idx_1, rand_n = chooseHost(TRANSITION, model, rand_n)
    pop_idx_2, rand_n = choosePopulationReceiveTransition(pop_idx_1, model, rand_n)

    addHostToPopulation!(model.populations[pop_idx_1].hosts[host_idx], model.populations[pop_idx_2], model)
    removeHostFromPopulation!(host_idx, model.populations[pop_idx_1], model)
end

const EVENT_FUNCTIONS = SA[
    establishMutant!, clearPathogen!, acquireResponse!,
    establishRecombinant!, hostContact!, loseResponse!,
    birth!, death!, transition!
]
