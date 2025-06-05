using DataFrames
using StatsBase
using Distances
using CSV

function saveCompartments(output::Output, file_name::String)
    pop_ids = collect(keys(output.compartment_vars))
    df = DataFrame(
        hcat([output.compartment_vars[id] for id in pop_ids]...)',
        COMPARTMENT_LABELS
    )
    df[:, "Population"] = repeat(pop_ids, inner=[length(output.time)])
    df[:, "Time"] = output.time

    df = df[:, ["Time", "Population", COMPARTMENT_LABELS...]]

    CSV.write(file_name, df)

    return df
end

function saveHistory(output::Output, file_name::String)
    pop_ids = collect(keys(output.host_samples))
    out_strs = repeat([""], inner=[length(pop_ids)])
    for i in eachindex(pop_ids)
        for t in eachindex(output.time)
            out_strs[i] = out_strs[i] * join(
                              [
                                  join(
                                      [
                                          string(output.time[t]), string(pop_ids[i]),
                                          output.host_samples[pop_ids[i]][h, t].id,
                                          output.host_samples[pop_ids[i]][h, t].sequence,
                                          "\"" * join([p.type.id for p in output.host_samples[pop_ids[i]][h, t].pathogens], WITHIN_HOST_SEPARATOR) * "\"",
                                          "\"" * join([p.sequence for p in output.host_samples[pop_ids[i]][h, t].pathogens], WITHIN_HOST_SEPARATOR) * "\"",
                                          "\"" * join([
                                                  join([isnothing(a) ? "None" : a.sequence for a in p.parents], PARENT_SEPARATOR)
                                                  for p in output.host_samples[pop_ids[i]
                                                  ][h, t].pathogens
                                              ], WITHIN_HOST_SEPARATOR) * "\"",
                                          "\"" * join([r.type.id for r in output.host_samples[pop_ids[i]][h, t].responses], WITHIN_HOST_SEPARATOR) * "\"",
                                          "\"" * join([r.imprinted_pathogen.sequence for r in output.host_samples[pop_ids[i]][h, t].responses], WITHIN_HOST_SEPARATOR) * "\"",
                                          "\"" * join([isnothing(r.matured_pathogen) ? "None" : r.matured_pathogen.sequence for r in output.host_samples[pop_ids[i]][h, t].responses], WITHIN_HOST_SEPARATOR) * "\"",
                                          "\"" * join([
                                                  join([isnothing(a) ? "None" : a.id for a in r.parents], PARENT_SEPARATOR)
                                                  for r in output.host_samples[pop_ids[i]][h, t].responses
                                              ], WITHIN_HOST_SEPARATOR) * "\"",
                                      ], ","
                                  )
                                  for h in 1:size(output.host_samples[pop_ids[i]])[1]
                              ], "\n"
                          ) * "\n"
        end
    end

    out = "Time,Population,Host,Host_sequence,Pathogen_id,Pathogen_sequence,Pathogen_parents,Response_id,Response_imprinted_sequence,Response_matured_sequence,Response_parents\n" * join(out_strs, "\n")

    open(file_name, "w") do file
        write(file, out)
    end

    return DataFrame(CSV.File(file_name))
end

function ancestors(output::Output, population_id::String, sequences::Vector{String})
    pats = [output.model.populations[output.model.population_dict[population_id]].pathogens[s] for s in sequences]
    seqs = Vector{String}(undef, 0)
    times = Vector{Float64}(undef, 0)

    pat = pats[1]

    function addPathogen!(p::Pathogen)
        push!(seqs, p.sequence)
        push!(times, p.birth_time)
        if !(isnothing(p.parents[1])) && !(p.parents[1].sequence in seqs)
            addPathogen!(p.parents[1])
        end
        if !(isnothing(p.parents[2])) && !(p.parents[2].sequence in seqs)
            addPathogen!(p.parents[2])
        end
    end

    for p in pats
        pat = p
        if !(isnothing(p)) && !(p.sequence in seqs)
            addPathogen!(p)
        end
    end

    return DataFrame(sequence=seqs, time=times)
end

function saveComposition(
    data::DataFrame, file_name::String;
    populations::Vector{String}=Vector{String}(undef, 0),
    type_of_composition::String="Pathogen_sequence", num_top_sequences::Int64=-1,
    genomic_positions::Vector{Int64}=Vector{Int64}(undef, 0), track_specific_sequences::Vector{String}=Vector{String}(undef, 0),
    top_host_sequence_function::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing)

    dat = coalesce.(data, "")
    if length(populations) > 0
        dat = dat[in.(dat["Population"], Ref(Set(populations))), :]
    end

    if length(genomic_positions) > 0
        dat[:, type_of_composition] = [
            join(
                [
                    join(
                        [
                            s[genomic_positions[i]] for i in genomic_positions
                        ], ""
                    ) for s in split(host_seqs, WITHIN_HOST_SEPARATOR)
                ], WITHIN_HOST_SEPARATOR
            ) for host_seqs in dat[:, type_of_composition]
        ]
    end

    if !isnothing(top_host_sequence_function)
        dat[:, "seqpop"] = string.(dat[:, type_of_composition], ":::", dat[:, "Population"])
        unique_seqpop = unique(dat[:, "patpop"])
        for combination in unique_seqpop
            gen_str = split(combination, ":::")[1]
            if gen_str != ""
                genomes = split(gen_str, WITHIN_HOST_SEPARATOR)
                values = [top_host_sequence_function(p) for p in genomes]
                top_pathogen = genomes[findmax(values)[2]]
                dat[:, "patpop"] = replace(dat[:, "patpop"], combination => top_pathogen)
            else
                dat[:, "patpop"] = replace(dat[:, "patpop"], combination => "")
            end
        end
        dat[:, type_of_composition] = dat[:, "patpop"]
    end

    all_seqs = split(
        join(dat[(dat[:, type_of_composition].!=""), type_of_composition], WITHIN_HOST_SEPARATOR),
        WITHIN_HOST_SEPARATOR
    )
    counts = countmap(all_seqs)
    unique_seqs = collect(keys(counts))
    top_sequences = DataFrame()
    top_sequences[:, "Sequence"] = unique_seqs
    top_sequences[:, "Frequency"] = [counts[seq] for seq in unique_seqs]
    if num_top_sequences < 0 || length(top_sequences[:, "Sequence"]) < num_top_sequences
        num_top_sequences = length(top_sequences[:, "Sequence"])
    end

    seqs_to_track = track_specific_sequences
    for seq in top_sequences[1:num_top_sequences, "Sequence"]
        if !(seq in seqs_to_track)
            push!(seqs_to_track, seq)
        end
    end

    out = DataFrame(
        zeros(Int64, length(unique(dat[:, "Time"])), length(seqs_to_track) + 1),
        [seqs_to_track..., "Other"]
    )
    out[:, "Time"] = unique(dat[:, "Time"])
    out = out[:, ["Time", seqs_to_track..., "Other"]]

    for i in 1:length(dat[:, "Time"])
        for seq in split(dat[i, type_of_composition], WITHIN_HOST_SEPARATOR)
            if seq in seqs_to_track
                out[(out[:, "Time"].==dat[i, "Time"]), seq] .+= 1
            elseif seq != ""
                out[(out[:, "Time"].==dat[i, "Time"]), "Other"] .+= 1
            end
        end
    end

    out[:, "Total"] = sum([out[:, c] for c in names(out)[2:end]])

    CSV.write(file_name, out)

    return out
end

function saveComposition(
    output::Output, file_name::String;
    populations::Vector{String}=Vector{String}(undef, 0),
    type_of_composition::String="Pathogen_sequence", num_top_sequences::Int64=-1,
    genomic_positions::Vector{Int64}=Vector{Int64}(undef, 0), track_specific_sequences::Vector{String}=Vector{String}(undef, 0),
    top_host_sequence_function::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing)
    return saveComposition(
        saveHistory(output, file_name * "_history.csv"), file_name,
        populations=populations,
        type_of_composition=type_of_composition, num_top_sequences=num_top_sequences,
        genomic_positions=genomic_positions, track_specific_sequences=track_specific_sequences,
        top_host_sequence_function=top_host_sequence_function
    )
end

function addToPhylogenyDicts!(
    nodes_dict::Dict, node_info::Dict, pop::Population,
    node::Union{Pathogen,Response,Host}, node_id::String, info_separator::String)
    if !haskey(nodes_dict, node)
        nodes_dict[node] = []
    end
    if !haskey(node_info, node)
        node_info[node] = pop.id * info_separator * node.type.id * info_separator * node_id
    end
    if !isnothing(node.parents[1]) # using only first parent for now because I'm not sure how to handle recombination
        if haskey(nodes_dict, node.parents[1])
            push!(nodes_dict[node.parents[1]], node)
        else
            nodes_dict[node.parents[1]] = [node]
        end
    end
    # for parent in node.parents
    #     if !isnothing(parent)
    #         if haskey(nodes_dict, parent)
    #             push!(nodes_dict[parent], node)
    #         else
    #             nodes_dict[parent] = [node]
    #         end
    #     end
    # end
end

function saveNewick(
    model::Model, file_name::String; type::String="Pathogens",
    info_separator::String="|", branch_length::String="Distance")
    nodes_dict = Dict()
    node_info = Dict()
    for pop in model.populations
        if type == "Pathogens"
            for (id, node) in pop.pathogens
                addToPhylogenyDicts!(nodes_dict, node_info, pop, node, id, info_separator)
            end
        elseif type == "Responses"
            for (id, node) in pop.responses
                addToPhylogenyDicts!(nodes_dict, node_info, pop, node, join(id, info_separator), info_separator)
            end
        elseif type == "Host"
            for node in pop.hosts
                addToPhylogenyDicts!(nodes_dict, node_info, pop, node, node.sequence, info_separator)
            end
        end
    end

    trees = Dict()
    for (node, children) in nodes_dict
        node_str = ""
        if length(children) > 0
            for child in children
                child_contents = ""
                if haskey(trees, node_info[child])
                    child_contents = trees[node_info[child]]
                    delete!(trees, node_info[child])
                end
                if branch_length == "Distance"
                    branch_length_number = hamming(node.sequence, child.sequence)
                elseif branch_length == "Time"
                    branch_length_number = child.birth_time - node.birth_time
                end
                node_str = node_str * child_contents * node_info[child] * ":" * string(branch_length_number) * ","
            end
            node_str = "(" * node_str[1:end-1] * ")"
        end

        found = false
        tree_i = 0
        roots = collect(keys(trees))
        while !found && tree_i < length(roots)
            tree_i += 1
            if occursin(node_info[node], trees[roots[tree_i]])
                found = true
            end
        end
        if length(children) > 0
            if found
                trees[roots[tree_i]] = trees[roots[tree_i]][1:findfirst(node_info[node], trees[roots[tree_i]])[1]-1] * node_str * trees[roots[tree_i]][findfirst(node_info[node], trees[roots[tree_i]])[1]:end]
            else
                trees[node_info[node]] = node_str
            end
        elseif !found
            trees[node_info[node]] = ""
        end
    end

    out = []

    for (root, stem) in trees
        push!(out, stem * root * ":0.0;")
        open(file_name, "w") do file
            write(file, out[end])
        end
    end

    return out
end

function saveNewick(output::Output, file_name::String; type::String="Pathogens", info_separator::String="|", branch_length::String="Distance")
    return saveNewick(output.model, file_name, type=type, info_separator=info_separator, branch_length=branch_length)
end
