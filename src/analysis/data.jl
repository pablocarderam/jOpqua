using DataFrames
using StatsBase
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

function saveComposition(
    data::DataFrame, file_name::String;
    populations::Vector{String}=Vector{String}(undef,0),
    type_of_composition::String="Pathogen_sequence", num_top_sequences::Int64=-1,
    genomic_positions::Vector{Int64}=Vector{Int64}(undef,0), track_specific_sequences::Vector{String}=Vector{String}(undef,0),
    top_host_sequence_function::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing)

    dat = deepcopy(data)
    if length(populations) > 0
        dat = dat[in.(dat["Population"], Ref(Set(populations))), :]
    end

    if length(genomic_positions) > 0
        dat[:,type_of_composition] = [
            join(
                [
                    join(
                        [
                            s[genomic_positions[i]] for i in genomic_positions
                        ], ""
                    ) for s in split(host_seqs, WITHIN_HOST_SEPARATOR)
                ], WITHIN_HOST_SEPARATOR
            ) for host_seqs in dat[:,type_of_composition]
        ]
    end

    if !isnothing(top_host_sequence_function)
        dat[:,"seqpop"] = string.(dat[:,type_of_composition], ":::", dat[:,"Population"])
        unique_seqpop = unique(dat[:,"patpop"])
        for combination in unique_seqpop
            gen_str = split(combination, ":::")[1]
            if gen_str != ""
                genomes = split(gen_str, WITHIN_HOST_SEPARATOR)
                values = [top_host_sequence_function(p) for p in genomes]
                top_pathogen = genomes[findmax(values)[2]]
                dat[:,"patpop"] = replace(dat[:,"patpop"], combination => top_pathogen)
            else
                dat[:,"patpop"] = replace(dat[:,"patpop"], combination => "")
            end
        end
        dat[:,type_of_composition] = dat[:,"patpop"]
    end

    all_seqs = split(
        join(dat[(dat[:,type_of_composition].!=""), type_of_composition], WITHIN_HOST_SEPARATOR),
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
        zeros(Int64, length(unique(dat[:,"Time"])), length(seqs_to_track) + 1),
        [seqs_to_track..., "Other"]
    )
    out[:,"Time"] = unique(dat[:,"Time"])
    out = out[:,["Time", seqs_to_track..., "Other"]]

    for i in 1:length(dat[:,"Time"])
        for seq in split(dat[i,type_of_composition], WITHIN_HOST_SEPARATOR)
            if seq in seqs_to_track
                out[(out[:,"Time"].==dat[i,"Time"]), seq] .+= 1
            elseif seq != ""
                out[(out[:,"Time"].==dat[i,"Time"]), "Other"] .+= 1
            end
        end
    end

    CSV.write(file_name, out)

    return out
end

function saveComposition(
    output::Output, file_name::String;
    populations::Vector{String}=Vector{String}(undef,0),
    type_of_composition::String="Pathogen_sequence", num_top_sequences::Int64=-1,
    genomic_positions::Vector{Int64}=Vector{Int64}(undef,0), track_specific_sequences::Vector{String}=Vector{String}(undef,0),
    top_host_sequence_function::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing)
    return saveComposition(
        saveHistory(output, file_name * "_history.csv"), file_name,
        populations=populations,
        type_of_composition=type_of_composition, num_top_sequences=num_top_sequences,
        genomic_positions=genomic_positions, track_specific_sequences=track_specific_sequences,
        top_host_sequence_function=top_host_sequence_function
    )
end
