using DataFrames
using CSV

function saveCompartments(output::Output, file_name::String)
    pop_ids = collect(keys(output.compartment_vars))
    df = DataFrame(
        hcat([output.compartment_vars[id] for id in pop_ids]...)',
        COMPARTMENT_LABELS
    )
    df[!, "Population"] = repeat(pop_ids, inner=[length(output.time)])
    df[!, "Time"] = output.time

    df = df[!, ["Time", "Population", COMPARTMENT_LABELS...]]

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
            )
        end
    end

    out = "Time, Population, Host, Pathogen_id, Pathogen_sequence, Pathogen_parents, Response_id, Response_imprinted_sequence, Response_matured_sequence, Response_parents\n" * join(out_strs, "\n")

    # fname = "foobar.csv"
    # dirpath = "/tmp"
    # fpath = joinpath(dirpath, fname)

    open(file_name, "w") do file
        write(file, out)
    end

    return DataFrame(CSV.File(file_name))
end

function saveComposition(
        data::DataFrame, file_name::String;
        type_of_composition::String="Pathogens", num_top_sequences::Int64=-1,
        genomic_positions::Vector{String}=[], track_specific_sequences::Vector{String}=[],
        fitness_function::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing)

    return -1
end

function saveComposition(
        output::Output, file_name::String;
        type_of_composition::String="Pathogens", num_top_sequences::Int64=-1,
        genomic_positions::Vector{String}=[], track_specific_sequences::Vector{String}=[],
        fitness_function::Union{Nothing,FunctionWrapper{Float64,Tuple{String}}}=nothing)
    return saveComposition(
        saveHistory(output, file_name*".csv"), file_name,
        type_of_composition=type_of_composition, num_top_sequences=num_top_sequences,
        genomic_positions=genomic_positions, track_specific_sequences=track_specific_sequences,
        fitness_function=fitness_function
    )
end
