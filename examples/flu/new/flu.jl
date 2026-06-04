# Run from base jOpqua directory as
# julia --project=. examples/flu/new/flu.jl
# unless viewing flamegraph, then run from console

# using Revise
using jOpqua

using StaticArrays
using Random

using Distances

using BenchmarkTools
using ProfileView


using SHA
using Statistics
using MultivariateStats

using DataFrames
using Plots
using Plots.PlotMeasures
using CSV

using Base.Threads

const StringOrSubString = Union{String,SubString{String}}

# model constants
const ha_sn89::String = "MKAKLLVLLCAFTATDADTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCRLKGIAPLQLGNCSIAGWILGNPECESLFSKESWSYIAETPNSENGTCYPGYFADYEELREQLSSVSSFERFEIFPKESSWPNHTVTKGVTAACSHNGKSSFYRNLLWLTEKNGLYPNLSKSYVNNKEKEVLVLWGVHHPSNIGDQRAIYHTENAYVSVVSSHYSRRFTPEIAKRPKVRGQEGRINYYWTLLEPGDTIIFEANGNLIAPWYAFALSRGFGSGIITSNASMDECDAKCQTPQGAINSSLPFQNVHPVTIGECPKYVRSTKLRMVTGLRNIPSVQSRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAINGITNKVNSVIEKMNTQFTAVGKEFNKLERRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVKNLYEKVKSQLKNNAKEIGNGCFEFYHKCNNECMESVKNGTYDYPKYSEESKLNREKIDGVKLESMGVYQILAIYSTVASSLVLLVSLGAISFWMCSNGSLQCRICI"
const ha_sn89_mut = "MKAKLLVLLCAFTATDADTICIGYHANTSTDTVDTVLEKNVTVTHNVNLLEDSHNGKLCRLKGIAPLQLGNCSIAGWILGNPECESLFSKESWSYIAETPNSENATCEYGTCADYEELREQLSSVSSFERFEPFPNESSWPNHLVTPGVTAACSHNGKASFYRNLLWLTEKNGLYPNLPKSYVNNKEKEVLVLWGVHHPDNIGDQRAIYHTENAYVSVVSSHYSRRFTPEIAKRPKVRGQEMEINYYWTLLEPGDTIIFEANGNLIAPWDAFALSRGFGSGIITSNASMDECDAKCQTPQGAINSSLPFQNVHPVTIGECPKYGRSTKLRMVTGLRNIPSDQSRGLFGAIAGFCEGGWTEMIDGWYGYHHQNEQGSGYAADQKSTQNGITGITNKVNSVIEKMNTQFTAQGKEFNKLERRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVKNCYEKVKSQLKNNAKEIGNGCFEFYHKCNNECMESVKNGTYDYPKYSEESKLNREKIDGVKLESMTVYQILAIYSTVASSLVLLVSAGAISFWMCSNGSLKCRICI"
# const rbs_residues = collect(range(140, 210))
# const rbs_epitope_key_residues = [144, 147, 159, 169, 170, 171, 172, 173, 203, 207]

const necessary_residues = collect(range(1, length(ha_sn89) - 25))
const evasion_residues = collect(range(length(ha_sn89) - 25 + 1, length(ha_sn89)))

function sequence_hamming(a::StringOrSubString, b::StringOrSubString; distance_function::Function=(x::Char, y::Char) -> x == y ? 0 : 1)
    distance = 0
    for i in eachindex(a)
        distance += distance_function(a[i], b[i])
    end
    return distance
end

function crossImmunity(
    seq1::StringOrSubString, seq2::StringOrSubString;
    epitope_residues::Array{Int64}=evasion_residues,
    distance_half_all_residues::Float64=5.0, hill_coef_all_residues=2.0,
    distance_half_epitope_residues::Float64=length(evasion_residues) * 4.4, hill_coef_epitope_residues=2.0)

    distance_key_residues = 0.0
    for p in epitope_residues
        distance_key_residues += seq1[p] != seq2[p]
    end

    # return (sequence_hamming(seq1, seq2) < 5) * 0.2 + 0.8

    return (1.0 - (1.0 -
                   (1.0 - jOpqua.hillFunction(
        Float64(sequence_hamming(seq1, seq2)),
        distance_half_all_residues, hill_coef_all_residues
    )) * (1.0 - jOpqua.hillFunction(
        distance_key_residues,
        distance_half_epitope_residues, hill_coef_epitope_residues
    ))
    ) *
                  1.0#genomeRand(seq2)
    ) *
           0.2 + 0.8
    # return seq1 == seq2
    # return 1.0
    # return 0.8
end

println(("TEST SAME: ", 1.0 - crossImmunity(ha_sn89, ha_sn89)))
println(("TEST MUT: ", 1.0 - crossImmunity(ha_sn89, ha_sn89_mut)))

function proteinFitness(
    seq::StringOrSubString, seq_wt::StringOrSubString;
    functional_site::Array{Int64}=necessary_residues,
    distance_half_all_residues::Float64=length(ha_sn89) * 0.5, hill_coef_all_residues=3.0,
    distance_half_functional_site_residues::Float64=length(necessary_residues) * 0.5, hill_coef_functional_site_residues=3.0)

    if occursin("*", seq)
        return 0.0
    else
        distance_key_residues = 0
        for p in functional_site
            distance_key_residues += seq[p] != seq_wt[p]
        end

        # return (1.0 - jOpqua.hillFunction(
        #     Float64(sequence_hamming(seq, seq_wt)),
        #     distance_half_all_residues, hill_coef_all_residues
        # )) * (1.0 - jOpqua.hillFunction(
        #     distance_key_residues,
        #     distance_half_functional_site_residues, hill_coef_functional_site_residues
        # ))
        # return genomeRand(seq)
        # return distance_key_residues < 1
        return 1.0
    end
end

function analyze(output::jOpqua.Output, t_vec::Vector{Float64}; ha_sn89::String="MKAKLLVLLCAFTATDADTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCRLKGIAPLQLGNCSIAGWILGNPECESLFSKESWSYIAETPNSENGTCYPGYFADYEELREQLSSVSSFERFEIFPKESSWPNHTVTKGVTAACSHNGKSSFYRNLLWLTEKNGLYPNLSKSYVNNKEKEVLVLWGVHHPSNIGDQRAIYHTENAYVSVVSSHYSRRFTPEIAKRPKVRGQEGRINYYWTLLEPGDTIIFEANGNLIAPWYAFALSRGFGSGIITSNASMDECDAKCQTPQGAINSSLPFQNVHPVTIGECPKYVRSTKLRMVTGLRNIPSVQSRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAINGITNKVNSVIEKMNTQFTAVGKEFNKLERRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVKNLYEKVKSQLKNNAKEIGNGCFEFYHKCNNECMESVKNGTYDYPKYSEESKLNREKIDGVKLESMGVYQILAIYSTVASSLVLLVSLGAISFWMCSNGSLQCRICI")

    # Data output and plots
    compartment_data = jOpqua.saveCompartments(output, "examples/flu/new/compartment_flu.csv")
    jOpqua.plotCompartments(compartment_data, ["pop"], "examples/flu/new/compartment_flu.png")

    his_dat = jOpqua.saveHistory(output, "examples/flu/new/history_flu.csv")

    # composition_data = jOpqua.saveComposition(
    #     his_dat, "examples/flu/new/composition_flu.csv",
    #     num_top_sequences=7, track_specific_sequences=[ha_sn89]
    # )
    # jOpqua.plotComposition(
    #     composition_data, "examples/flu/new/composition_flu.png",
    #     # normalized=true, ylabel="Fraction",
    #     normalized=false, ylabel="Number", legend=false
    # )

    # composition_data = jOpqua.saveComposition(
    #     his_dat, "examples/flu/new/flu_responses.csv",
    #     num_top_sequences=7, track_specific_sequences=[ha_sn89],
    #     type_of_composition="Response_imprinted_sequence"
    # )
    # jOpqua.plotComposition(
    #     composition_data, "examples/flu/new/flu_responses.png",
    #     # normalized=true, ylabel="Fraction",
    #     normalized=false, ylabel="Number", legend=false
    # )

    # println("Plot phylogeny")
    # nwks = jOpqua.saveNewick(output, "examples/flu/new/pathogen_newick_flu.nwk")
    # for nwk in nwks
    #     jOpqua.plotPhylogeny(nwk, "examples/flu/new/pathogen_newick_flu.png")
    # end

    println("Calculating distances")
    his_dat = coalesce.(his_dat, "")
    distances = []
    for t in t_vec
        if length(his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]) > 0
            push!(distances, mean([
                sequence_hamming(s, ha_sn89)
                for seqs in his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]
                # for seqs in subset(his_dat, :Time => ByRow(Time -> Time == t))[!, "Pathogen_sequence"]
                for s in split(seqs, jOpqua.WITHIN_HOST_SEPARATOR)
            ]))
        else
            push!(distances, 0.0)
        end
    end
    font_scale = 2.0
    Plots.scalefontsizes(font_scale)
    plot(t_vec, distances, xlabel="Time", ylabel="Mean distance from SN89 (aa)",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(1000, 600), margin=20mm)
    savefig("examples/flu/new/mean_distances.png")

    distances = []
    for t in t_vec
        if length(his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]) > 0
            push!(distances, maximum([
                sequence_hamming(s, ha_sn89)
                for seqs in his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]
                # for seqs in subset(his_dat, :Time => ByRow(Time -> Time == t))[!, "Pathogen_sequence"]
                for s in split(seqs, jOpqua.WITHIN_HOST_SEPARATOR)
            ]))
        else
            push!(distances, 0.0)
        end
    end
    # Plots.scalefontsizes(2.0)
    plot(t_vec, distances, xlabel="Time", ylabel="Max distance from SN89 (aa)",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(1000, 600), margin=20mm)
    savefig("examples/flu/new/max_distances.png")
    Plots.scalefontsizes(1 / font_scale)

    distances = []
    for t in t_vec
        if length(his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]) > 0
            push!(distances, mean([
                crossImmunity(s, ha_sn89)
                for seqs in his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]
                # for seqs in subset(his_dat, :Time => ByRow(Time -> Time == t))[!, "Pathogen_sequence"]
                for s in split(seqs, jOpqua.WITHIN_HOST_SEPARATOR)
            ]))
        else
            push!(distances, 0.0)
        end
    end
    font_scale = 2.0
    Plots.scalefontsizes(font_scale)
    plot(t_vec, distances, xlabel="Time", ylabel="Mean cross-immunity to SN89",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(1000, 600), margin=20mm)
    savefig("examples/flu/new/mean_distances_antigenic.png")

    distances = []
    for t in t_vec
        if length(his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]) > 0
            push!(distances, minimum([
                crossImmunity(s, ha_sn89)
                for seqs in his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]
                # for seqs in subset(his_dat, :Time => ByRow(Time -> Time == t))[!, "Pathogen_sequence"]
                for s in split(seqs, jOpqua.WITHIN_HOST_SEPARATOR)
            ]))
        else
            push!(distances, 0.0)
        end
    end
    # Plots.scalefontsizes(2.0)
    plot(t_vec, distances, xlabel="Time", ylabel="Min cross-immunity to SN89",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(1000, 600), margin=20mm)
    savefig("examples/flu/new/max_distances_antigenic.png")
    Plots.scalefontsizes(1 / font_scale)

    ## HERE WE START GENETIC AND ANTIGENIC MAPS.
    println("Maps")

    # distances = []
    seqs_list = Vector{String}(undef, 0)
    times_dict = Dict()
    times_first = Vector{Float64}(undef, 0)
    occurrences = Dict()

    for t in t_vec#[end]
        if length(his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]) > 0
            for seqs in his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]
                for s in split(seqs, jOpqua.WITHIN_HOST_SEPARATOR)
                    if !haskey(times_dict, s)
                        times_dict[s] = [t]
                        push!(seqs_list, s)
                        push!(times_first, t)
                        occurrences[s] = 1
                    else
                        occurrences[s] += 1
                    end
                end
            end
        end
    end

    println("Obtained sequences, " * string(length(seqs_list)))

    # ancestors = jOpqua.ancestors(output, "pop", seqs_list)
    # seqs_list = ancestors.sequence
    # times_first = ancestors.time

    println("Obtained ancestors, " * string(length(seqs_list)))

    correction_number = 0.0000001
    dist_mat_ant = zeros((length(seqs_list), length(seqs_list)))
    dist_mat_gen = zeros((length(seqs_list), length(seqs_list)))
    Threads.@threads for i in 1:length(seqs_list)
        for j in i:length(seqs_list)
            dist_mat_ant[i, j] = (1.0 / (crossImmunity(seqs_list[i], seqs_list[j]) + correction_number)) - (1.0 / (1.0 + correction_number))
            dist_mat_gen[i, j] = sequence_hamming(seqs_list[i], seqs_list[j])
            dist_mat_ant[j, i] = dist_mat_ant[i, j]
            dist_mat_gen[j, i] = dist_mat_gen[i, j]
        end
    end

    println("Calculated pairwise distances")

    pca_ant = fit(PCA, dist_mat_ant; maxoutdim=2)
    pca_gen = fit(PCA, dist_mat_gen; maxoutdim=2)
    proj_ant = projection(pca_ant)
    proj_gen = projection(pca_gen)

    println("Antigenic PCA variance fractions: " * string(principalvars(pca_ant) / var(pca_ant)))
    println("Genetic PCA variance fractions: " * string(principalvars(pca_gen) / var(pca_gen)))

    # proj_ant = umap(dist_mat_ant, 2; metric=:precomputed)'
    # proj_gen = umap(dist_mat_gen, 2; metric=:precomputed)'

    occurrences_list = Vector{Int64}(undef, 0)
    for s in seqs_list
        if haskey(occurrences, s)
            push!(occurrences_list, occurrences[s])
        else
            push!(occurrences_list, 0)
        end
    end

    # occurrences_list_sizes = 10.0 * occurrences_list / maximum(occurrences_list)
    occurrences_list_sizes = 20.0 * min.(2.0, (occurrences_list .+ 0) / 10)
    # occurrences_list_sizes = 20.0 * (occurrences_list .+ 0) / 10

    # ant = DataFrame(proj_ant, ["Component 1", "Component 2"])
    # ant[:, "Time"] = times_first
    # ant[:, "Occurrences"] = occurrences_list
    # ant[:, "Occurrences_mod"] = occurrences_list_sizes
    # ant[:, "Alpha"] = ones(length(seqs_list))
    # ant[:, "Linewidth"] = ones(length(seqs_list))
    # ant = sort(ant, ("Occurrences_mod"), rev=true)
    # ant[(ant[:, "Occurrences"].==0), "Alpha"] .= 0.4
    # ant[(ant[:, "Occurrences"].==0), "Occurrences_mod"] .= 4
    # ant[(ant[:, "Occurrences"].==0), "Linewidth"] .= 0

    # font_scale = 2.0
    # Plots.scalefontsizes(font_scale)
    # scatter(
    #     ant[:, "Component 1"], ant[:, "Component 2"], zcolor=ant[:, "Time"], ms=ant[:, "Occurrences_mod"], c=:viridis,
    #     xlabel="Component 1", ylabel="Component 2",
    #     alpha=ant[:, "Alpha"],
    #     legend=false, colorbar=true,
    #     markerstrokewidth=ant[:, "Linewidth"],
    #     linewidth=3.0, thickness_scaling=1.0,
    #     grid=false,
    #     colorbar_title="Time",
    # )
    # plot!(size=(1000, 800), margin=20mm)
    # savefig("examples/flu/new/pca_antigenic.png")
    # Plots.scalefontsizes(1 / font_scale)

    # CSV.write("examples/flu/new/components_antigenic.csv", ant)
    # CSV.write("examples/flu/new/distances_antigenic.csv", DataFrame(dist_mat_ant, :auto))

    gen = DataFrame(proj_gen, ["Component 1", "Component 2"])
    gen[:, "Time"] = times_first
    gen[:, "Occurrences"] = occurrences_list
    gen[:, "Occurrences_mod"] = occurrences_list_sizes
    gen[:, "Alpha"] = ones(length(seqs_list))
    gen[:, "Linewidth"] = ones(length(seqs_list))
    gen = sort(gen, ("Occurrences_mod"), rev=true)

    gen[(gen[:, "Occurrences"].==0), "Alpha"] .= 0.4
    gen[(gen[:, "Occurrences"].==0), "Occurrences_mod"] .= 4
    gen[(gen[:, "Occurrences"].==0), "Linewidth"] .= 0

    font_scale = 2.0
    Plots.scalefontsizes(font_scale)
    scatter(
        gen[:, "Component 1"], gen[:, "Component 2"], zcolor=gen[:, "Time"], ms=gen[:, "Occurrences_mod"], c=:viridis,
        xlabel="Component 1", ylabel="Component 2",
        alpha=gen[:, "Alpha"],
        legend=false, colorbar=true,
        markerstrokewidth=gen[:, "Linewidth"],
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,
        colorbar_title="Time",
    )
    plot!(size=(1000, 800), margin=20mm)
    savefig("examples/flu/new/pca_genetic.png")
    Plots.scalefontsizes(1 / font_scale)

    CSV.write("examples/flu/new/components_genetic.csv", gen)
    CSV.write("examples/flu/new/distances_genetic.csv", DataFrame(dist_mat_gen, :auto))

    CSV.write("examples/flu/new/sequences.csv", DataFrame(Sequences=seqs_list))
end

# Model setup
function run(seed::Int64, t_vec::Vector{Float64})


    # Parameters
    start_genome = ha_sn89
    optimal_genome = ha_sn89#"BBBB"
    genome_length = length(optimal_genome)

    pat_type = jOpqua.newPathogenType(
        "pat_type",
        num_loci=genome_length,
        possible_alleles="AB",
        contactSpecificCoefficient=s::String -> 1.0 + (0.1 * (genome_length - hamming(s, optimal_genome)) / genome_length),
        receiveContactHostwideCoefficient=s::String -> 0.0, # makes infected hosts immune to superinfection
    )

    res_type = jOpqua.newResponseType(
        "res_type",
        reactivityCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> crossImmunity(imp_g, pat_g),
        # transmissionEfficiencyInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 0.0 : 1.0,
        clearanceInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> crossImmunity(imp_g, pat_g)*0.75 + 1.0,
        # responseAcquisitionInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 0.0 : 1.0,
    )

    hos_type = jOpqua.newHostType("hos_type")

    pop_type = jOpqua.newPopulationType(
        "pop_type",
        clearance_coefficient=1.0,
        contact_coefficient=1.05,
        # response_acquisition_coefficient=0.5,
        response_acquisition_upon_clearance_coefficient=0.2,
        response_loss_coefficient=0.01,
        receive_contact_coefficient=1.0,
        mutations_per_generation_coefficient=5e-4,#0.0005,
        inoculum_coefficient=1.0,
        # death_coefficient=0.001,
        # birth_coefficient=0.001,
        pathogenFractions=jOpqua.pathogenFractionsProportionalFitness,
        response_types=Dict{String,jOpqua.ResponseType}([(res_type.id => res_type)]),
        developResponses=(
            pathogen::jOpqua.Pathogen, host::jOpqua.Host,
            existing_responses::Dict{Tuple{String,String,String,String},jOpqua.Response},
            response_types::Dict{String,jOpqua.ResponseType},
            birth_time::Float64
        ) -> [jOpqua.deNovoResponse(
            pathogen, host, existing_responses, response_types, birth_time;
            response_type_id="res_type"
        )],
    )

    num_hosts = 10000
    num_infected = Int(num_hosts * 0.05)
    host_genome = ""

    # Setup
    model = jOpqua.newModel()
    pop = jOpqua.newPopulation!("pop", pop_type, model)
    jOpqua.addHostsToPopulation!(num_hosts, host_genome, hos_type, pop, model)
    pat = jOpqua.newPathogen!(start_genome, pop, pat_type)

    for h in 1:num_infected
        jOpqua.addPathogenToHost!(pat, h, pop, model)
    end

    # Simulate
    Random.seed!(seed)
    model, output = jOpqua.simulate!(
        model, t_vec, population_host_samples=Dict("pop" => 200)
    )

    # Data output and plots
    compartment_data = jOpqua.saveCompartments(output, "examples/flu/new/compartment_flu.csv")
    jOpqua.plotCompartments(compartment_data, ["pop"], "examples/flu/new/compartment_flu.png")

    his_dat = jOpqua.saveHistory(output, "examples/flu/new/history_flu.csv")
    composition_data = jOpqua.saveComposition(
        his_dat, "examples/flu/new/composition_flu.csv",
        num_top_sequences=7, #track_specific_sequences=["AAAA", "BBBB"]legend=false
    )
    jOpqua.plotComposition(
        composition_data, "examples/flu/new/composition_flu.png",
        normalized=true, ylabel="Fraction", legend=false
    )

    composition_data = jOpqua.saveComposition(
        his_dat, "examples/flu/new/composition_flu_responses.csv",
        num_top_sequences=7, #track_specific_sequences=["AAAA", "BBBB"],
        type_of_composition="Response_imprinted_sequence"
    )
    jOpqua.plotComposition(
        composition_data, "examples/flu/new/composition_flu_responses.png",
        # normalized=true, ylabel="Fraction",
        normalized=false, ylabel="Number", legend=false
    )

    # nwks = jOpqua.saveNewick(output, "examples/flu/new/pathogen_newick_flu.nwk", branch_length="Time", info_separator=" ")
    # for nwk in nwks
    #     jOpqua.plotPhylogeny(nwk, "examples/flu/new/pathogen_newick_flu.png")
    # end
    #

    return output
end



# run(1, collect(0.0:2.0:4.0)) # compile
# @time run(10, collect(0.0:2.0:1500.0))
# Total events: 3342817
#   7.751170 seconds (82.72 M allocations: 14.894 GiB, 20.86% gc time, 7.99% compilation time: <1% of which was recompilation)
# 25 May 2026 Julia 1.12.6 Apple M3 Max 128 GB RAM

# @profview run(2, collect(0.0:2.0:1500.0))

max_time = 1500
t_steps = 500
println("Num threads: " * string(nthreads()))

@time out::jOpqua.Output = run(10, collect(0.0:(max_time/t_steps):max_time))
analyze(out, collect(0.0:(max_time/t_steps):max_time))
