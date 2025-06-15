# Run from base jOpqua directory as
# julia --project=. --threads=auto examples/flu/flu.jl
# unless viewing flamegraph, then run from console

using Revise
using jOpqua

using StaticArrays
using Random

using Distances

using Statistics
using MultivariateStats
# using UMAP
using DataFrames
using Plots
using Plots.PlotMeasures
using CSV

using Base.Threads

using BenchmarkTools
using ProfileView

# Model setup
function run(seed::Int64, t_vec::Vector{Float64})
    # Parameters
    ha_sn89 = "MKAKLLVLLCAFTATDADTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCRLKGIAPLQLGNCSIAGWILGNPECESLFSKESWSYIAETPNSENGTCYPGYFADYEELREQLSSVSSFERFEIFPKESSWPNHTVTKGVTAACSHNGKSSFYRNLLWLTEKNGLYPNLSKSYVNNKEKEVLVLWGVHHPSNIGDQRAIYHTENAYVSVVSSHYSRRFTPEIAKRPKVRGQEGRINYYWTLLEPGDTIIFEANGNLIAPWYAFALSRGFGSGIITSNASMDECDAKCQTPQGAINSSLPFQNVHPVTIGECPKYVRSTKLRMVTGLRNIPSVQSRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAINGITNKVNSVIEKMNTQFTAVGKEFNKLERRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVKNLYEKVKSQLKNNAKEIGNGCFEFYHKCNNECMESVKNGTYDYPKYSEESKLNREKIDGVKLESMGVYQILAIYSTVASSLVLLVSLGAISFWMCSNGSLQCRICI"
    rbs_residues = collect(range(140, 210))
    rbs_epitope_key_residues = [144, 147, 159, 169, 170, 171, 172, 173, 203, 207]

    ha_sn89_mut = "MKAKLLVLLCAFTATDADTICIGYHANTSTDTVDTVLEKNVTVTHNVNLLEDSHNGKLCRLKGIAPLQLGNCSIAGWILGNPECESLFSKESWSYIAETPNSENATCEYGTCADYEELREQLSSVSSFERFEPFPNESSWPNHLVTPGVTAACSHNGKASFYRNLLWLTEKNGLYPNLPKSYVNNKEKEVLVLWGVHHPDNIGDQRAIYHTENAYVSVVSSHYSRRFTPEIAKRPKVRGQEMEINYYWTLLEPGDTIIFEANGNLIAPWDAFALSRGFGSGIITSNASMDECDAKCQTPQGAINSSLPFQNVHPVTIGECPKYGRSTKLRMVTGLRNIPSDQSRGLFGAIAGFCEGGWTEMIDGWYGYHHQNEQGSGYAADQKSTQNGITGITNKVNSVIEKMNTQFTAQGKEFNKLERRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVKNCYEKVKSQLKNNAKEIGNGCFEFYHKCNNECMESVKNGTYDYPKYSEESKLNREKIDGVKLESMTVYQILAIYSTVASSLVLLVSAGAISFWMCSNGSLKCRICI"

    # ha_sn89 = "AAAAAAAAAA"
    # ha_sn89_mut = "*AAAAAAAAA"
    # rbs_residues = collect(range(1,10))
    # rbs_epitope_key_residues = [1,2,3,4]

    function crossImmunity(
        seq1, seq2;
        epitope_residues=rbs_epitope_key_residues,
        distance_half_all_residues=length(ha_sn89) * 0.1, hill_coef_all_residues=3.0,
        distance_half_epitope_residues=length(rbs_epitope_key_residues) * 0.1, hill_coef_epitope_residues=2.0)

        distance_key_residues = 0.0
        for p in epitope_residues
            distance_key_residues += seq1[p] != seq2[p]
        end

        return (1.0 - jOpqua.hillFunction(
            Float64(hamming(seq1, seq2)),
            distance_half_all_residues, hill_coef_all_residues
        )) * (1.0 - jOpqua.hillFunction(
            distance_key_residues,
            distance_half_epitope_residues, hill_coef_epitope_residues
        ))
        # return seq1 == seq2
    end

    println(("TEST SAME: ", 1.0 - crossImmunity(ha_sn89, ha_sn89)))
    println(("TEST MUT: ", 1.0 - crossImmunity(ha_sn89, ha_sn89_mut)))

    function proteinFitness(
        seq, seq_wt;
        functional_site=rbs_residues,
        distance_half_all_residues=length(ha_sn89) * 0.5, hill_coef_all_residues=3.0,
        distance_half_functional_site_residues=length(rbs_residues) * 0.5, hill_coef_functional_site_residues=3.0)

        if occursin("*", seq)
            return 0.0
        else
            distance_key_residues = 0.0
            for p in functional_site
                distance_key_residues += seq[p] != seq_wt[p]
            end

            # return (1.0 - jOpqua.hillFunction(
            #     Float64(hamming(seq, seq_wt)),
            #     distance_half_all_residues, hill_coef_all_residues
            # )) * (1.0 - jOpqua.hillFunction(
            #     distance_key_residues,
            #     distance_half_functional_site_residues, hill_coef_functional_site_residues
            # ))
            return 1.0
        end
    end

    println(("TEST SAME fit: ", proteinFitness(ha_sn89, ha_sn89)))
    println(("TEST MUT fit: ", proteinFitness(ha_sn89, ha_sn89_mut)))

    println(("TEST SAME hamm: ", hamming(ha_sn89, ha_sn89)))
    println(("TEST MUT hamm: ", hamming(ha_sn89, ha_sn89_mut)))

    max_immunity = 1.0

    pat_type = jOpqua.newPathogenType(
        "pat_type",
        num_loci=length(ha_sn89),
        possible_alleles="ARNDCEQGHILKMFPSTWYV*",
        # mean_effective_inoculum=1.0,
        # mean_mutations_per_replication=0.00003,
        contactSpecificCoefficient=s::String -> proteinFitness(s, ha_sn89),
        receiveContactHostwideCoefficient=s::String -> 0.0,
    )

    res_type_spe = jOpqua.newResponseType(
        "Specific",
        reactivityCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> crossImmunity(imp_g, pat_g),
        transmissionEfficiencyInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> (1.0 - max_immunity * crossImmunity(imp_g, pat_g)),
        clearanceInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> 100.0 * 7.0e5 * crossImmunity(imp_g, pat_g) + 1.0,
        responseLossStaticSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String) -> 1.0,
        # responseAcquisitionInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 0.0 : 1.0,
    )

    res_type_bro = jOpqua.newResponseType(
        "Broad_temp",
        reactivityCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> 1.0,
        transmissionEfficiencyInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> 1.0 - max_immunity,
        clearanceInteractionHostwideCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> 100.0 * 7.0e5,
        responseLossStaticSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String) -> 100.0, # total immunity for ~4 months
        # responseAcquisitionInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 0.0 : 1.0,
    )

    function developResponse(
        pathogen::jOpqua.Pathogen, host::jOpqua.Host,
        existing_responses::Dict{Tuple{String,String,String,String},jOpqua.Response},
        response_types::Dict{String,jOpqua.ResponseType},
        birth_time::Float64)

        for response in host.responses
            if response.type.id == "Broad_temp"
                return [
                    jOpqua.deNovoResponse(
                        pathogen, host, existing_responses, response_types, birth_time;
                        response_type_id="Specific"
                    )
                ]
            end
        end
        return [
            jOpqua.deNovoResponse(
                pathogen, host, existing_responses, response_types, birth_time;
                response_type_id="Specific"
            ),
            jOpqua.deNovoResponse(
                pathogen, host, existing_responses, response_types, birth_time;
                response_type_id="Broad_temp"
            )
        ]
    end

    hos_type = jOpqua.newHostType("hos_type")

    pop_type = jOpqua.newPopulationType(
        "pop_type",
        clearance_coefficient=7.0e-5, # 2 weeks
        contact_coefficient=0.125 * 5.0, # R_0 of ~5.0
        response_acquisition_coefficient=0.125, # ~8 days infectiousness
        response_loss_coefficient=0.0,#14e-3 / 365, # 3.3e-5 birth rate
        receive_contact_coefficient=1.0,
        mutations_upon_infection_coefficient=0.161 * 1.0, # 0.161 = 566 aa * 3 nt/codon * (1-1/(21 mut aa - 1 WT)) * (1-(1-(1/100000 mut per site per replication))^(10 rounds of replication before transmission) )
        inoculum_coefficient=1.0,
        pathogenFractions=jOpqua.pathogenFractionsProportionalFitness,
        response_types=Dict{String,jOpqua.ResponseType}([(res_type_spe.id => res_type_spe), (res_type_bro.id => res_type_bro)]),
        developResponses=developResponse
    )

    num_hosts = 100000
    num_infected = 5 #Int(num_hosts * 0.01)
    host_genome = ""

    println("Creating model")

    # Setup
    model = jOpqua.newModel()
    pop = jOpqua.newPopulation!("pop", pop_type, model)
    println("Adding hosts")
    jOpqua.addHostsToPopulation!(num_hosts, host_genome, hos_type, pop, model)
    pat = jOpqua.newPathogen!(ha_sn89, pop, pat_type)
    pat_mut = jOpqua.newPathogen!(ha_sn89_mut, pop, pat_type)

    println("Adding pathogens")

    for h in 1:num_infected
        # println((h,num_hosts))
        jOpqua.addPathogenToHost!(pat, h, pop, model)
        # jOpqua.addPathogenToHost!(pat_mut, h, pop, model)
        # jOpqua.acquireResponse!(1, h, 1, model, rand())
    end

    println("Setup complete")

    # for h in 1:num_infected
    #     jOpqua.addPathogenToHost!(pat_mut, h, pop, model)
    # end

    # println(("Host responses: ", [h.responses for h in pop.hosts]))

    # Simulate
    Random.seed!(seed)
    model, output = jOpqua.simulate!(
        model, t_vec, population_host_samples=Dict("pop" => 200)
    )
    # println(("Host responses: ", [h.responses for h in pop.hosts]))

    # Data output and plots
    compartment_data = jOpqua.saveCompartments(output, "examples/flu/compartment_flu.csv")
    jOpqua.plotCompartments(compartment_data, ["pop"], "examples/flu/compartment_flu.png")

    his_dat = jOpqua.saveHistory(output, "examples/flu/history_flu.csv")
    composition_data = jOpqua.saveComposition(
        his_dat, "examples/flu/composition_flu.csv",
        num_top_sequences=7, track_specific_sequences=[ha_sn89]
    )
    jOpqua.plotComposition(
        composition_data, "examples/flu/composition_flu.png",
        # normalized=true, ylabel="Fraction",
        normalized=false, ylabel="Number", legend=false
    )

    composition_data = jOpqua.saveComposition(
        his_dat, "examples/flu/flu_responses.csv",
        num_top_sequences=7, track_specific_sequences=[ha_sn89],
        type_of_composition="Response_imprinted_sequence"
    )
    jOpqua.plotComposition(
        composition_data, "examples/flu/flu_responses.png",
        # normalized=true, ylabel="Fraction",
        normalized=false, ylabel="Number", legend=false
    )

    # nwks = jOpqua.saveNewick(output, "examples/flu/pathogen_newick_flu.nwk")
    # for nwk in nwks
    #     jOpqua.plotPhylogeny(nwk, "examples/flu/pathogen_newick_flu.png")
    # end


    distances = []
    for t in t_vec
        if length(his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]) > 0
            push!(distances, mean([
                hamming(s, ha_sn89)
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
    plot(t_vec, distances, xlabel="Time", ylabel="Mean distance from SN98 (aa)",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(1000, 600), margin=20mm)
    savefig("examples/flu/mean_distances.png")

    distances = []
    for t in t_vec
        if length(his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]) > 0
            push!(distances, maximum([
                hamming(s, ha_sn89)
                for seqs in his_dat[(his_dat.Time.==t).&(his_dat.Pathogen_sequence.!=""), "Pathogen_sequence"]
                # for seqs in subset(his_dat, :Time => ByRow(Time -> Time == t))[!, "Pathogen_sequence"]
                for s in split(seqs, jOpqua.WITHIN_HOST_SEPARATOR)
            ]))
        else
            push!(distances, 0.0)
        end
    end
    # Plots.scalefontsizes(2.0)
    plot(t_vec, distances, xlabel="Time", ylabel="Max distance from SN98 (aa)",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(1000, 600), margin=20mm)
    savefig("examples/flu/max_distances.png")
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
    plot(t_vec, distances, xlabel="Time", ylabel="Mean cross-immunity to SN98",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(1000, 600), margin=20mm)
    savefig("examples/flu/mean_distances_antigenic.png")

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
    plot(t_vec, distances, xlabel="Time", ylabel="Max cross-immunity to SN98",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(1000, 600), margin=20mm)
    savefig("examples/flu/max_distances_antigenic.png")
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

    ancestors = jOpqua.ancestors(output, "pop", seqs_list)
    seqs_list = ancestors.sequence
    times_first = ancestors.time

    println("Obtained ancestors, " * string(length(seqs_list)))

    correction_number = 0.0000001
    dist_mat_ant = zeros((length(seqs_list), length(seqs_list)))
    dist_mat_gen = zeros((length(seqs_list), length(seqs_list)))
    Threads.@threads for i in 1:length(seqs_list)
        for j in i:length(seqs_list)
            dist_mat_ant[i, j] = (1.0 / (crossImmunity(seqs_list[i], seqs_list[j]) + correction_number)) - (1.0 / (1.0 + correction_number))
            dist_mat_gen[i, j] = hamming(seqs_list[i], seqs_list[j])
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

    ant = DataFrame(proj_ant, ["Component 1", "Component 2"])
    ant[:, "Time"] = times_first
    ant[:, "Occurrences"] = occurrences_list
    ant[:, "Occurrences_mod"] = occurrences_list_sizes
    ant[:, "Alpha"] = ones(length(seqs_list))
    ant[:, "Linewidth"] = ones(length(seqs_list))
    ant = sort(ant, ("Occurrences_mod"), rev=true)

    gen = DataFrame(proj_gen, ["Component 1", "Component 2"])
    gen[:, "Time"] = times_first
    gen[:, "Occurrences"] = occurrences_list
    gen[:, "Occurrences_mod"] = occurrences_list_sizes
    gen[:, "Alpha"] = ones(length(seqs_list))
    gen[:, "Linewidth"] = ones(length(seqs_list))
    gen = sort(gen, ("Occurrences_mod"), rev=true)

    ant[(ant[:, "Occurrences"].==0), "Alpha"] .= 0.4
    gen[(gen[:, "Occurrences"].==0), "Alpha"] .= 0.4
    ant[(ant[:, "Occurrences"].==0), "Occurrences_mod"] .= 4
    gen[(gen[:, "Occurrences"].==0), "Occurrences_mod"] .= 4
    ant[(ant[:, "Occurrences"].==0), "Linewidth"] .= 0
    gen[(gen[:, "Occurrences"].==0), "Linewidth"] .= 0

    font_scale = 2.0
    Plots.scalefontsizes(font_scale)
    scatter(
        ant[:, "Component 1"], ant[:, "Component 2"], zcolor=ant[:, "Time"], ms=ant[:, "Occurrences_mod"], c=:viridis,
        xlabel="Component 1", ylabel="Component 2",
        alpha=ant[:, "Alpha"],
        legend=false, colorbar=true,
        markerstrokewidth=ant[:, "Linewidth"],
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,
        colorbar_title="Time",
    )
    plot!(size=(1000, 800), margin=20mm)
    savefig("examples/flu/pca_antigenic.png")

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
    savefig("examples/flu/pca_genetic.png")
    Plots.scalefontsizes(1 / font_scale)

    CSV.write("examples/flu/components_antigenic.csv", ant)
    CSV.write("examples/flu/components_genetic.csv", gen)
    # CSV.write("examples/flu/distances_antigenic.csv", DataFrame(dist_mat_ant, :auto))
    # CSV.write("examples/flu/distances_genetic.csv", DataFrame(dist_mat_gen, :auto))
    CSV.write("examples/flu/sequences.csv", DataFrame(Sequences=seqs_list))
end

println("Num threads: " * string(nthreads()))

# run(1, collect(0.0:2.0:4.0)) # compile
@time run(1, collect(0.0:2.0:365.0))
