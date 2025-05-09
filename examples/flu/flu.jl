# Run from base jOpqua directory as
# julia --project=. examples/flu/acquired_immunity.jl
# unless viewing flamegraph, then run from console

using Revise
using jOpqua

using StaticArrays
using Random

using Distances

using Statistics
using DataFrames
using Plots

using BenchmarkTools
using ProfileView

# Model setup
function run(seed::Int64, t_vec::Vector{Float64})
    # Parameters
    ha_sn89 = "MKAKLLVLLCAFTATDADTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCRLKGIAPLQLGNCSIAGWILGNPECESLFSKESWSYIAETPNSENGTCYPGYFADYEELREQLSSVSSFERFEIFPKESSWPNHTVTKGVTAACSHNGKSSFYRNLLWLTEKNGLYPNLSKSYVNNKEKEVLVLWGVHHPSNIGDQRAIYHTENAYVSVVSSHYSRRFTPEIAKRPKVRGQEGRINYYWTLLEPGDTIIFEANGNLIAPWYAFALSRGFGSGIITSNASMDECDAKCQTPQGAINSSLPFQNVHPVTIGECPKYVRSTKLRMVTGLRNIPSVQSRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAINGITNKVNSVIEKMNTQFTAVGKEFNKLERRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVKNLYEKVKSQLKNNAKEIGNGCFEFYHKCNNECMESVKNGTYDYPKYSEESKLNREKIDGVKLESMGVYQILAIYSTVASSLVLLVSLGAISFWMCSNGSLQCRICI*"
    rbs_residues = collect(range(140, 210))
    rbs_epitope_key_residues = [144, 147, 159, 169, 170, 171, 172, 173, 203, 207]

    function crossImmunity(
        seq1, seq2;
        epitope_residues=rbs_epitope_key_residues,
        distance_half_all_residues=length(ha_sn89) * 0.000025, hill_coef_all_residues=5.0,
        distance_half_epitope_residues=length(rbs_epitope_key_residues) * 0.25, hill_coef_epitope_residues=5.0)

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

    function proteinFitness(
        seq1, seq2;
        functional_site=rbs_residues,
        distance_half_all_residues=length(ha_sn89) * 0.5, hill_coef_all_residues=5.0,
        distance_half_functional_site_residues=length(rbs_residues) * 0.5, hill_coef_functional_site_residues=5.0)

        distance_key_residues = 0.0
        for p in functional_site
            distance_key_residues += seq1[p] != seq2[p]
        end

        return (1.0 - jOpqua.hillFunction(
            Float64(hamming(seq1, seq2)),
            distance_half_all_residues, hill_coef_all_residues
        )) * (1.0 - jOpqua.hillFunction(
            distance_key_residues,
            distance_half_functional_site_residues, hill_coef_functional_site_residues
        ))
        # return 1.0
    end

    pat_type = jOpqua.newPathogenType(
        "pat_type",
        num_loci=length(ha_sn89),
        possible_alleles="ARNDCEQGHILKMFPSTWYV*",
        # mean_effective_inoculum=1.0,
        # mean_mutations_per_replication=0.00003,
        contactCoefficient=s::String -> proteinFitness(s, ha_sn89),
        receiveContactCoefficient=s::String -> 0.0,
    )

    res_type = jOpqua.newResponseType(
        "res_type",
        reactivityCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g, #crossImmunity(imp_g, pat_g),
        transmissionEfficiencySpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> 1.0 - crossImmunity(imp_g, pat_g),
        clearanceSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> 1.0e3 * crossImmunity(imp_g, pat_g),
        # responseAcquisitionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 0.0 : 1.0,
    )

    hos_type = jOpqua.newHostType("hos_type")

    pop_type = jOpqua.newPopulationType(
        "pop_type",
        clearance_coefficient=1.0e-3,
        contact_coefficient=1.05,
        response_acquisition_coefficient=0.5,
        response_loss_coefficient=0.0,
        receive_contact_coefficient=1.0,
        mutations_upon_infection_coefficient=0.1,
        inoculum_coefficient=1.0,
        pathogenFractions=jOpqua.pathogenFractionsProportionalFitness,
        response_types=Dict{String,jOpqua.ResponseType}([(res_type.id => res_type)]),
        developResponses=(
            pathogen::jOpqua.Pathogen, host::jOpqua.Host,
            existing_responses::Dict{Tuple{String,String,String,String},jOpqua.Response},
            response_types::Dict{String,jOpqua.ResponseType},
            birth_time::Float64
        ) -> jOpqua.deNovoResponse(
            pathogen, host, existing_responses, response_types, birth_time;
            response_type_id="res_type"
        ),
    )

    num_hosts = 5000
    num_infected = Int(num_hosts * 0.05)
    host_genome = ""

    # Setup
    model = jOpqua.newModel()
    pop = jOpqua.newPopulation!("pop", pop_type, model)
    jOpqua.addHostsToPopulation!(num_hosts, host_genome, hos_type, pop, model)
    pat = jOpqua.newPathogen!(ha_sn89, pop, pat_type)

    for h in 1:num_infected
        jOpqua.addPathogenToHost!(pat, h, pop, model)
    end

    # Simulate
    Random.seed!(seed)
    model, output = jOpqua.simulate!(
        model, t_vec, population_host_samples=Dict("pop" => 200)
    )

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
        end
    end
    font_scale = 2.0
    Plots.scalefontsizes(font_scale)
    plot(distances, xlabel="Time", ylabel="Distance from SN98 (aa)",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(800, 700))
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
        end
    end
    # Plots.scalefontsizes(2.0)
    plot(distances, xlabel="Time", ylabel="Distance from SN98 (aa)",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(800, 700))
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
        end
    end
    font_scale = 2.0
    Plots.scalefontsizes(font_scale)
    plot(distances, xlabel="Time", ylabel="Cross-immunity to SN98",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(800, 700))
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
        end
    end
    # Plots.scalefontsizes(2.0)
    plot(distances, xlabel="Time", ylabel="Cross-immunity to SN98",
        legend=false,
        linewidth=3.0, thickness_scaling=1.0,
        grid=false,)
    plot!(size=(800, 700))
    savefig("examples/flu/max_distances_antigenic.png")
    Plots.scalefontsizes(1 / font_scale)

end

run(1, collect(0.0:2.0:4.0)) # compile
@time run(0, collect(0.0:2:100.0))
