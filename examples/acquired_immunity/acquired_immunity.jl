# Run from base jOpqua directory as
# julia --project=. examples/acquired_immunity/acquired_immunity.jl
# unless viewing flamegraph, then run from console

using Revise
using jOpqua

using StaticArrays
using Random

using Distances

using BenchmarkTools
using ProfileView

# Model setup
function run(seed::Int64, t_vec::Vector{Float64})
    # Parameters
    start_genome = "AAAA"
    optimal_genome = "BBBB"

    pat_type = jOpqua.newPathogenType(
        "pat_type",
        num_loci=4,
        possible_alleles="AB",
        mean_effective_inoculum=1.0,
        mean_mutations_per_replication=0.00003,
        contactCoefficient=s::String -> 1.0 + (0.1 * (4.0 - hamming(s, optimal_genome)) / 4.0),
        receiveContactCoefficient=s::String -> 0.0,
    )

    res_type = jOpqua.newResponseType(
        "res_type",
        reactivityCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 1.0 : 0.0,
        infectionCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 0.0 : 1.0,
        clearanceSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 1.0e3 : 1.0,
        # responseAcquisitionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 0.0 : 1.0,
    )

    pop_type = jOpqua.newPopulationType(
        "pop_type",
        clearance_coefficient=1.0e-3,
        contact_coefficient=1.05,
        response_acquisition_coefficient=0.5,
        response_loss_coefficient=0.0,
        receive_contact_coefficient=1.0,
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

    num_hosts = 50000
    num_infected = Int(num_hosts * 0.05)
    host_genome = ""

    # Setup
    model = jOpqua.newModel()
    pop = jOpqua.newPopulation!("pop", pop_type, model)
    jOpqua.addHostsToPopulation!(num_hosts, host_genome, pop, model)
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
    compartment_data = jOpqua.saveCompartments(output, "examples/acquired_immunity/compartment_acquired_immunity.csv")
    jOpqua.plotCompartments(compartment_data, ["pop"], "examples/acquired_immunity/compartment_acquired_immunity.png")

    his_dat = jOpqua.saveHistory(output, "examples/acquired_immunity/history_acquired_immunity.csv")
    composition_data = jOpqua.saveComposition(
        his_dat, "examples/acquired_immunity/composition_acquired_immunity.csv",
        num_top_sequences=7, track_specific_sequences=["AAAA", "BBBB"]
    )
    jOpqua.plotComposition(
        composition_data, "examples/acquired_immunity/composition_acquired_immunity.png",
        # normalized=true, ylabel="Fraction",
        normalized=false, ylabel="Number",
    )

    composition_data = jOpqua.saveComposition(
        his_dat, "examples/acquired_immunity/composition_acquired_immunity_responses.csv",
        num_top_sequences=7, track_specific_sequences=["AAAA", "BBBB"],
        type_of_composition="Response_imprinted_sequence"
    )
    jOpqua.plotComposition(
        composition_data, "examples/acquired_immunity/composition_acquired_immunity_responses.png",
        # normalized=true, ylabel="Fraction",
        normalized=false, ylabel="Number",
    )

    nwks = jOpqua.saveNewick(output, "examples/acquired_immunity/pathogen_newick_acquired_immunity.nwk")
    for nwk in nwks
        jOpqua.plotPhylogeny(nwk, "examples/acquired_immunity/pathogen_newick_acquired_immunity.png")
    end
end

run(1, collect(0.0:2.0:4.0)) # compile
@time run(0, collect(0.0:0.20:100.0))
