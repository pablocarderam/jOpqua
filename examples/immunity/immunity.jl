# Run from base jOpqua directory as
# julia --project=. examples/immunity/immunity.jl
# unless viewing flamegraph, then run from console

using Revise
using jOpqua

using StaticArrays
using Random

using Distances

using BenchmarkTools
using ProfileView

# Model setup
function run(seed::Int64)
    # Parameters
    start_genome = "AAAA"
    optimal_genome = "BBBB"

    pat_type = jOpqua.newPathogenType(
        "pat_type",
        num_loci=4,
        possible_alleles="AB",
        mean_effective_inoculum=1.0,
        mean_mutations_per_replication=0.001,
        contactCoefficient=s::String -> 1.0 + (0.1 * (4.0 - hamming(s, optimal_genome)) / 4.0),
        receiveContactCoefficient=s::String -> 0.0,
    )

    pop_type = jOpqua.newPopulationType(
        "pop_type",
        clearance_coefficient=1.0,
        contact_coefficient=1.05,
        receive_contact_coefficient=1.0,
        pathogenFractions=jOpqua.pathogenFractionsProportionalFitness,
    )

    num_hosts = 10000
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

    t_vec = collect(0.0:2.0:1500.0)

    # Simulate
    Random.seed!(seed)
    model, output = jOpqua.simulate!(
        model, t_vec, population_host_samples=Dict("pop" => 200)
    )

    # Data output and plots
    compartment_data = jOpqua.saveCompartments(output, "examples/immunity/compartment_immunity.csv")
    jOpqua.plotCompartments(compartment_data, ["pop"], "examples/immunity/compartment_immunity.png")

    his_dat = jOpqua.saveHistory(output, "examples/immunity/history_immunity.csv")
    composition_data = jOpqua.saveComposition(
        his_dat, "examples/immunity/composition_immunity.csv",
        num_top_sequences=7, track_specific_sequences=["AAAA", "BBBB"]
    )
    jOpqua.plotComposition(
        composition_data, "examples/immunity/composition_immunity.png",
        normalized=true, ylabel="Fraction",
    )

    nwks = jOpqua.saveNewick(output, "examples/immunity/pathogen_newick_immunity.nwk")
    for nwk in nwks
        jOpqua.plotPhylogeny(nwk, "examples/immunity/pathogen_newick_immunity.png")
    end
end

run(0)
