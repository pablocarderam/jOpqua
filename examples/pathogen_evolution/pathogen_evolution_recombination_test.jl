# Run from base jOpqua directory as
# julia --project=. examples/pathogen_evolution/pathogen_evolution.jl
# unless viewing flamegraph, then run from console

# using Revise
using jOpqua

using StaticArrays
using Random

using Distances

using BenchmarkTools
using ProfileView

# Model setup
function run(seed::Int64, t_vec::Vector{Float64})
    # Parameters
    start_genome = "AAAAAAAA"
    optimal_genome = "BBBBBBBB"

    pat_type = jOpqua.newPathogenType(
        "pat_type",
        num_loci=length(optimal_genome),
        possible_alleles="AB",
        # contactSpecificCoefficient=s::String -> 1.0 + (0.1 * (length(optimal_genome) - hamming(s, optimal_genome)) / length(optimal_genome)),
        # receiveContactHostwideCoefficient=s::String -> 0.0,

    )

    hos_type = jOpqua.newHostType("hos_type")

    pop_type = jOpqua.newPopulationType(
        "pop_type",
        clearance_coefficient=1.0,
        contact_coefficient=1.05,
        receive_contact_coefficient=1.0,
        # mutations_upon_infection_coefficient=0.0005,
        # mutant_establishment_coefficient=0.0004,
        recombinant_establishment_coefficient=0.1,
        inoculum_coefficient=1.0,
        # death_coefficient=0.001,
        # birth_coefficient=0.001,
        pathogenFractions=jOpqua.pathogenFractionsProportionalFitness,
    )

    num_hosts = 10000
    num_infected = Int(num_hosts * 0.05)
    host_genome = ""

    # Setup
    model = jOpqua.newModel()
    pop = jOpqua.newPopulation!("pop", pop_type, model)
    jOpqua.addHostsToPopulation!(num_hosts, host_genome, hos_type, pop, model)
    pat = jOpqua.newPathogen!(start_genome, pop, pat_type)
    pat2 = jOpqua.newPathogen!(optimal_genome, pop, pat_type)

    for h in 1:num_infected
        jOpqua.addPathogenToHost!(pat, h, pop, model)
    end
    for h in num_infected+1:2*num_infected
        jOpqua.addPathogenToHost!(pat2, h, pop, model)
    end

    # Simulate
    Random.seed!(seed)
    model, output = jOpqua.simulate!(
        model, t_vec, population_host_samples=Dict("pop" => 200)
    )

    # Data output and plots
    compartment_data = jOpqua.saveCompartments(output, "examples/pathogen_evolution/compartment_pathogen_evolution.csv")
    jOpqua.plotCompartments(compartment_data, ["pop"], "examples/pathogen_evolution/compartment_pathogen_evolution.png")

    his_dat = jOpqua.saveHistory(output, "examples/pathogen_evolution/history_pathogen_evolution.csv")
    composition_data = jOpqua.saveComposition(
        his_dat, "examples/pathogen_evolution/composition_pathogen_evolution.csv",
        num_top_sequences=7, track_specific_sequences=["AAAA", "BBBB"]
    )
    jOpqua.plotComposition(
        composition_data, "examples/pathogen_evolution/composition_pathogen_evolution.png",
        normalized=true, ylabel="Fraction",
    )

    nwks = jOpqua.saveNewick(output, "examples/pathogen_evolution/pathogen_newick_pathogen_evolution.nwk", branch_length="Time", info_separator=" ")
    for nwk in nwks
        jOpqua.plotPhylogeny(nwk, "examples/pathogen_evolution/pathogen_newick_pathogen_evolution.png")
    end
end

run(1, collect(0.0:2.0:4.0)) # compile
@time run(0, collect(0.0:2.0:1500.0))
# Total events: 3173812
# 7.122476 seconds (71.25 M allocations: 8.838 GiB, 16.15% gc time, 7.82% compilation time: <1% of which was recompilation)
# 24 May 2026 Julia 1.12.6 Apple M3 Max 128 GB RAM

# @profview run(2, collect(0.0:2.0:1500.0))
