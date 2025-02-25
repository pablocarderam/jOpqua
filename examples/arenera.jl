# Run from base jOpqua directory as
# julia --project=. examples/arenera.jl
# unless viewing flamegraph, then run from console

using Revise
using jOpqua

using StaticArrays
using Random

using Distances

using BenchmarkTools
using ProfileView

# Model setup
function testRun(seed::Int64)
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

    res_type = jOpqua.newResponseType("res_type")
    pop_type = jOpqua.newPopulationType(
        "pop_type",
        clearance_coefficient=1.0,
        contact_coefficient=1.05,
        receive_contact_coefficient=1.0,
        pathogenFractions=jOpqua.pathogenFractionsProportionalFitness,
    )

    num_hosts = 10000
    num_infected = Int(num_hosts * 0.05)
    num_immune = 0
    host_genome = ""

    model = jOpqua.newModel()
    pop = jOpqua.newPopulation!("pop", pop_type, model)
    jOpqua.addHostsToPopulation!(num_hosts, host_genome, pop, model)
    pat = jOpqua.newPathogen!(start_genome, pop, pat_type)
    res = jOpqua.newResponse!(pat, pat, host_genome, pop, res_type)

    for h in 1:num_infected
        jOpqua.addPathogenToHost!(pat, h, pop, model)
    end
    for h in 1:num_immune
        jOpqua.addResponseToHost!(res, h, pop, model)
    end

    t_vec = collect(0.0:2.0:1500.0)

    Random.seed!(seed)

    # @profview , @time
    model, output = jOpqua.simulate!(
        model, t_vec, population_host_samples=Dict("pop" => 200)
    )
    println(output.compartment_vars["pop"][:, end])
    # println(output.host_samples["pop"][:, end][1:3])

    compartment_data = jOpqua.saveCompartments(output, "examples/compartment_test.csv")
    jOpqua.plotCompartments(compartment_data, ["pop"], "examples/compartment_test.png")

    his_dat = jOpqua.saveHistory(output, "examples/history_test.csv")
    composition_data = jOpqua.saveComposition(
        his_dat, "examples/composition_test.csv",
        num_top_sequences=7, track_specific_sequences=["AAAA", "BBBB"]
    )
    jOpqua.plotComposition(
        composition_data, "examples/composition_test.png",
        normalized=true, ylabel="Fraction",
    )

    nwks = jOpqua.saveNewick(output, "examples/pathogen_newick_test.nwk")
    for nwk in nwks
        jOpqua.plotPhylogeny(nwk, "examples/pathogen_newick_test.png")
    end
end

@time testRun(1)

# @profview testRun(0)
@time testRun(0)

# Result M3 Max 64 GB 9 Feb (second run) seed 0:
# 94438
#   2.023602 seconds (4.98 M allocations: 10.808 GiB, 18.23% gc time)
# [0.0, 465.0, 0.0, 0.0, 488.2500000000037, 0.0, 0.0, 0.0, 0.0]
#
# Same machine, same code, different day—what changed???:
# 94438
#   1.926932 seconds (4.98 M allocations: 10.808 GiB, 17.70% gc time)
# [0.0, 465.0, 0.0, 0.0, 488.2500000000037, 0.0, 0.0, 0.0, 0.0]

# Result M3 Max 64 GB 11 Feb (second run) seed 0:
# 94438
#   0.550135 seconds (4.55 M allocations: 170.740 MiB, 1.45% gc time)
# [0.0, 465.0, 0.0, 0.0, 488.2500000000037, 0.0, 0.0, 0.0, 0.0]
