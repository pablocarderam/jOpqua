# Run from base jOpqua directory as
# julia --project=. examples/arenera.jl

using Revise
using jOpqua

using StaticArrays
using Random

using BenchmarkTools
using ProfileView

# Model setup
function testRun(seed::Int64)
    # Parameters
    pat_type = jOpqua.newPathogenType("pat_type")
    res_type = jOpqua.newResponseType("res_type")
    pop_type = jOpqua.newPopulationType("pop_type")

    num_hosts = 10000
    num_infected = Int(num_hosts * 0.5)
    num_immune = 0

    model = jOpqua.newModel()
    pop = jOpqua.newPopulation!("pop", pop_type, model)
    jOpqua.addHostsToPopulation!(num_hosts, pop, model)
    pat = jOpqua.newPathogen!("AAAA", pop, pat_type)
    res = jOpqua.newResponse!(pat, pat, pop, res_type)

    for h in 1:num_infected
        jOpqua.addPathogenToHost!(pat, h, pop, model)
    end
    for h in 1:num_immune
        jOpqua.addResponseToHost!(res, h, pop, model)
    end

    t_vec = collect(0.0:0.1:50.0)

    Random.seed!(seed)

    # @profview , @time
    model, output = jOpqua.simulate!(
        model, t_vec, host_samples_population=Dict("pop"=>100)
    )
    println(output.compartment_vars["pop"][:,end])
    println(output.host_samples["pop"][:,end][1:3])

    jOpqua.plotCompartments(output, ["pop"], "examples/compartment_test.png")
end

@time testRun(1)
@time testRun(0)
@profview testRun(0)

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
