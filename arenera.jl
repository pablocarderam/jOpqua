using Revise
using jOpqua

using StaticArrays
using Random

# Parameters
pa_type = jOpqua.PathogenType(
    "pa_type",
    10,
    "ARNDCEQGHILKMFPSTWYV*",
    1.0,
    0.0,
    0.0,
    0.0,
    SA[ # pathogen_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # Each element takes seq argument, returns Float64
        g->1.0, g->1.0, g->1.0, g->1.0,
        g->1.0, g->1.0, g->1.0, g->1.0,
        g->1.0, g->1.0, g->1.0, g->1.0,
    ],
    g -> 1.0, g -> 1.0, g -> 1.0, g -> 1.0,
)

re_type = jOpqua.ResponseType(
    "re_type", # id::String,
    SA[ # static_coefficient_functions::SVector{NUM_COEFFICIENTS,Function},
        (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0,
        (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0,
        (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0,
    ],
    SA[ # specific_coefficient_functions::SVector{NUM_COEFFICIENTS,Function},
        (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0,
        (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0,
        (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0,
    ],
    0.0,
    (imp_g, mat_g, pat_g) -> 1.0,
    (imp_g, mat_g, pat_g) -> 1.0,
)

pop_parameters = jOpqua.PopulationType(
    "TestParams", # id::String
    SA[ # base_coefficients::SVector{NUM_COEFFICIENTS,Float64}
        0.00, 1.0, 0.0, 0.0,
        1.05, 0.0, 0.0, 0.0,
        0.00, 0.0, 1.0, 0.0,
    ],
    true, false,
    jOpqua.pathogenFractionsWinnerTakesAll,
    jOpqua.weightedResponseArithmeticMean,
    jOpqua.infectionProbabilityArithmeticMean,
    (p, h, c) -> Nothing, # developResponses::Function,
    1.0, # inoculum_coefficient
    1.0, # mutation_coefficient
    1.0, # recombination_coefficient
)

# Model setup
num_populations = 1
num_hosts = 10000
num_infected = Int(num_hosts * 0.5)
num_immune = 0

model = jOpqua.newModel()
pop = jOpqua.newPopulation!("pop", pop_parameters, model)
for i in 1:num_hosts
    # println(i)
    host = jOpqua.newHost!(pop, model)
end
pat = jOpqua.newPathogen!("AAAA", pop, pa_type)
res = jOpqua.newResponse!(pat, pat, (pat.sequence, pat.sequence, re_type.id), pop, re_type)

for h in 1:num_infected
    # println(h)
    jOpqua.addPathogenToHost!(pat, h, pop, model)
end
for h in 1:num_immune
    jOpqua.addResponseToHost!(pat, h, pop, model)
end

t_vec = collect(0.0:50.0)

@time jOpqua.simulate!(model, t_vec)
model.event_rates

# jOpqua.establishMutant!(model, rand())

# jOpqua.removePathogenFromHost!(1, 1, pop, model)
# jOpqua.removeResponseFromHost!(1, 1, pop, model)

# rand_n = rand()
# x, rand_n = jOpqua.choosePathogen(1, 1, 1, 1, model, rand_n)
# x, rand_n = jOpqua.chooseResponse(1, 1, 1, 7, model, rand_n)
# x, rand_n = jOpqua.chooseHost(1, 1, 7, model, rand_n)
# x, rand_n = jOpqua.choosePopulation(7, model, rand_n)
# x, rand_n = jOpqua.chooseEvent(model, rand_n)

# jOpqua.addPathogenToHost!(pat, 1, pop, model)
# jOpqua.addResponseToHost!(res, 1, pop, model)

# x, rand_n = jOpqua.choosePathogen(1, 1, 1, 1, model, rand_n)
# x, rand_n = jOpqua.chooseResponse(1, 1, 1, 7, model, rand_n)
# x, rand_n = jOpqua.chooseHost(1, 1, 7, model, rand_n)
# x, rand_n = jOpqua.choosePopulation(7, model, rand_n)
# x, rand_n = jOpqua.chooseEvent(model, rand_n)

# jOpqua.establishMutant!(model, rand())
# jOpqua.clearPathogen!(model, rand())
# jOpqua.acquireResponse!(model, rand())
# jOpqua.establishRecombinant!(model, rand())
# jOpqua.intraPopulationContact!(model, rand())
