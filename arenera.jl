using Revise
using jOpqua

using StaticArrays
using Random

# Parameters
pa_type = jOpqua.PathogenType(
    "pa_type",
    10,
    "ARNDCEQGHILKMFPSTWYV*",
    0,
    1.0,
    SA[ # pathogen_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # Each element takes seq argument, returns Float64
        g->1.0, g->1.0, g->1.0, g->1.0, g->1.0,
        g->1.0, g->1.0, g->1.0, g->1.0, g->1.0,
        g->1.0, g->1.0, g->1.0, g->1.0, g->1.0,
        g->1.0, g->1.0,
    ],
)

re_type = jOpqua.ResponseType(
    "re_type", # id::String,
    SA[ # static_coefficient_functions::SVector{NUM_COEFFICIENTS,Function},
        (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0,
        (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0,
        (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0,
        (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0,
    ],
    SA[ # specific_coefficient_functions::SVector{NUM_COEFFICIENTS,Function},
        (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0,
        (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0,
        (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0,
        (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0,
    ],
    (imp_g, mat_g, pat_g)->1.0
)

class_parameters = jOpqua.ClassParameters(
    "TestParams", # id::String
    SA[ # base_coefficients::SVector{NUM_COEFFICIENTS,Float64}
        0.0, 1.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0,
    ],
    Dict(), # class_change_fractions::Dict{String,Float64} # size CLASSES, must sum to 1
    Dict(), # inter_population_contact_fractions::Dict{String,Float64} # size POPULATIONS, must sum to 1
    Dict(), # migration_fractions::Dict{String,Float64} # size POPULATIONS, must sum to 1
    Dict(pa_type.id => pa_type), # response_types::Dict{String,ResponseType}
    Dict(re_type.id => re_type), # response_types::Dict{String,ResponseType}
    (p, h, c) -> Nothing # developResponses::Function
)

# Model setup
num_populations = 1
num_hosts = 10
num_infected = 5
num_immune = 0

model = jOpqua.newModel()
pop = jOpqua.newPopulation!("pop", model)
class = jOpqua.newClass!("class", class_parameters, pop)
for _ in 1:num_hosts
    host = jOpqua.newHost!(class, pop, model)
end
pat = jOpqua.newPathogen!("AAAA", class, pa_type)
res = jOpqua.newResponse!(pat, pat, (pat.sequence, pat.sequence, re_type.id), class, re_type)

for h in 1:num_infected
    jOpqua.addPathogenToHost!(pat, h, class, pop, model)
end
for h in 1:num_immune
    jOpqua.addResponseToHost!(pat, h, class, pop, model)
end

jOpqua.establishMutant!(model, rand())

# jOpqua.removePathogenFromHost!(1, 1, class, pop, model)
# jOpqua.removeResponseFromHost!(1, 1, class, pop, model)

# rand_n = rand()
# x, rand_n = jOpqua.choosePathogen(1, 1, 1, 1, model, rand_n)
# x, rand_n = jOpqua.chooseResponse(1, 1, 1, 7, model, rand_n)
# x, rand_n = jOpqua.chooseHost(1, 1, 7, model, rand_n)
# x, rand_n = jOpqua.chooseClass(1, 7, model, rand_n)
# x, rand_n = jOpqua.choosePopulation(7, model, rand_n)
# x, rand_n = jOpqua.chooseEvent(model, rand_n)

# jOpqua.addPathogenToHost!(pat, 1, class, pop, model)
# jOpqua.addResponseToHost!(res, 1, class, pop, model)

# x, rand_n = jOpqua.choosePathogen(1, 1, 1, 1, model, rand_n)
# x, rand_n = jOpqua.chooseResponse(1, 1, 1, 7, model, rand_n)
# x, rand_n = jOpqua.chooseHost(1, 1, 7, model, rand_n)
# x, rand_n = jOpqua.chooseClass(1, 7, model, rand_n)
# x, rand_n = jOpqua.choosePopulation(7, model, rand_n)
# x, rand_n = jOpqua.chooseEvent(model, rand_n)

# jOpqua.establishMutant!(model, rand())
# jOpqua.clearPathogen!(model, rand())
# jOpqua.acquireResponse!(model, rand())
# jOpqua.establishRecombinant!(model, rand())
# jOpqua.intraPopulationContact!(model, rand())
