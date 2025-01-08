using Revise
using jOpqua

using StaticArrays

# Parameters
im_type = jOpqua.ImmunityType(
    "im_type", # id::String,
    SA[ # static_coefficient_functions::SVector{NUM_COEFFICIENTS,Function},
        (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0,
        (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0,
        (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0, (imp_g, mat_g)->1.0,
    ],
    SA[ # specific_coefficient_functions::SVector{NUM_COEFFICIENTS,Function},
        (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0,
        (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0,
        (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0, (imp_g, mat_g, pat_g)->1.0,
    ],
    (imp_g, mat_g, pat_g) -> 1.0, # immunodominance::Function,
)

class_parameters = jOpqua.ClassParameters(
    "TestParams", # id::String
    SA[ # base_coefficients::SVector{NUM_COEFFICIENTS,Float64}
        0.0, 1.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 0.0
    ],
    Dict(), # class_change_fractions::Dict{String,Float64} # size CLASSES, must sum to 1
    Dict(), # inter_population_contact_fractions::Dict{String,Float64} # size POPULATIONS, must sum to 1
    Dict(), # migration_fractions::Dict{String,Float64} # size POPULATIONS, must sum to 1
    SA[ # pathogen_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # Each element takes seq argument, returns Float64
        g->1.0, g->1.0, g->1.0, g->1.0, g->1.0,
        g->1.0, g->1.0, g->1.0, g->1.0, g->1.0,
        g->1.0, g->1.0, g->1.0, g->1.0, g->1.0,
    ],
    Dict(im_type.id => im_type), # immunity_types::Dict{String,ImmunityType}
    (p, h, c) -> Nothing # acquireImmunities::Function
)

# Model setup
num_populations = 1

model = jOpqua.newModel()
pop = jOpqua.newPopulation!("pop", model)
class = jOpqua.newClass!("class", class_parameters, pop)
host = jOpqua.newHost!(class, pop, model)
pat = jOpqua.newPathogen!("ATCG", class)
imm = jOpqua.newImmunity!(pat, pat, class, im_type)

jOpqua.infect!(pat, host, class, pop, model)
jOpqua.immunize!(imm, host, class, pop, model)

jOpqua.clear!(1, host, class, pop, model)
jOpqua.deimmunize!(1, host, class, pop, model)
