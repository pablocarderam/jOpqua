using StaticArrays

const DEFAULT_PATHOGEN_TYPE = PathogenType(
    "Default",
    10,
    "ARNDCEQGHILKMFPSTWYV*",
    1.0,
    0.0,
    0.0,
    0.0,
    g->1.0, g->1.0, g->1.0, g->1.0,
    SA[ # order defined in COEFFICIENTS
        g->1.0, g->1.0, g->1.0,
        g->1.0, g->1.0, g->1.0,
        g->1.0, g->1.0, g->1.0,
        g->1.0, g->1.0, g->1.0
    ],
)

const DEFAULT_RESPONSE_TYPE = ResponseType(
    "Default",
    0.0,
    (imp_g, mat_g, pat_g) -> 1.0,
    (imp_g, mat_g, pat_g) -> 1.0,
    SA[ # order defined in COEFFICIENTS
        (imp_g, mat_g) -> 1.0, (imp_g, mat_g) -> 1.0, (imp_g, mat_g) -> 1.0,
        (imp_g, mat_g) -> 1.0, (imp_g, mat_g) -> 1.0, (imp_g, mat_g) -> 1.0,
        (imp_g, mat_g) -> 1.0, (imp_g, mat_g) -> 1.0, (imp_g, mat_g) -> 1.0,
        (imp_g, mat_g) -> 1.0, (imp_g, mat_g) -> 1.0, (imp_g, mat_g) -> 1.0,
    ],
    SA[ # order defined in COEFFICIENTS
        (imp_g, mat_g, pat_g) -> 1.0, (imp_g, mat_g, pat_g) -> 1.0, (imp_g, mat_g, pat_g) -> 1.0,
        (imp_g, mat_g, pat_g) -> 1.0, (imp_g, mat_g, pat_g) -> 1.0, (imp_g, mat_g, pat_g) -> 1.0,
        (imp_g, mat_g, pat_g) -> 1.0, (imp_g, mat_g, pat_g) -> 1.0, (imp_g, mat_g, pat_g) -> 1.0,
        (imp_g, mat_g, pat_g) -> 1.0, (imp_g, mat_g, pat_g) -> 1.0, (imp_g, mat_g, pat_g) -> 1.0
    ],
)

const DEFAULT_POPULATION_TYPE = PopulationType(
    "Default", # id::String
    true, false,
    1.0, # inoculum_coefficient
    1.0, # mutation_coefficient
    1.0, # recombination_coefficient
    SA[ # base_coefficients::SVector{NUM_COEFFICIENTS,Float64}
        0.0, 1.0, 0.0,
        0.0, 1.05, 0.0,
        0.0, 0.0, 0.0,
        0.0, 1.0, 0.0
    ],
    pathogenFractionsWinnerTakesAll,
    weightedResponseArithmeticMean,
    infectionProbabilityArithmeticMean,
    (p, h, c) -> Nothing, # developResponses::Function,
)
