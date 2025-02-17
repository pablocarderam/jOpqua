using StaticArrays

function newPathogenType(id::String;)
    return PathogenType(
        id,
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
end

function newResponseType(id::String;)
    return ResponseType(
        id, # id::String,
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
end

function newPopulationType(id::String;)
    return PopulationType(
        id, # id::String
        SA[ # base_coefficients::SVector{NUM_COEFFICIENTS,Float64}
            0.00, 1.0, 0.0, 0.0,
            1.05, 0.0, 0.0, 0.0,
            0.00, 0.0, 1.0, 0.0,
        ],
        true, false,
        pathogenFractionsWinnerTakesAll,
        weightedResponseArithmeticMean,
        infectionProbabilityArithmeticMean,
        (p, h, c) -> Nothing, # developResponses::Function,
        1.0, # inoculum_coefficient
        1.0, # mutation_coefficient
        1.0, # recombination_coefficient
    )
end
