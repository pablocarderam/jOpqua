using StaticArrays

const DEFAULT_PATHOGEN_TYPE = PathogenType(
    "Default",
    10,
    "ARNDCEQGHILKMFPSTWYV*",
    1.0,
    0.0,
    0.0,
    0.0,
    g -> 1.0, g -> 1.0, g -> 1.0, g -> 1.0,
    SA[ # order defined in COEFFICIENTS
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
    ],
)

const DEFAULT_RESPONSE_TYPE = ResponseType(
    "Default",
    0.0,
    (imp_g::String, mat_g::String, pat_g::String) -> 1.0,
    (imp_g::String, mat_g::String, pat_g::String) -> 1.0,
    SA[ # order defined in COEFFICIENTS
        (imp_g::String, mat_g::String)->1.0, (imp_g::String, mat_g::String)->1.0, (imp_g::String, mat_g::String)->1.0,
        (imp_g::String, mat_g::String)->1.0, (imp_g::String, mat_g::String)->1.0, (imp_g::String, mat_g::String)->1.0,
        (imp_g::String, mat_g::String)->1.0, (imp_g::String, mat_g::String)->1.0, (imp_g::String, mat_g::String)->1.0,
        (imp_g::String, mat_g::String)->1.0, (imp_g::String, mat_g::String)->1.0, (imp_g::String, mat_g::String)->1.0,
    ],
    SA[ # order defined in COEFFICIENTS
        (imp_g::String, mat_g::String, pat_g::String)->1.0, (imp_g::String, mat_g::String, pat_g::String)->1.0,
        (imp_g::String, mat_g::String, pat_g::String)->1.0, (imp_g::String, mat_g::String, pat_g::String)->1.0,
        (imp_g::String, mat_g::String, pat_g::String)->1.0, (imp_g::String, mat_g::String, pat_g::String)->1.0,
        (imp_g::String, mat_g::String, pat_g::String)->1.0, (imp_g::String, mat_g::String, pat_g::String)->1.0,
        (imp_g::String, mat_g::String, pat_g::String)->1.0, (imp_g::String, mat_g::String, pat_g::String)->1.0,
        (imp_g::String, mat_g::String, pat_g::String)->1.0, (imp_g::String, mat_g::String, pat_g::String)->1.0,
    ],
)

const DEFAULT_POPULATION_TYPE = PopulationType(
    "Default",
    true, false,
    1.0,
    1.0,
    1.0,
    SA[ # order defined in COEFFICIENTS
        0.0, 1.0, 0.0,
        0.0, 1.05, 0.0,
        0.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
    ],
    pathogenFractionsWinnerTakesAll,
    weightedResponseWinnerTakesAll,
    infectionProbabilityWinnerTakesAll,
    (p::Pathogen, h::Host, r::Vector{Response}) -> Vector{Response}(undef, 0),
)
