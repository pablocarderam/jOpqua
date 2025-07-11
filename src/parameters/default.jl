using StaticArrays

const DEFAULT_PATHOGEN_TYPE = PathogenType(
    "Default",
    0,
    "", #"ARNDCEQGHILKMFPSTWYV*",
    SA[ # order defined in COEFFICIENTS
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
    ],
    SA[ # order defined in COEFFICIENTS
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
    ],
)

const DEFAULT_RESPONSE_TYPE = ResponseType(
    "Default",
    (hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> 1.0,
    SA[ # order defined in COEFFICIENTS
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0,
    ],
    SA[ # order defined in COEFFICIENTS
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
    ],
    SA[ # order defined in COEFFICIENTS
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String)->1.0,
    ],
    SA[ # order defined in COEFFICIENTS
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0, (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
        (hos_g::String, imp_g::String, mat_g::String, pat_g::String)->1.0,
    ],
)

const DEFAULT_HOST_TYPE = HostType(
    "Default",
    0,
    "", #"ARNDCEQGHILKMFPSTWYV*",
    SA[ # order defined in COEFFICIENTS
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
        g::String->1.0, g::String->1.0, g::String->1.0,
    ],
)

const DEFAULT_POPULATION_TYPE = PopulationType(
    "Default",
    true,
    false,
    true,
    (h_1::String, h_2::String) -> true,
    SA[ # order defined in COEFFICIENTS
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
    ],
    pathogenFractionsWinnerTakesAll,
    weightedInteractionPathogenWinnerTakesAll,
    weightedInteractionResponseWinnerTakesAll,
    weightedInteractionHostwideProduct,
    (
        pathogen::Pathogen, host::Host,
        existing_responses::Dict{Tuple{String,String,String,String},Response},
        response_types::Dict{String,ResponseType},
        birth_time::Float64
    ) -> [deNovoResponse(
        pathogen, host, existing_responses, response_types, birth_time
    )],
    Dict{String,ResponseType}([(DEFAULT_RESPONSE_TYPE.id => DEFAULT_RESPONSE_TYPE)])
)
