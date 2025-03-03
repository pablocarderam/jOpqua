using StaticArrays

# 1. Sampling variables and events
# 1.1 Pathogen events, zooming out in scale (order matters)
# 1.1.1 Single pathogen, single host
const MUTANT_ESTABLISHMENT = 1
const CLEARANCE = 2
const RESPONSE_ACQUISITION = 3
# 1.1.2 Two pathogens (symmetrical), single host
const RECOMBINANT_ESTABLISHMENT = 4
# 1.1.3 Single pathogen, two hosts (asymmetrical)
const CONTACT = 5

# 1.2 Response events
# 1.2.1 Single response, single host
const RESPONSE_LOSS = 6

# 1.3 Host events
# 1.3.1 Single host, single population
const BIRTH = 7
const DEATH = 8
# 1.3.2 Single host, two populations (asymmetrical)
const TRANSITION = 9

# 1.4 Choice modifiers and second entity samplers, zooming in in scale (order matters)
# 1.4.1 Population choice
const RECEIVE_TRANSITION = 10

# 1.4.2 Host choice
const RECEIVE_CONTACT = 11

# 1.4.3 Pathogen choice
const INTRAHOST_FITNESS = 12

# 2. Non-sampling variables and sub-events, zooming out in scale
# 2.1 Pathogen-specific sub-events
const MUTATIONS_UPON_INFECTION = 13
const RECOMBINATIONS_UPON_INFECTION = 14
const INOCULUM = 15
const TRANSMISSION_EFFICIENCY = 16
const VERTICAL_TRANSMISSION = 17
const RESPONSE_ACQUISITION_UPON_CLEARANCE = 18

# 2.2 Response-specific sub-events
const RESPONSE_INHERITANCE = 19

# 2.3 Host-specific sub-events
const HOST_MUTATIONS_UPON_BIRTH = 20
const HOST_RECOMBINATIONS_UPON_BIRTH = 21

# Global trackers
const PATHOGEN_EVENTS = SA[
    MUTANT_ESTABLISHMENT, CLEARANCE, RESPONSE_ACQUISITION,
    RECOMBINANT_ESTABLISHMENT, CONTACT
]
const NUM_PATHOGEN_EVENTS = length(PATHOGEN_EVENTS)

const RESPONSE_EVENTS = SA[RESPONSE_LOSS]
const NUM_RESPONSE_EVENTS = length(RESPONSE_EVENTS)

const HOST_EVENTS = SA[BIRTH, DEATH, TRANSITION]
const NUM_HOST_EVENTS = length(HOST_EVENTS)

const EVENTS = SA[PATHOGEN_EVENTS..., RESPONSE_EVENTS..., HOST_EVENTS...]
const NUM_EVENTS = length(EVENTS)

const CHOICE_MODIFIERS = SA[
    RECEIVE_TRANSITION, RECEIVE_CONTACT, INTRAHOST_FITNESS
]
const NUM_CHOICE_MODIFIERS = length(CHOICE_MODIFIERS)

const SAMPLING_COEFFICIENTS = SA[EVENTS..., CHOICE_MODIFIERS...]
const NUM_SAMPLING_COEFFICIENTS = length(SAMPLING_COEFFICIENTS)

const PATHOGEN_NONSAMPLING_COEFFICIENTS = SA[
    MUTATIONS_UPON_INFECTION, RECOMBINATIONS_UPON_INFECTION, INOCULUM,
    TRANSMISSION_EFFICIENCY, VERTICAL_TRANSMISSION,
    RESPONSE_ACQUISITION_UPON_CLEARANCE
]

const RESPONSE_NONSAMPLING_COEFFICIENTS = SA[RESPONSE_INHERITANCE]

const HOST_NONSAMPLING_COEFFICIENTS = SA[
    HOST_MUTATIONS_UPON_BIRTH, HOST_RECOMBINATIONS_UPON_BIRTH
]

const NON_SAMPLING_COEFFICIENTS = SA[
    PATHOGEN_NONSAMPLING_COEFFICIENTS...,
    RESPONSE_NONSAMPLING_COEFFICIENTS..., HOST_NONSAMPLING_COEFFICIENTS...
]

const COEFFICIENTS = SA[SAMPLING_COEFFICIENTS..., NON_SAMPLING_COEFFICIENTS...]
const NUM_COEFFICIENTS = length(COEFFICIENTS)

# Starter coefficients
const START_COEFFICIENTS = SVector{NUM_COEFFICIENTS,Float64}([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    ])

# History variable trackers
const UNINFECTED_NAIVE = 1
const UNINFECTED_IMMUNE = 2
const INFECTED_NAIVE = 3
const INFECTED_IMMUNE = 4
const DEAD = 5
const COMPARTMENTS = SA[
    UNINFECTED_NAIVE, UNINFECTED_IMMUNE,
    INFECTED_NAIVE, INFECTED_IMMUNE, DEAD
]
const NUM_COMPARTMENTS = length(COMPARTMENTS)

# Misc model constants
const CHROMOSOME_SEPARATOR = "/"
const HOMOLOGOUS_CHROMOSOME_SEPARATOR = "|"

# Output constants
const WITHIN_HOST_SEPARATOR = ";"
const PARENT_SEPARATOR = "+"
const COMPARTMENT_LABELS = [
    "Uninfected Naive", "Uninfected Immune", "Infected Naive",
    "Infected Immune", "Dead"
]

# Simulation constants
const ERROR_TOLERANCE = 1.0e-9
