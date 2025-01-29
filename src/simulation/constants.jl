using StaticArrays

# Pathogen events, zooming out in scale (order matters)
# Single pathogen, single host
const MUTANT_ESTABLISHMENT = 1
const CLEARANCE = 2
const RESPONSE_ACQUISITION = 3
# Two pathogens, single host
const RECOMBINANT_ESTABLISHMENT = 4
# Single pathogen, two hosts
const CONTACT = 5

# Response events
# Single response, single host
const RESPONSE_LOSS = 6

# Host events
# Single host, single population
const BIRTH = 7
const DEATH = 8
# Single host, two populations
const MIGRATION = 9

# Choice modifiers, zooming in in scale (order matters)
# Population choice
const RECEIVE_MIGRATION = 10
# Host choice
const RECEIVE_CONTACT = 11
# Pathogen choice
const INTRAHOST_FITNESS = 12

# Global trackers
const PATHOGEN_EVENTS = SA[
    MUTANT_ESTABLISHMENT, CLEARANCE, RESPONSE_ACQUISITION,
    RECOMBINANT_ESTABLISHMENT, CONTACT
]
const NUM_PATHOGEN_EVENTS = length(PATHOGEN_EVENTS)

const RESPONSE_EVENTS = SA[RESPONSE_LOSS]
const NUM_RESPONSE_EVENTS = length(RESPONSE_EVENTS)

const HOST_EVENTS = SA[BIRTH, DEATH, MIGRATION]
const NUM_HOST_EVENTS = length(HOST_EVENTS)

const EVENTS = SA[PATHOGEN_EVENTS..., RESPONSE_EVENTS..., HOST_EVENTS...]
const NUM_EVENTS = length(EVENTS)

const CHOICE_MODIFIERS = SA[
    RECEIVE_MIGRATION, RECEIVE_CONTACT, INTRAHOST_FITNESS
]
const NUM_CHOICE_MODIFIERS = length(CHOICE_MODIFIERS)

const COEFFICIENTS = SA[EVENTS..., CHOICE_MODIFIERS...]
const NUM_COEFFICIENTS = length(COEFFICIENTS)

# Starter coefficients
const START_COEFFICIENTS = SVector{NUM_COEFFICIENTS,Float64}(
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0]
)

# Misc constants
const CHROMOSOME_SEPARATOR = "/"
