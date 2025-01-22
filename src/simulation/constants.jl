using StaticArrays

# Pathogen events, zooming out in scale (order matters)
# Single pathogen, single host
const MUTANT_ESTABLISHMENT = 1
const CLEARANCE = 2
const RESPONSE_ACQUISITION = 3
# Two pathogens, single host
const RECOMBINANT_ESTABLISHMENT = 4
# Single pathogen, two hosts
const INTRA_POPULATION_CONTACT = 5
const INTER_POPULATION_CONTACT = 6

# Response events
# Single response, single host
const RESPONSE_LOSS = 7

# Host events
# Single host, single class
const BIRTH = 8
const DEATH = 9
# Single host, two classes
const CLASS_CHANGE = 10
# Single host, two populations
const MIGRATION = 11

# Choice modifiers, zooming in in scale (order matters)
# Population choice
const RECEIVE_MIGRATION = 12
# Class choice
const RECEIVE_CLASS_CHANGE = 13
# Host choice
const RECEIVE_CONTACT = 14
# Pathogen choice
const INTRAHOST_FITNESS = 15

# Global trackers
const PATHOGEN_EVENTS = SA[
    MUTANT_ESTABLISHMENT, CLEARANCE, RESPONSE_ACQUISITION,
    RECOMBINANT_ESTABLISHMENT, INTRA_POPULATION_CONTACT, INTER_POPULATION_CONTACT
]
const NUM_PATHOGEN_EVENTS = length(PATHOGEN_EVENTS)

const RESPONSE_EVENTS = SA[RESPONSE_LOSS]
const NUM_RESPONSE_EVENTS = length(RESPONSE_EVENTS)

const HOST_EVENTS = SA[BIRTH, DEATH, CLASS_CHANGE, MIGRATION]
const NUM_HOST_EVENTS = length(HOST_EVENTS)

const EVENTS = SA[PATHOGEN_EVENTS..., RESPONSE_EVENTS..., HOST_EVENTS...]
const NUM_EVENTS = length(EVENTS)

const CHOICE_MODIFIERS = SA[
    RECEIVE_MIGRATION, RECEIVE_CLASS_CHANGE, RECEIVE_CONTACT, INTRAHOST_FITNESS
]
const NUM_CHOICE_MODIFIERS = length(CHOICE_MODIFIERS)

const COEFFICIENTS = SA[EVENTS..., CHOICE_MODIFIERS...]
const NUM_COEFFICIENTS = length(COEFFICIENTS)

# Starter coefficients
const START_COEFFICIENTS = SVector{NUM_COEFFICIENTS,Float64}(
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0]
# class change, inter-population contact and migration numbers per class or
# population are fractions that sum to one, so no need to account for in here
)

# Misc constants
const CHROMOSOME_SEPARATOR = "/"
