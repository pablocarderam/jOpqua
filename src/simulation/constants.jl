using StaticArrays

# Pathogen events, zooming out in scale (order matters)
# Single pathogen, single host
const MUTANT_ESTABLISHMENT = 1
const CLEARANCE = 2
const RESPONSE_ACQUISITION = 3
# Two pathogens, single host
const RECOMBINATION = 4
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
# Event choice
const MUTATION_AT_CONTACT = 12
# Population choice
const RECEIVE_MIGRATION = 13
# Class choice
const RECEIVE_CLASS_CHANGE = 14
# Host choice
const RECEIVE_CONTACT = 15
# Pathogen choice
const INTRAHOST_FITNESS = 16

# Global trackers
const PATHOGEN_EVENTS = SA[
    MUTANT_ESTABLISHMENT, CLEARANCE, RESPONSE_ACQUISITION,
    RECOMBINATION, INTRA_POPULATION_CONTACT, INTER_POPULATION_CONTACT
]
const NUM_PATHOGEN_EVENTS = length(PATHOGEN_EVENTS)

const RESPONSE_EVENTS = SA[RESPONSE_LOSS]
const NUM_RESPONSE_EVENTS = length(RESPONSE_EVENTS)

const HOST_EVENTS = SA[BIRTH, DEATH, CLASS_CHANGE, MIGRATION]
const NUM_HOST_EVENTS = length(HOST_EVENTS)

const EVENTS = SA[PATHOGEN_EVENTS..., RESPONSE_EVENTS..., HOST_EVENTS...]
const NUM_EVENTS = length(EVENTS)

const CHOICE_MODIFIERS = SA[
    MUTATION_AT_CONTACT, RECEIVE_MIGRATION,
    RECEIVE_CLASS_CHANGE, RECEIVE_CONTACT, INTRAHOST_FITNESS
]
const NUM_CHOICE_MODIFIERS = length(CHOICE_MODIFIERS)

const COEFFICIENTS = SA[EVENTS..., CHOICE_MODIFIERS...]
const NUM_COEFFICIENTS = length(COEFFICIENTS)
