using StaticArrays

# Pathogen events, zooming out in scale (order matters)
# Single pathogen, single host
const MUTANT_ESTABLISHMENT = 1
const CLEARANCE = 2
# Two pathogens, single host
const RECOMBINATION = 3
# Single pathogen, two hosts
const INTRA_POPULATION_CONTACT = 4
const INTER_POPULATION_CONTACT = 5

# Immunity events
# Single immunity, single host
const IMMUNITY_LOSS = 6

# Host events
# Single host, single class
const BIRTH = 7
const DEATH = 8
# Single host, two classes
const CLASS_CHANGE = 9
# Single host, two populations
const MIGRATION = 10

# Choice modifiers, zooming in in scale (order matters)
# Event choice
const MUTATION_AT_CONTACT = 11
# Population choice
const RECEIVE_MIGRATION = 12
# Class choice
const RECEIVE_CLASS_CHANGE = 13
# Host choice
const RECEIVE_CONTACT = 14
# Pathogen choice
const INTRAHOST_FITNESS = 15

# Global trackers
const PATHOGEN_EVENTS = SA[1, 2, 3, 4, 5]
const NUM_PATHOGEN_EVENTS = length(PATHOGEN_EVENTS)

const IMMUNITY_EVENTS = SA[6]
const NUM_IMMUNITY_EVENTS = length(IMMUNITY_EVENTS)

const HOST_EVENTS = SA[7, 8, 9, 10]
const NUM_HOST_EVENTS = length(HOST_EVENTS)

const EVENTS = SA[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
const NUM_EVENTS = length(EVENTS)

const CHOICE_MODIFIERS = SA[11, 12, 13, 14, 15]
const NUM_CHOICE_MODIFIERS = length(CHOICE_MODIFIERS)

const COEFFICIENTS = SA[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
const NUM_COEFFICIENTS = length(COEFFICIENTS)

# Array sizes
# const MAX_PATHOGENS = 100
# const MAX_IMMUNITIES = 1000
# const MAX_HOSTS = 100
# const CLASSES = 100
# const POPULATIONS = 100
