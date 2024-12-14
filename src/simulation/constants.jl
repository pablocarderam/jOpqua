using StaticArrays

# Pathogen events
const MUTANT_ESTABLISHMENT = 1
const MUTANT_TRANSMISSION = 2
const RECOMBINATION = 3
const TRANSMISSION = 4
const RECOVERY = 5

# Immunity events
const IMMUNITY_LOSS = 6

# Host events
const BIRTH = 7
const DEATH = 8
const MIGRATION = 9
const CLASS_CHANGE = 10

# Event modifiers
const SUSCEPTIBILITY = 11

# Global trackers
const PATHOGEN_EVENTS = SA[1, 2, 3, 4, 5]
const NUM_PATHOGEN_EVENTS = length(PATHOGEN_EVENTS)

const IMMUNITY_EVENTS = SA[6]
const NUM_IMMUNITY_EVENTS = length(IMMUNITY_EVENTS)

const HOST_EVENTS = SA[7, 8, 9, 10]
const NUM_HOST_EVENTS = length(HOST_EVENTS)

const EVENTS = SA[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
const NUM_EVENTS = length(EVENTS)

const COEFFICIENTS = SA[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
const NUM_COEFFICIENTS = length(COEFFICIENTS)

# Array sizes
# const MAX_PATHOGENS = 100
# const MAX_IMMUNITIES = 1000
# const MAX_HOSTS = 100
# const CLASSES = 100
# const POPULATIONS = 100
