# jOpqua Changelog

## 12 February 2025
- Removed `total_hosts`, replaced with instances of `length(population.hosts)` (PCR)
- Optimized `Host` initialization at startup (CLM)
- Made `event_functions` into constant `EVENT_FUNCTIONS` (PCR)

## 11 February 2025
- Added interventions (PCR)
- Added flamegraphs to sandbox file (PCR)
- Defined a couple of static array types within `events.jl` (PCR)
- Made new matrices in `Population` with sums of host weights pre-multiplied by
coefficients, to be used in `chooseHost` replacing the inline, broadcasted,
element-wise multiplication by the base coefficient; resulted in a ~4X speedup and
garbage collection time going from ~18% to ~1.5%, even with the speedup--the problem
was found with the dual wonders of flame graphs and CLM (PCR)
- Defined type of `AbstractArray` in `randChoose` (PCR)

There was a bug in the use of parenthesis for `@views` in `choice.jl`, very
minor effect on performance. The fix was made obsolete by removing the broadcasted
multiplication anyway.

Next is saving history--a clean implementation of this would take advantage of
the intervention framework to save whatever is needed whenever needed, but might
result in performance issues given the type instability of interventions, so maybe
a standalone implementation of sampling and saving might be best?

## 9 February 2025
- Small bug fix in `transition!` (PCR)

## 5 February 2025
- Added pathogen fraction to calculation of vertical transmission probability (PCR)
- Moved code to add a `Host` to a `Population` into a separate function (PCR)
- Optimizations including: add types to global variables; match literal types to
their context; extend use of `@views`; cache emitting `Host` in `hostContact!` (CLM)
- Backtrack 2 optimizations (one `@views` use with no slicing, one comparison between
ints that did not require using a float zero) (PCR)
- Added `removeHostFromPopulation` and death event (PCR)
- Added transition event (PCR)

All events implemented!!! but not tested lol

TODO:
- function for running and benchmarking repeated simulations; current setup
has a global `model` variable which I (CLM) believe retains its state after
simulation _n_ at the beginning of simulation _n+1_, affecting both simulation
trajectory and runtime
- Test and debug the **** out of everything (all parameters and events,
particularly events other than clearance and contact)

Dev roadmap:
- implement interventions
- simulation history saving
- input parameters
- output data + plotting
- immunity DLC

## 4 February 2025
- Added `constant_contact_density` as a parameter of `PopulationType` to specify
whether or not to normalize by host population size in calculation of
receive contact weights, equivalently, added `constant_transition_density`
for receive transition weights (`constant_transition_density` should by
default be `true`; can't imagine otherwise)
- Add response loss event (not debugged)
- Add birth event, along with response inheritance and vertical transmission
parameters (not debugged)

## 3 February 2025
Significant overhaul of weights calculations in order to accommodate correct
inter- and self-population contact

More than just this, but some notes I took:
- Change weights calculation for contacts at population level so that it
incorporates infected and total host numbers as well as inter-population
contacts
- Change weight propagation so that event weights are not propagated above
population level if the number of hosts is zero, and if they are zero, the
corresponding population weight is set to zero and propagated
- Implement contacts as inter- or intra-population contacts

## 31 January 2025
- Moved `hostContact!` down to section with events
- Renamed `MIGRATION` to `TRANSITION` and `migration_fractions` to
`transition_rates`
- Changed name of `PopulationParameters` to `PopulationType`
- `pathogenFractions!` is now `pathogenFractionsWinnerTakesAll`,
a parameter of `PopulationType`; default function is stored in
new file, `intrahost.jl`
- Made `weightedResponse` and `infectionProbability` into
`weightedResponseArithmeticMean` and `infectionProbabilityArithmeticMean`,
parameters of `PopulationType`; default functions are stored in
`immunity.jl`
- Moved `immunity.jl` and `intrahost.jl` into their corresponding
directories within the new `extensions` directory
- Add inoculum, mutation, recombination population-level coefficient
parameters
- Rename `population_contact_sum` to `contact_sum`

TODO: Implement inter-population contacts such that self-contact is an option
and we just have a single contacts function for all cases; receive contact
weights should be handled like contact weights, one per possible population
contact relationship


## 29 January 2025
- Removed `Class` and everything associated to it
- Made a single `CONTACT` event

## 23 January 2025
- Small bug fix in `newHost!`
- Track `intra_population_contact_sum` within `Population`
- Made intra-population contact rate at the `Population` level account for
receive contact and total host populations
- Changes in event modifier weights are now computed before events so that
receive weights are computed prior to contact rates
- Make changes in receive weights propagate into contact rates at the
`Population` level (different from what I described yesterday, but I like
this better)
- Fix weight propagation

#TODO: continue debugging contacts, not convinced simulation trajectories are
reasonable... actually, I think it works
If we can make inter-population contact rates depend on receive rates the way
intra-population contact works, we should consider two major changes in
architecture: (1) getting rid of `Classes` completely since `Populations` can
model anything a `Class` does, and (2) getting rid of the distinction between
intra- and inter-population contacts entirely, since an inter-population contact
with the same population achieves the same objective.

## 22 January 2025
- Make `Class`, `Population`, and `Model` receive weight matrices include weights
for all coefficients except for fitness
- Change `intraPopulationContact!` such that host 2 is sampled using receive
contact coefficients (instead of intrapopulation contact coefficients)
- Minor debug in `hostWeightsReceive!`, missing a couple lines
- Created `attemptInfection` function to handle checking for infection likelihod,
etc. before successfully adding a `Pathogen` to a `Host` so that it can be
called by `intraPopulationContact`, `interPopulationContact`, `establishMutant`,
`establishRecombinant`, and `birth` (if there is vertical transmission)
- Move `RECOMBINATION_PER_REPLICATION` and `MUTATION_PER_REPLICATION`
to named fields of `PathogenType`, add corresponding coefficient functions (and a
function for effective inoculum), make corresponding fields within `Pathogen`;
don't bother with the specific or static coefficient functions  in `Response`,
they don't matter; considered moving fitness to a separate field too but because
fitness is actually used in sampling, this is incorrect (rationale is the other
two coefficients are not used to sample entities)
- Changed computation of mutation and recombination probability during
transmission to be based on Poisson process means
- Moved contact algorithm into its own function to be called by
`intraPopulationContact` and `interPopulationContact` with a parameter specifying
whether or not to resample populations (#TODO: sampling of population 2 needs to be
done according to contact rates from population 1)
- Added propagation functions that go from the population level to the event level

KNOWN ISSUE: contact rate computation is incorrect; we have to multiply the sum of
contact rates by the sum of receive weights inside the `Class` level and divide by
the total number of hosts inside the `Population` level. This is a problem for the
current propagation architecture because it can handle changes in a single
variable at a time--to account for a change in contact weight and a change in
receive weight stemming from a single event, it currently has to go through
propagation cascades twice. To solve this, we will add a series of propagation
functions specifically for contact rates that take not only the change in contact
rate as a parameter, but also the simultaneous change in receive rate.

- Relatedly, we added `total_hosts` and `receive_contact_sum` fields to
`Population` to be used in the computation of contact rate and modified
`newHost!` to change both `total_hosts` and contact weights at the `Population`
level correspondingly

## 21 January 2025
- Some simple syntax debugging in events functions
- Important debugging in hostWeights computation, some weight matrix indexes were
missing and weight changes for propagation were not computed correctly for pathogen
and response events
- Make populations retain `RECEIVE_CONTACT` weights in `class_weights`
- Added `immunity.jl` to contain functions to simulate immune responses using
`Response` elements, added `immunityProbability` as a function that computes
probability that a `Host` is immune to a specific infecting `Pathogen`
- Added `immunityProbability` field to `ResponseType` that computes the
probability that a `Host` is immune to a specific infecting `Pathogen` thanks to
a `Response` of this type
- Changed `intraPopulationContact!` to check for immunity with `immunityProbability`
before adding a `Pathogen`
- Fix weight/rate computation to modify rate sum (wasn't happening)
- Replaced `sample()` with `chooseRandom()` in `intraPopulationContact`, Connor
benchmarked and said it was ~5X faster, at least sometimes?
- Renamed `immunityProbability` field within `ResponseType` as `infectionCoefficient`,
changed mechanism such that a coefficient of `1.0` is complete susceptibility and
`0.0` is complete immunity
- Added `reactivityCoefficient` as a field of `ResponseType` to handle how much a
`Response` contributes to a reaction against a specific pathogen
- Change sums of `Response` specific coefficients to weighted arithmetic means,
added `reactivityCoefficient` to calculation as weight, created a `weightedResponse`
function to compute all specific weighted mean responses (`infectionCoefficient`
weighted mean is computed separately); thought about whether these should be
weighted arithmetic or harmonic means, I think arithmetic (the work is being done
by the `reactivityCoefficient` function defining weights)
- Added logic to make recombination weights/probability zero when less than two
strains

Considered renaming `pathogen_fractions` to `pathogen_populations` and
`pathogenFractions` to `pathogenPopulations`, change the latter's behavior such
that the top pathogen's population is reported as its coefficient product, not just
`1.0`. Decided for sure against it; `pathogen_populations` captures competition,
actual pathogen-specific strength/influence on events is captured by other functions.

KNOWN ISSUES: infection events still don't seem to result in infection. This is a
problem for tomorrow. First time simulation runs though!

## 20 January 2025
- Rename `acquireResponse` to `developResponse` to avoid ambiguity about when a
`Response` is actually added to a `Host`
- Added a quasi zero-truncated Poisson distribution function as used in Python Opqua;
worth looking into whether an actual zero-truncated Poisson is preferrable
- Added `num_recombination_crossovers` and `mean_effective_inoculum` parameters to
`PathogenType`
(as an aside, I have decided to use Poisson distributions to model inoculum size
rather than negative binomials for simplicity and due to the fact that the increased
variability in mean inoculum size seen in Sobel Leonard _et al._ (2017) may have been
inflated due to the technical contaminations in Poon _et al._ (2016),
(https://doi.org/10.1128/jvi.00936-19)[https://doi.org/10.1128/jvi.00936-19])
- Change `pathogens` dictionary in `Class` to use genomes as keys and
correspondingly changed `responses` dictionary to use tuples of genomes as keys;
removed `id` field from both `Pathogen` and `Response` and `pathogens_idx` from
`Class`
- Rename `RECOMBINATION` event as `RECOMBINANT_ESTABLISHMENT`, rename
`MUTATION_AT_CONTACT` as `MUTATION_PER_REPLICATION`, add a new event choice
coefficient `RECOMBINATION_PER_REPLICATION`
- Added `binomial()` to generate binomial random numbers, faster than
`Distributions.jl`
- Added event functions for mutant establishment, pathogen clearance, response
acquisition, recombinant establishment, and intra-population contact--debugging
needed for these

## 19 January 2025
- Changed `pathogens` and `responses` vectors in `Class` to be dictionaries instead,
with keys corresponding to the assigned `Pathogen` integer code in the case of
`pathogens` and a tuple with the codes for imprinted and matured pathogens
in the case of `responses`
- Add `parent` as a new field of `Response`, storing a tuple referencing the
`Response` that this `Response` was derived from within its `Class`; this is only
useful for response lineage tracing, but not the simulation itself (at least at
the moment), so we may remove it
- Added base function to create a mutant `Pathogen`
- Added vector to track population receive weight sums in `Model`
- Changed weight calculation to update sum of weights in next level of hierarchy as
you iterate through the current level
- Added `Choice.jl` with functions to sample different entities based on their weights

## 18 January 2025
- Changed language such that all concepts relating to "immunity" now use "response"
- Renamed functions to add and remove pathogens and immunities from hosts
- Added new event, `RESPONSE_ACQUISITION`
- Removed `immunodominance` from `ResponseType`; this can be handled through
user-defined functions specifying how `Response` elements are acquired during
an `RESPONSE_ACQUISITION` event
- Added `PathogenType` to store `Pathogen`-specific parameters

## 17 Jan 2025
- Change `type` field within `Immunity` to be a reference to the `ImmunityType`,
not an index

## 12 Jan 2025
- Fixed many weight computation bugs

## 8 Jan 2025
- Added infection, immunization, clearance, deimmunization functions
- Changed all `Pathogen`, `Immunity`, and `ImmuneType` arrays to be arrays of
(references to) the entities, rather than arrays of their integer indexes
TODO: revise coefficient/weight/event nomenclature, clean up `weights.jl` in particular

## 5 Jan 2025
- Debugged `Host` weight initialization
TODO: Adding infection and immunization functions

## 4 Jan 2025
- Reorganized functions into new file structure
- Add and first debug struct creation methods
- Change `Pathogen` and `Immunity` entity storage within `Class` to use a simple vector
of entities, rather than a dictionary with their associated ID
- Change `Pathogen` and `Immunity` entity storage within `Host` to use vectors of entities
(references to them), rather than their indexes in another matrix
To do: Code to sample events and entities based on rates and weights

## 31 Dec 2024
Central change here is implementing a system for handling immunity through
`Immunity` entities that store references to the `Pathogen` sequences that were
imprinted and affinity matured to generate them and a reference to an `ImmunityType`
entity that defines how the immunity affects coefficients based on its imprinted
and matured sequences, as well as any infecting pathogens.
More changes:
- Moved struct definitions into `structs.jl` to compile before methods
- Changed `Pathogen` and `Immunity` vectors within `Hosts` to be vectors of keys
 used to access those objects within `Class` dictionaries instead
- Remove immune priority and the previous scheme used for weighting the effect of
different immunities; now, effect of immunities is not dependent on other pathogens
present
- Added `ImmunityType` struct to store "epitope types"; `Immunity` entities now
 have an associated `ImmunityType` that defines their coefficient functions and
 immunodominance; `ImmunityType` entities are stored inside `ClassParameters`
- Added `acquireImmunities` function field to `ClassParameters` to define which and
 how many `Immunity` entities (with corresponding `ImmunityTypes`) are acquired during
 immunization
- Changed `Immunity` to store the index of a `Pathogen` within the pathogen vector
in `Class` rather than storing a sequence
- Added coefficient vector to `Immunity`, make those coefficients affect weights
(these are coefficients that are "statically" active, i.e. they don't compare
sequences)
- Changed `Immunity` to store both imprinted and matured `Pathogen` indexes, then
change immunity functions such that they use up to three pathogen sequences (imprinted,
matured, and infecting) instead of only comparing two (immunized and infecting); in
this logic, imprinted immunity would have an imprinted pathogen sequence, but not a
matured one


## 16 Dec 2024
- Some quick bug fixes
- Added quick rundown of model structure to README
To do: `Immunity` is currently just a sequence, can be just a reference to the corresponding pathogen.

## 15 Dec 2024
- Made fitness into a pathogen-level coefficient
- Renamed "transmission" language to "contact"
- Distinguished between `INTRA_POPULATION_CONTACT` and `INTER_POPULATION_CONTACT`
- Changed transmission-based mutation into an event modifier coefficient rather
  than its own independent event

First draft of the following three:
- Added `Class` matrices for sampling pathogen, immunity, and host event weights according to event
- Added `Population` matrix for sampling class event weights according to event
- Added `Model` matrix for sampling population event weights according to event

Next: Code `Model` vector for sampling total event rates with correct rate computation,
add code for sampling class change and migration destination using corresponding rates in class.

## 13 Dec 2024
First draft of `Pathogen`, `Immunity`, and `Host` coefficients, fractions, and weights.
Next up are `Class` and real rates.

## 10 Dec 2024
First attempt at structs and event vector structure.

## 9 Dec 2024
Created directory and file structure.
