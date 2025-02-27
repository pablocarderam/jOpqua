# jOpqua Changelog

TODO:
- Add response acquisition upon clearance as preferred alternative to response
acquisition during infection (as in mutations upon infection vs. mutation
establishment)
- Add option for Newick to show real time rather than Hamming distance

TODO: Not debugged:
- Mutant establishment
- Recombinant establishment
- Recombination upon contact
- Inter-population contact
- Birth
- Death
- Transition

## 26 February 2025
- Fixed bug where changes in contact receive weights were not correctly
propagating into contact rates (missing base population coefficient
multiplication) (PCR)
- Fixed floating point error that leads to very small but nonzero contact rate
(when it should be zero) (PCR)
- Fixed response acquisition bug (wrong constant) (PCR)
- Fixed response loss bug (wrong constant) (PCR)
- Fixed bugs in winner takes all immunity functions, which weren't taking the
`Response` with highest reactivity (PCR)
- Added `reactivityCoefficient`, `infectionCoefficient`, and
`responseStaticCoefficient` functions to streamline immunity and weight
functions (PCR)
- Set error tolerance constant to 1e-9 (PCR)
- Added error tolerance checking at model rate level (PCR)
- Moved `generateGamete` and `generateZygote` functions up in file (PCR)

## 25 February 2025
- Removed intrahost fitness coefficient from `PopulationType` since it was
never used (PCR)
- Changed default parameters to be all zeros for base coefficient/rates and
all maps to one for coefficient functions; Booleans were left as the most
common usage case (PCR)
- Reorganized examples into folders (PCR)
- Separated initializer functions for Type structs into a new file, `setup.jl`
- Renamed `parameters.jl` to `default.jl` (PCR)
- Fixed bug renormalizing the same random number with recombination events
during transmission (PCR)
- Changed structure of `developImmunity` arguments (PCR)
- Added logic to verify whether `Response` entities have actual `Pathogen`
entities associated, or just `nothing` (PCR)
- Added `deNovoResponse` as default `developImmunity` function (PCR)

KNOWN ISSUE: runtime didn't finish when I tried with
`response_acquisition_coefficient=1.0,`!!!

## 24 February 2025
- Added proportional fitness function for pathogen fractions (PCR)
- Debugged problems with complex weight calculation (contact, transition,
receive contact, receive transition), particularly (1) when adding hosts in
bulk, (2) when removing hosts, (3) in the main `hostWeights!` function (PCR)
- Write test function for `flexle.jl` user functions (CLM)
- Fix `flexle.jl` bugs mostly involving improper handling of zero-weight entries,
plus undeclared variable in `addToFlexleSampler!` as noted by PCR (CLM)

## 23 February 2025
- Remove function that plots compartments from output object to standardize
all plotting functions as working with `DataFrames` (PCR)
- Fix bugs in `compartmentPlot` function so that it works with `DataFrame`
input (PCR)
- Add total calculation to `saveComposition` (PCR)
- Make composition plots (PCR)
- Make Newick format phylogenies of pathogens, hosts, and responses (PCR)
- Plot phylogenies using the `NewickTree.jl` package--might be easier to plot
the newicks online though (PCR)
- Fix variable type bug in `zeroTruncatedPoisson` (PCR)
- Fix missing `transition!` function in event function constant vector (PCR)
- Change weight calculation at top level to minimize imprecision (overflow?)
errors (PCR)
- Removed `FunctionWrapper` type definitions in arguments of setup functions
because they can't handle passing functions as arguments apparently... hopefully
the type definitions at the struct level keep performance okay; the flamegraph
still doesn't show any non-compiled code at least (PCR)

## 22 February 2025
- Moved `logbounds` to `flexle.jl` (PCR)
- Type-defined `zeroTruncatedPoisson` and `catcol` (PCR)
- Fixed small bug in host weight calculation from adding host genomes (PCR)
- Separated out `mutateSequence` function (PCR)
- Add sexual reproduction option for `birth` (PCR)
- Added parameters for mutation and recombination of `Host` genome sequences
(PCR)
- Changed `mutateSequence` to account for chromosome separators (PCR)

Note: I think flexle's `addToFlexleSampler!` has an undeclared variable `i`

## 21 February 2025
- decided on `flexle.jl` API, functions written but not yet fully tested (CLM)
- Added genome sequence to `Host` (PCR)
- Made `Response` entities store the `Host` genome sequence they were derived in
(PCR)
- Made functions in `ResponseType` take host genome sequences into account (PCR)

## 20 February 2025
- `flexle.jl` additions including sampling; updating sampler `weights`; and
sample testing (CLM)
- Removed commented code for harmonic means and alternate algorithms in
immunity (PCR)
- Small bug fix in `developResponse` to pass correct arguments (PCR)
- Added a commented out `birth!` function to be used as a starting point for
sexual reproduction
- Separated out `recombineSequences` function (PCR)
- Changed `weightedResponse` and `infectionProbability` defaults to not be
the arithmetic means, but rather winner takes all style functions (PCR)

I decided that these `weightedResponse` and `infectionProbability` functions
must be imprinting-like functions. The simplest is winner takes all, but other
options could include top N respondents or some sort of a competition-based
algorithm. This could be something like this: you have a set amount of
resource for your total response (or response to a specific epitope,
corresponding to a `ResponseType`) and you allocate that resource to each
`Response` (of that `ResponseType`) starting from the one with the strongest
`weightedResponse` and going as far down the list as you still have resource
left to allocate. I don't love that specific idea though, not very mechanistic.
Regardless, as a default, it's best to do winner takes all, just like we're doing
with the default for intrahost pathogen fractions. For the "advanced" immunity
expansion/DLC, we might have different classes of `ResponseTypes` corresponding
to naive B cells, memory B cells, and long-lived plasmablasts, and at that
point we can rethink the best algorithm to use.

Based on intuition from this paper Amanda Elyssa Ruiz posted on Slack:
[10.1016/j.celrep.2025.115334](https://doi.org/10.1016/j.celrep.2025.115334)
which (copying from my Slack message) "vibes with the point Rachel \[Hecht]
raised during lab meeting that (if I'm remembering correctly) during
flu imprinting, when people look at GCs, they see a lot of B cell diversity,
but when they look at serum antibodies, they just see imprinted responses.
It vibes in the sense that this also suggests imprinting doesn't occur at
the level of selecting B cells that undergo maturation, but at the level of
selecting B cells to produce plasmablasts? If I understand things correctly?"

TODO:
- Add `Host` genome sequences, pass them to `developResponse`
- Add `Host` sexual reproduction, with second parent selection and
recombination during `Host` birth events

## 19 February 2025
- Finished composition dataframe (PCR)
- Bug fix on history dataframe saving (PCR)

## 18 February 2025
- Added data csv export of compartments data (PCR)
- Added csv plotting of compartments data (PCR)
- Added option of `nothing` for `Response.matured.pathogen` (PCR)
- Added data csv export of full data, can't believe it worked easily (PCR)
- Change `Host` sampling for history to a random sample of hosts fixed over
time (PCR)
- Started work on composition dataframe creation (PCR)

Next: composition dataframe and plot, then immunity

## 17 February 2025
- Small bug fix: `mean_recombination_crossovers` type corrected (PCR)
- Changed argument order in Type structs to put coefficients last within
their section (PCR)
- Finish model setup parameter initialization functions, default values (PCR)
- Extend `FlexleSampler` levels to include higher and lower bounds (CLM)
- `FlexleSampler` testing (CLM)
- Replace abstract `Function` type definitions in structs and function
arguments with `FunctionWrapper` functions with defined outputs and inputs (PCR)
- Type-define default argument functions (PCR)
- Identified performance issue using flamegraph in which mutable static vector
defining `inocula` within `hostContact!` resulted in runtime function resolution;
made `inocula` into a regular `Vector` and this removed the runtime resolution
of the calls (but addded orange bars on the flamegraph in what I'd consider
unrelated parts of the code, the return lines of `choosePathogen!` and
`chooseHost!`--worth keeping an eye on, maybe with `inocula` of length > 1 we
will be able to spot performance problems and `MVector` will be better) (PCR)

At this point, the flamegraph shows no runtime (red) calls!

## 16 February 2025
- Bug fix in compartment variable handling (PCR)
- New `Output` struct for simulation output (PCR)
- Population compartment plot (PCR)
- Moved sandbox file to examples folder (PCR)
- Default model setup functions (PCR)
- Add/remove/move indices to/from/between `FlexleSampler` levels (CLM)

Next:
- Finish model setup functions
- Immunity draft

## 15 February 2025
- Moved FlexLev/Flexle code into a new file for modularization of the technique
(PCR)
- Added `StaticHost` for more lightweight storage of host info in simulation
history (PCR)
- Added history tracking in two ways: accurate counts of compartment variables
(uninfected naive, infected naive, uninfected immune, infected immune, and dead
hosts) in each population at each indicated time point, and a sample of hosts
(each stored as `StaticHost`) from each population at each time point (PCR)
- Reworked `Flexle` sampling to use arrays as described below (CLM)

## 14 February 2025
- Added references to actual parent entities in `Pathogen` and `Response` using
CLM's `Union{X, Nothing}` trick (PCR)
- Implemented flexible level + rejection sampling (current working titles "FlexLev
sampling", or perhaps "Flexle sampling") as described in PCR's note yesterday.
However, in discussions today re: the relative benefits of using a linked list
structure for `FlexLevel` entities in a `FlexlevSampler` as opposed to an array,
we decided that the primary LL advantage—namely the ability to easily add new
`FlexLevel` without copying data—is likely to be outweighed by both (1) that arrays
allow for constant-time indexing, meaning adding, removing, and moving weights
to/from/between levels (actions which are all likely to be far more common than
creating new levels) are much faster for arrays than for linked lists, and (2) that
the "add level" operation is actually still _O(n)_ in the number of levels for linked
lists when levels are added in the middle of the list, thereby partially negating the
linked list's primary benefit. (As a side note, arrays actually lose their constant-time
indexing property for our use case if you allow for non-contiguous `FlexLevel`s in a
`FlexlevSampler`. However, we have decided for now at least that the additional
space/time to enforce level contiguousness is an appropriate price to pay for
constant-time indexing.) All of that is to say that
*flexle sampling as currently implemented using linked lists will not remain*.
CLM will rework the implementation to use arrays instead. Nevertheless, this iteration
of flexle sampling is being committed for posterity. (CLM)

## 13 February 2025
No changes yet, but based on (unnecessarily, sorry) lengthy discussion, we conclude
that the optimal implementation of random weighted sampling of `Host` entities is
a modified version of the binary level rejection sampling technique described by
[Aaron DeFazio](https://www.aarondefazio.com/tangentially/?p=58) based on work
from Donald Knuth and others ("The technique used is not novel, indeed it is
based on publications from the 1960s" says DeFazio--but who?!).

We will implement rejection sampling with binary levels as described by DeFazio &
1960s, including optimizations such as random number renormalizing as done in Opqua
already and also explained by DeFazio in the "Enhancements" section. We discussed
whether we should include in the algorithm a check on whether to add a new level
at the top or bottom of the vector of levels, in case the ratio between the largest
and smallest weights exceeded the 2^20-fold recommended by DeFazio (who argues for
20 intervals, given that the dynamic range of weight values is usually less than
1 million, but nothing guarantees this). We decided to do this, and we might as
well also remove unnecessary levels at the top or bottom of the level list if they
become empty of weights. In order to most efficiently add and remove levels both at
the start and end of the list of levels, CLM suggested we use a linked list. We
also discussed whether rejection sampling of the list of levels itself would
outperform linear searching our way through it, as implemented by DeFazio, but we're
unsure, given the potentially lopsided differences in the sums of weights between
different levels. We decided we will stick with linear search across levels.

We are pretty confident this is a good idea for `Host` sampling given the number of
hosts in a `Population`. We suspect the current linear search algorithm implemented in
`randChoose` is optimal for sampling event types, given there are only 9 events:
according to tests done by CLM in Julia, a uniformly-balanced weights vector is
sampled faster by linear search than by simple rejection sampling if there are fewer
than 10 weights being sampled* (PCR is now unsure of this, see note below). We are
unsure whether sampling of `Pathogen` and `Response` entities within a `Host` is best
done with linear search or rejection sampling given that the distribution of the
number of pathogens and responses is both variable and skewed in favor of smaller
numbers. Because of the latter, we will remain on linear search for the time being.

*Bear in mind that sampling from a uniform weight array represents the best possible
scenario for simple rejection sampling, since the algorithm reduces to indexing into
a particular position in the vector based on a random integer and then always
accepting the given element regardless of the second random number that accepts or
rejects the chosen weight/element. However, I now remember that this second random
number can be obtained by renormalization instead of fresh regeneration. CLM did not
implement this detail in the test, so maybe the vector size at which rejection
sampling becomes faster than linear search is closer to 5 than 10? This makes me less
confident that linear search is preferrable now; perhaps testing binary level
rejection sampling for events and even pathogens and responses is actually worth it.

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
