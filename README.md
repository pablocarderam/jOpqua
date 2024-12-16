# jOpqua
Julia implementation of [Opqua](https://github.com/pablocarderam/opqua).
This is work in progress.

## Here's a rundown of the model.
The model is built on the following hierarchy of entities,
which are represented as structs:
`Pathogen` and `Immunity` << `Host` << `Class` << `Population` << `Model`

Each one has its own file. Besides the entities of the model, there's a list of
key processes that the simulation runs on. The likelihood of each process can
be modified by multiplication by coefficient of the corresponding type. Each
one of the entities in the model contains coefficients that modify some or all
of the processes in the model. The types of these coefficients (corresponding to
the processes) are defined as constants (`Int` numbers that are used to index into
matrices) in the file `constants.jl`. You'll see there are constants defined for
15 different coefficients. Each of these is a process within the simulation whose
likelihood can be affected by the genetic sequence of the pathogens infecting any
hosts that are involved, as well as any immunities that those hosts might be
carrying. The constants have been grouped according to what they represent and
what entities are involved. You'll see there are 10 constants referring to events
in the simulation (which will have their own independent rates) and 5 referring
to "choice modifiers", which are used to choose additional entities or events
that are needed when certain events happen.

The simulation works by randomly choosing which event type happens using
`event_rates` that are stored at the `Model` level and updated based on changes
to the model (currently working on that part). With an event chosen, now we need
to choose the specific entities that will be involved. To do this, each entity
within the hierarchy contains a matrix of "weights" with events on the rows and
sub-entities on the columns. These weights represent the relative likelihood of
the corresponding sub-entity participating in that particular event. For
instance, the `Model` entity has a matrix with 10 rows (one for each event) and
_p_ columns, where _p_ is the number of populations in the simulation. These are
used to randomly sample a population whenever a given event occurs using the
weights across the corresponding row. Similarly, `Population` entities have
equivalent matrices with _c_ columns, one per `Class`, and each `Class` has a
matrix with _h_ columns, one per `Host`. At the `Host` level, we have two
different matrices with different dimensions: one corresponding to events
involving `Pathogen` entities, and one corresponding to `Immunity` entities.
Some events involve `Host` entities but not specific `Pathogen` or `Immunity`
entities, those events are handled at the level of `Class` and are not
represented within `Host` entities.

Now, in addition to these weight matrices that allow us to sample entities from
their parents, we also have similar weight matrices that we use to sample
entities when they are "receiving" an event. Examples of this are contact
(formerly transmission, I think contact is more accurate depending on how
immunity is handled but I might revisit this nomenclature) events within a
`Population`, in which a `Pathogen` within a `Host` within a `Class` within a
`Population` (all sampled using their respective weight matrices along the
`INTRA_POPULATION_CONTACT` row) attempts to infect a different `Host` within
some `Class` within the same `Population` (the `Class` and `Host` sampled using
their respective `weight_receive` matrices along the `INTRA_POPULATION_CONTACT`
row).

At the bottom-most level, sampling of `Pathogen` entities within a `Host` is
computed based on the coefficients of that `Pathogen` genome, which are
functions of the genome sequence and are computed using user-defined functions
(in `pathogen_coefficient_functions`). These sequence-defined effects cascade
upwards through the hierarchy and influence the overall rates of events in the
simulation. `Immunity` also modifies `Pathogen` sampling at the bottom layer
through user-defined functions (`immunity_coefficient_effect_functions`) that
describe the effect of an existing `Immunity` sequence within a host on a given
`Pathogen` sequence for the corresponding coefficient type.

However, sampling of `Pathogen` entities also depends on the fraction of the total
intrahost population that each `Pathogen` occupies. Calculating this fraction is
tricky because it's not static over time. In the original Opqua algorithm, the
fraction was proportional to intrahost fitness. In the way I'm currently
implementing jOpqua, it's winner-takes-all according to intrahost fitness and the
effect of `Immunity` entities on that fitness—probably most accurate for acute
infections with strong population bottlenecks, like flu. In the future version of
Opqua I keep talking about, we'll use population genetics approximations to track
changes in intrahost populations. An exception to sampling `Pathogen` entities using
population fraction is the event of clearance (formerly recovery): the likelihood
of pathogen clearance is not necessarily proportional to the pathogen's population
fraction, in fact it may be the opposite. Since it'll depend on many things, I
leave it as equal weighting for all pathogens, regardless of representation.

Similarly, each `Immunity`'s contribution to `Pathogen` sampling weights depends
on the fraction of activity of that `Immunity` within the `Host`. That fraction is
computed based on each `Immunity`'s immunodominance within the host, which is a
function of both the `Immunity` sequence and the `Pathogen` entities present at
that given time within the `Host` (but not their population fractions, otherwise
we have a circular definition).

The only `Immunity` event currently in the simulation is loss of immunity. The
likelihood of `Immunity` events depends on each immunity's sequence, as defined in
`immunity_coefficient_functions`. Immunodominance doesn't increase the likelihood
of loss (probably decreases it, but we'll leave it as equal weighting for all
immunities, as done with pathogen clearance).

The coefficient modifier functions that operate on `Pathogen` and `Immunity`
sequences along with all other static base parameters describing event rates and
likelihoods are defined at the level of the `Class`. This means that `Class`
weights within a `Population` not only are weighted according to the sum of `Host`
weights, but are also multiplied by `Class`-specific rates.

This model structure is pretty flexible: we can accommodate metapopulation
structures, host compartments with different epidemiological parameters (not
possible in original Opqua), vector-borne dynamics (using different populations
for hosts and vectors, more elegantly done than original Opqua), and flexible and
powerful effects of immunity (much more so than in the original Opqua). I'm still
not sure how to go about distinguishing imprinted vs. affinity matured immunity
using this structure—I think there's a way to code it entirely through the
user-defined functions, or maybe just with some light changes. But I'll worry
about that later.

I haven't started debugging so it's probably full of issues. I'm waiting to have
a working prototype.
