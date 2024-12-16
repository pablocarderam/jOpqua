# jOpqua Changelog

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
