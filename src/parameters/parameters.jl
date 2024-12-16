using StaticArrays

struct ClassParameters
    id::String

    base_coefficients::SVector{NUM_COEFFICIENTS,Float64}

    class_change_fractions::Dict{String,Float64} # size CLASSES, must sum to 1

    inter_population_contact_fractions::Dict{String,Float64} # size POPULATIONS, must sum to 1
    migration_fractions::Dict{String,Float64} # size POPULATIONS, must sum to 1

    pathogen_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # Each element takes seq argument, returns Float64
    immunity_coefficient_functions::SVector{NUM_COEFFICIENTS,Function} # Each element takes seq argument, returns Float64

    immunityDominance::Function # Takes two arguments (vector of all pathogen sequences in host, immunity sequence), returns Float64
    immunity_coefficient_effect_functions::SVector{NUM_COEFFICIENTS,Function} # Each element takes two seq arguments (infecting pathogen, immunity), returns Float64
end
