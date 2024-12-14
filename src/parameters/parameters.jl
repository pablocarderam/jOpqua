using StaticArrays

struct ModelParameters
    id::String

    num_populations::Int64
end

struct PopulationParameters
    id::String

    num_classes::Int64
end

struct ClassParameters
    id::String

    base_coefficients::SVector{NUM_COEFFICIENTS,Float64}

    pathogenFitness::Function
    pathogen_coefficient_functions::SVector{NUM_COEFFICIENTS,Function}
    immunityDominance::Function
    immunity_coefficient_functions::SVector{NUM_COEFFICIENTS,Function}
    immunityEffectOnFitness::Function
    immunity_coefficient_effect_functions::SVector{NUM_COEFFICIENTS,Function}
end
