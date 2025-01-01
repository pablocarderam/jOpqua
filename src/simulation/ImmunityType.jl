function newImmunityType!(
    id::String,
    static_coefficient_functions::SVector{NUM_COEFFICIENTS,Function},
    specific_coefficient_functions::SVector{NUM_COEFFICIENTS,Function},
    immunodominance::Function,
    class::Class)
    class.parameters.immunity_types[id] = ImmunityType(
        id, static_coefficient_functions, specific_coefficient_functions, immunodominance,
    )
end
