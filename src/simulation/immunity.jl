using StaticArrays

struct Immunity
    id::Int64
    event_rates::SVector{NUM_IMMUNITY_EVENTS,Float64}
    coefficients::SVector{NUM_COEFFICIENTS,Int64}
    # Vector of IDs (in master user function dict) of functions that take the
    # sequences of the infecting pathogen and the immunized sequence
    # and return coefficients of the corresponding event
end
