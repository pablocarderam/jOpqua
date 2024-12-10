using StaticArrays

struct Pathogen
    sequence::String
    rates::SVector{NUM_PATHOGEN_EVENTS,Float64}
    coefficients::SVector{NUM_COEFFICIENTS,Float64}
end
