using StaticArrays

struct Immunity
    sequence::String
    # priority::Float64
    # event_rates::SVector{NUM_IMMUNITY_EVENTS,Float64}
    # coefficient_functions::SVector{NUM_COEFFICIENTS,Int64}
    # # Vector of IDs (in master user function dict) of functions that take the
    # # sequences of the infecting pathogen and the immunized sequence
    # # and return coefficients of the corresponding event
end

function newImmunity!(class::Class, sequence::String)
    immunity = Immunity(
        sequence,
        # sequencePriority(class, sequence)
        # getEventRatesImmunity(class, sequence),
    )
    push!(class.immunities, immunity)
end
