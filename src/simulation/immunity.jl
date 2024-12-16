struct Immunity
    sequence::String
end

function newImmunity!(sequence::String, class::Class)
    immunity = Immunity(sequence)
    push!(class.immunities, immunity)
end
