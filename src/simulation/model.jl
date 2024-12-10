struct Model
    id::Int64
    populations::Vector{POPULATIONS,Population}
    rates::Vector{POPULATIONS,Float64}
end
