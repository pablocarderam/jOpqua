struct Model
    id::Int64
    populations::Vector{POPULATIONS,Population}
    pathogen_rates::Vector{POPULATIONS,Float64}
    immunity_rates::Vector{POPULATIONS,Float64}
    host_rates::Vector{POPULATIONS,Float64}
end
