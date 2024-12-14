# using StaticArrays

# struct Model
#     id::Int64
#     populations::MVector{POPULATIONS,Population}
#     pathogen_rates::MVector{POPULATIONS,Float64}
#     immunity_rates::MVector{POPULATIONS,Float64}
#     host_rates::MVector{POPULATIONS,Float64}
# end

mutable struct Model
    id::Int64
    parameters::ModelParameters
    populations::Vector{Population} # size POPULATIONS
    pathogen_rates::Vector{Float64} # size POPULATIONS
    immunity_rates::Vector{Float64} # size POPULATIONS
    host_rates::Vector{Float64} # size POPULATIONS
end

function newModel(parameters::ModelParameters)

    # return Model(populations, pathogen_rates, immunity_rates, host_rates)
end
