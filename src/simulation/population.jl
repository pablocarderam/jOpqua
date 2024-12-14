mutable struct Population
    id::String
    parameters::PopulationParameters
    classes::Vector{Class} # size CLASSES
    pathogen_rates::Vector{Float64} # size CLASSES
    immunity_rates::Vector{Float64} # size CLASSES
    host_rates::Vector{Float64} # size CLASSES
end

# function id(c::Population)
#     return (c.population_id)
# end
