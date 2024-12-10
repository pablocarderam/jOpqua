struct Population
    id::Int64
    classes::Vector{CLASSES,Class}
    pathogen_rates::Vector{CLASSES,Float64}
    immunity_rates::Vector{CLASSES,Float64}
    host_rates::Vector{CLASSES,Float64}
end
