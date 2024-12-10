struct Class
    id::Int64
    hosts::Vector{MAX_HOSTS,Host}
    pathogen_rates::Vector{MAX_HOSTS,Float64}
    immunity_rates::Vector{MAX_HOSTS,Float64}
    host_rates::Vector{MAX_HOSTS,Float64}
end
