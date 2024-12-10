struct Class
    id::Int64
    hosts::Vector{MAX_HOSTS,Host}
    rates::Vector{MAX_HOSTS,Float64}
end
