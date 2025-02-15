using Printf

# Flexible binary level rejection sampling

# Structs

mutable struct FlexLevel
    bounds::Tuple{Float64,Float64}
    sum::Float64
    max::Float64
    indices::Vector{Int64}
    next::Union{FlexLevel,Nothing}
end

mutable struct FlexlevSampler
    min_level::Float64
    max_level::Float64
    levels::FlexLevel
    weights::AbstractVector{Float64}
end

# Methods

# Initialization

function newFlexlevSampler(weights::AbstractVector{Float64})
    if length(weights) == 0
        throw("Cannot create FlexlevSampler from AbstractVector of length 0.")
    end

    head = newFlexLevel(1, weights[1])
    min, max = head.bounds
    for i in firstindex(weights)+1:lastindex(weights)
        w = weights[i]

        # if need new level at beginning
        if w >= head.bounds[2]
            new_head = newFlexLevel(i, w, next=head)
            head = new_head
            max = head.bounds[2]
            continue
        end

        curr::Union{FlexLevel,Nothing} = head
        while true
            if curr.bounds[1] <= w < curr.bounds[2]     # belongs in current level
                addToFlexLevel!(i, w, curr)
                break
            elseif isnothing(curr.next)                 # need new level at end
                curr.next = newFlexLevel(i, w)
                min = curr.next.bounds[1]
                break
            elseif w >= curr.next.bounds[2]             # need new level between curr and next
                curr.next = newFlexLevel(i, w, next=curr.next)
                break
            end
            curr = curr.next
        end
    end

    return FlexlevSampler(min, max, head, weights)
end

function newFlexLevel(i::Int64, w::Float64; next::Union{FlexLevel,Nothing}=nothing)
    return FlexLevel(logBounds(w), w, w, [i], next)
end

function addToFlexLevel!(i::Int64, w::Float64, level::FlexLevel)
    push!(level.indices, i)
    level.sum += w
    if w > level.max
        level.max = w
    end
end

# Testing

function verifyFlexlevSampler(sampler::FlexlevSampler)
    errors = 0

    # ]all weights in sampler.weights are represented in some level?
    for i in eachindex(sampler.weights)
        if !inSampler(i, sampler)
            errors += 1
            @printf "Error %i: index %i (weight %f) not in any level\n" errors i sampler.weights[i]
        end
    end

    # all weights in levels are represented in sampler.weights, and if so, in the correct level?
    curr = sampler.levels
    num_weights = length(sampler.weights)
    while !isnothing(curr)
        for i in curr.indices
            if !(1 <= i <= num_weights)
                errors += 1
                @printf "Error %i: index %i present in level (%f, %f) but not in weights\n" errors i curr.bounds[1] curr.bounds[2]
            elseif !(curr.bounds[1] <= sampler.weights[i] < curr.bounds[2])
                errors += 1
                @printf "Error %i: index %i (weight %f) present in level (%f, %f)\n" errors i sampler.weights[i] curr.bounds[1] curr.bounds[2]
            end
        end
        curr = curr.next
    end

    # all level sums/maxes correct?
    curr = sampler.levels
    while !isnothing(curr)
        expected_sum = sum(sampler.weights[i] for i in curr.indices)
        if expected_sum != curr.sum
            errors += 1
            @printf "Error %i: level (%f, %f) sum incorrect (expected %f, got %f)\n" errors curr.bounds[1] curr.bounds[2] expected_sum curr.sum
        end
        expected_max = maximum(sampler.weights[i] for i in curr.indices)
        if expected_max != curr.max
            errors += 1
            @printf "Error %i: level (%f, %f) max incorrect (expected %f, got %f)\n" errors curr.bounds[1] curr.bounds[2] expected_max curr.max
        end
        curr = curr.next
    end

    @printf"Final error count: %i" errors
    return errors
end

function inSampler(i::Int64, sampler::FlexlevSampler)
    curr = sampler.levels
    while !isnothing(curr)
        if i in curr.indices
            return true
        end
        curr = curr.next
    end
    return false
end

function printFlexlevSampler(sampler::FlexlevSampler)
    @printf "\nFlexlevSampler (%f, %f)\n" sampler.min_level sampler.max_level
    show(sampler.weights)
    println()
    curr = sampler.levels
    while !isnothing(curr)
        @printf "\nLevel (%f, %f)\n" curr.bounds[1] curr.bounds[2]
        @printf "sum=%f\n" curr.sum
        @printf "max=%f\n" curr.max
        println(curr.indices)
        curr = curr.next
    end
end
