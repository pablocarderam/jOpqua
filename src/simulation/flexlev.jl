using Printf

# Flexible binary level rejection sampling

# Structs

mutable struct FlexLevel
    bounds::Tuple{Float64,Float64}
    sum::Float64
    max::Float64
    indices::Vector{Int64}
end

mutable struct FlexlevSampler
    levels::Vector{FlexLevel}
    weights::AbstractVector{Float64}
    min::Float64
    max::Float64
    sum::Float64
end

# Methods

"""
    levelIndex(w, u)

Given a weight `w`, returns the index of the level in some `FlexlevSampler.levels` with maximum upper bound `2^u` where `w` belongs.

Returns `0` if `w` is `0.0`, indicating that `w` belongs in no level.

# Examples

`levelIndex(14.2, 6)` ==>  `3`

The value `14.2` belongs in the `(8.0, 16.0)` level, which in a `FlexlevSampler` that
starts with a level of bounds `(32.0, 64.0)` (`64` being `2^6`) is at `levels[3]`.

`levelIndex(8.0, 4)` ==> `1`

`levelIndex(8.0, 5)` ==> `2`
"""
function levelIndex(w::Float64, u::Int64)
    return iszero(w) ? 0 : u - Int64(floor(log2(w)))
end

# Initialization

function newFlexlevSampler(weights::AbstractVector{Float64})
    if length(weights) == 0
        throw("Cannot create FlexlevSampler from AbstractVector of length 0.")
    end

    w_sum = 0.0
    weights_nonzero = filter(x->!iszero(x), weights)
    all_zero = isempty(weights_nonzero)
    
    if all_zero
        w_min, w_max = 0.0, 0.0
        num_levels = 0
    else
        w_min, w_max = Inf, 0.0
        for w in weights_nonzero
            (w < w_min) && (w_min = w)
            (w > w_max) && (w_max = w)
        end

        max_u_log   = Int64(ceil(log2(w_max)))
        max_l_log   = Int64(floor(log2(w_max)))

        min_l_log_f = floor(log2(w_min))
        min_l_log   = Int64(min_l_log_f)

        l_bound     = 2.0^min_l_log_f
        num_levels = max_l_log - min_l_log + 1     # e.g. -2,4 ==> 7 levels [4,3,2,1,0,-1,-2]
    end

    levels = Vector{FlexLevel}(undef, num_levels)   # add check for unreasonable number of levels before allocating space?
    
    if !all_zero
        for i in num_levels:-1:1
            u_bound = l_bound * 2.0
            levels[i] = FlexLevel((l_bound, u_bound), 0.0, 0.0, Vector{Int64}()) 
            l_bound = u_bound
        end

        for i in eachindex(weights)
            w = weights[i]
            if !iszero(w)
                l = levels[levelIndex(w, max_u_log)]
                push!(l.indices, i)
                (w > l.max) && (l.max = w)
                l.sum += w
                w_sum += w
            end
        end
    end

    return FlexlevSampler(levels, weights, w_min, w_max, w_sum)
end

function newFlexLevel(i::Int64, w::Float64)
    return FlexLevel(logBounds(w), w, w, [i])
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

    # all weights in sampler.weights are represented in some level?
    for i in eachindex(sampler.weights)
        if !iszero(sampler.weights[i]) && !inSampler(i, sampler)
            errors += 1
            @printf "Error %i: index %i (weight %f) not in any level\n" errors i sampler.weights[i]
        end
    end

    # all weights in levels are represented in sampler.weights, and if so, in the correct level?
    num_weights = length(sampler.weights)
    for level in sampler.levels
        for i in level.indices
            if !(1 <= i <= num_weights)
                errors += 1
                @printf "Error %i: index %i present in level (%f, %f) but not in weights\n" errors i level.bounds[1] level.bounds[2]
            elseif !(level.bounds[1] <= sampler.weights[i] < level.bounds[2])
                errors += 1
                @printf "Error %i: index %i (weight %f) present in level (%f, %f)\n" errors i sampler.weights[i] level.bounds[1] level.bounds[2]
            end
        end
    end

    # all level sums/maxes correct?
    for level in sampler.levels
        empty = isempty(level.indices)
        expected_sum = empty ? 0.0 : sum(sampler.weights[i] for i in level.indices)
        if expected_sum != level.sum
            errors += 1
            @printf "Error %i: level (%f, %f) sum incorrect (expected %f, got %f)\n" errors level.bounds[1] level.bounds[2] expected_sum level.sum
        end
        expected_max = empty ? 0.0 : maximum(sampler.weights[i] for i in level.indices)
        if expected_max != level.max
            errors += 1
            @printf "Error %i: level (%f, %f) max incorrect (expected %f, got %f)\n" errors level.bounds[1] level.bounds[2] expected_max level.max
        end
    end

    @printf"Final error count: %i" errors
    return errors
end

function inSampler(i::Int64, sampler::FlexlevSampler)
    for level in sampler.levels
        if i in level.indices
            return true
        end
    end
    return false
end

function printFlexlevSampler(sampler::FlexlevSampler)
    @printf "\nFlexlevSampler (%f, %f)\n" sampler.min sampler.max
    show(sampler.weights)
    @printf "\nSum %f" sampler.sum
    println()
    if isempty(sampler.levels)
        println("(sampler contains no levels)")
    else
        for l in sampler.levels
            @printf "\nLevel (%f, %f)\n" l.bounds[1] l.bounds[2]
            @printf "sum=%f\n" l.sum
            @printf "max=%f\n" l.max
            println(l.indices)
        end
    end
end
