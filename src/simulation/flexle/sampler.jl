# Structs and maintenance functions internal to the module

const EXPONENT_MASK_FLOAT64::Int64   = 0x7FF0000000000000
const EXPONENT_SHIFT_FLOAT64::Int64  = 52
const EXPONENT_OFFSET_FLOAT64::Int64 = 1023

# Structs

mutable struct FlexLevel
    bounds::Tuple{Float64,Float64}
    sum::Float64
    max::Float64
    indices::Vector{Int64}
end

mutable struct FlexleSampler
    levels::Vector{FlexLevel}
    weights::Vector{Float64}
    sum::Float64
    index_positions::Vector{Int64}
end

# Initialization

function FlexLevel(i::Int64, w::Float64)
    return FlexLevel(logBounds(w), w, w, [i])
end

# Utility methods

import Base.&

Base.:&(f::Float64, i::Int64) = reinterpret(Float64, (reinterpret(Int64, f) & i))

"""
    approxeq(a, b, t=1e-9)

Return a `Bool` indicating whether `a` and `b` are equal to within `t`.
"""
function approxeq(a::Float64, b::Float64; t::Float64=1e-9)
    return abs(a-b) < t
end

"""
    lowerLog2Bound(n)

Get the largest power of 2 less than or equal to a non-negative `n`.

# Examples

`lowerLog2Bound(33.0)` ==> `32.0`

`lowerLog2Bound(32.0)` ==> `32.0`

`lowerLog2Bound(0.75)` ==> `0.5`
"""
function lowerLog2Bound(n::Float64)
    return n & EXPONENT_MASK_FLOAT64
end

"""
    floorLog2(n)

Get the `Int64` floor of the log2 of `n`.

# Examples

`floorLog2(33.0)` ==> `5`

`floorLog2(32.0)` ==> `5`

`floorLog2(0.75)` ==> `-1`
"""
function floorLog2(n::Float64)
    return ((reinterpret(Int64, n) & EXPONENT_MASK_FLOAT64) >> EXPONENT_SHIFT_FLOAT64) - EXPONENT_OFFSET_FLOAT64
end

"""
    logBounds(n)

Return a tuple `l,u` giving two adjacent powers of 2 such that `l <= n < u`.
"""
function logBounds(n::Float64)
    l = lowerLog2Bound(n)
    return l, l*2.0
end

"""
    levelIndex(w, u)

Given a weight `w`, return the index of the level in some `FlexleSampler.levels` with maximum upper bound `2^u` where `w` would belong.

Returns `0` if `w` is `0.0`, indicating that `w` belongs in no level.

# Examples

`levelIndex(14.2, 6)` ==>  `3`

The value `14.2` belongs in the `(8.0, 16.0)` level, which in a `FlexleSampler` that
starts with a level of bounds `(32.0, 64.0)` (`64` being `2^6`) is at `levels[3]`.

`levelIndex(8.0, 4)` ==> `1`

`levelIndex(8.0, 5)` ==> `2`
"""
function levelIndex(w::Float64, u::Int64)
    return iszero(w) ? 0 : u - floorLog2(w)
end

"""
    levelIndex(w, levels)

Get the index of the `FlexLevel` in `levels` where a weight `w` would belong.

Returns `0` if no such level exists.
"""
function levelIndex(w::Float64, levels::Vector{FlexLevel})
    isempty(levels) && return 0
    idx = floorLog2(levels[1].bounds[2]) - floorLog2(w)
    return idx > length(levels) ? 0 : idx
end

"""
    levelIndex(bounds, levels)

Get the index of the `FlexLevel` in `levels` with bounds given by `bounds`.

Returns `0` if no such level exists.
"""
function levelIndex(bounds::Tuple{Float64,Float64}, levels::Vector{FlexLevel})
    return levelIndex(bounds[1], levels)
end

"""
    getLevel(bounds, levels)

Return the `FlexLevel` in `levels` with bounds given by `bounds`.

Returns `nothing` if no such level exists.
"""
function getLevel(bounds::Tuple{Float64,Float64}, levels::Vector{FlexLevel})
    l = levelIndex(bounds, levels)
    return l==0 ? nothing : levels[l]
end

"""
    getLevel(w, levels)

Return the `FlexLevel` in `levels` where a weight `w` would belong.

Returns `nothing` if no such level exists.
"""
function getLevel(w::Float64, levels::Vector{FlexLevel})
    l = levelIndex(w, levels)
    return l==0 ? nothing : levels[l]
end

"""
    logDist(a, b)

Return the floor of the log2 of `a/b`. 
"""
function logDist(a::Float64, b::Float64)
    return floorLog2(b) - floorLog2(a)
end

"""
    inSampler(i, sampler)

Return a `Bool` indicating whether an index `i` is present at some `FlexLevel` in `sampler`.
"""
function inSampler(i::Int64, sampler::FlexleSampler)
    for level in sampler.levels
        if i in level.indices
            return true
        end
    end
    return false
end

"""
    inSampler(bounds, sampler)

Return a `Bool` indicating whether a `FlexLevel` with bounds `bounds` is present in `sampler`.
"""
function inSampler(bounds::Tuple{Float64,Float64}, sampler::FlexleSampler)
    return (sampler.levels[begin].bounds[1] >= bounds[1]) && (bounds[1] >= sampler.levels[end].bounds[1])     # bounds between largest and smallest levels' bounds (inclusive)
end

"""
    levelIsPopulated(level)

Return a `Bool` indicating whether `level` contains any elements.
"""
function levelIsPopulated(level::FlexLevel)
    return !isempty(level.indices)
end

"""
    maxLevelWeight(level, sampler)

Get the largest weight in `sampler` of any element in `level`.
"""
function maxLevelWeight(level::FlexLevel, sampler::FlexleSampler)
    m::Float64 = 0.0
    for i in level.indices
        w = sampler.weights[i]
        (w > m) && (m = w)
    end
    return m
end

"""
    printFlexleSampler(sampler, name="")

Print `sampler` with optional label `name`.
"""
function printFlexleSampler(sampler::FlexleSampler; name::String="")
    @printf "\nFlexleSampler %s\n" name
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

# Maintenance methods

"""
    addToFlexLevel!(i, level, sampler, update_sampler_sum=true)

Place the element of index `i` in `sampler` into `level`.

Also updates `level.sum` and (iff `update_sampler_sum`) `sampler.sum` to reflect the addition of the weight of element `i`
to `level` and `sampler`, respectively.

# Note on `update_sampler_sum`

Keyword argument `update_sampler_sum` should only be `false` when `addToFlexLevel!` is used in conjunction with
[`removeFromFlexLevel!`](@ref) to update the weight of an existing element `i` in `sampler`, as in the
following call pattern: 
```
removeFromFlexLevel!(i, old_level, sampler, update_sampler_sum=false)
addToLevelFlexLevel(i, new_level, sampler, update_sampler_sum=false)
sampler.sum += new_i_weight - old_i_weight
``` 
This option is provided for performance purposes, as reading and writing to `sampler.sum` can be expensive. The caller MUST update
`sampler.sum` themselves if calling with `update_sampler_sum=false`.
"""
function addToFlexLevel!(i::Int64, level::FlexLevel, sampler::FlexleSampler; update_sampler_sum::Bool=true)
    push!(level.indices, i)
    w::Float64 = sampler.weights[i]
    level.sum += w
    (update_sampler_sum) && (sampler.sum += w)
    if w > level.max
        level.max = w
    end
    if i <= length(sampler.index_positions)
        sampler.index_positions[i] = length(level.indices)
    else
        push!(sampler.index_positions, length(level.indices))
    end
end

"""
    removeFromFlexLevel!(i, level, sampler, update_sampler_sum=true)

Remove the element of index `i` in `sampler` from `level`.

Also updates `level.sum` and (iff `update_sampler_sum`) `sampler.sum` to reflect the removal of the weight of element `i`
from `level` and `sampler`, respectively.

# Note on `update_sampler_sum`

Keyword argument `update_sampler_sum` should only be `false` when `removeFromFlexLevel!` is used in conjunction with
[`addToFlexLevel!`](@ref) to update the weight of an existing element `i` in `sampler`, as in the following call pattern: 
```
removeFromFlexLevel!(i, old_level, sampler, update_sampler_sum=false)
addToLevelFlexLevel(i, new_level, sampler, update_sampler_sum=false)
sampler.sum += new_i_weight - old_i_weight
``` 
This option is provided for performance purposes, as reading and writing to `sampler.sum` can be expensive. The caller MUST update
`sampler.sum` themselves if calling with `update_sampler_sum=false`.
"""
function removeFromFlexLevel!(i::Int64, level::FlexLevel, sampler::FlexleSampler; update_sampler_sum::Bool=true)
    w::Float64 = sampler.weights[i]
    len = length(level.indices)
    idx = sampler.index_positions[i]
    last = pop!(level.indices)
    if idx != len   # take last index and put it in the place of the one to be removed, unless the last one is itself to be removed
        level.indices[idx] = last
        sampler.index_positions[last] = idx
    end
    sampler.index_positions[i] = 0
    level.sum -= w
    (update_sampler_sum) && (sampler.sum -= w)
    if !levelIsPopulated(level)
        level.max = 0.0
    elseif w == level.max
        level.max = maxLevelWeight(level, sampler)
    end
end

"""
    extendLevels!(bounds, levels)

Extend `levels` to contain all appropriate `FlexLevel`s up to and including that specified by `bounds`.

Throws an error if a level with such bounds already exists in `levels`.
"""
function extendLevels!(bounds::Tuple{Float64,Float64}, levels::Vector{FlexLevel})
    if bounds[1] * 2.0 != bounds[2]
        throw("Invalid bounds - must be two adjacent powers of 2")
    end

    l_bound = bounds[1]
    extend_up = l_bound > levels[begin].bounds[1]
    extend_down = l_bound < levels[end].bounds[1]
    if extend_up
        num_new_levels = logDist(levels[begin].bounds[1], l_bound)
        pre = Vector{FlexLevel}(undef, num_new_levels)
        for i in 1:num_new_levels
            u_bound = l_bound * 2.0
            pre[i] = FlexLevel((l_bound, u_bound), 0.0, 0.0, Vector{Int64}())
            l_bound /= 2.0
        end
        prepend!(levels, pre)
    elseif extend_down
        num_new_levels = logDist(l_bound, levels[end].bounds[1])
        post = Vector{FlexLevel}(undef, num_new_levels)
        for i in num_new_levels:-1:1
            u_bound = l_bound * 2.0
            post[i] = FlexLevel((l_bound, u_bound), 0.0, 0.0, Vector{Int64}())
            l_bound = u_bound
        end
        append!(levels, post)
    else
        throw("levels already contains FlexLevel of specified bounds")
    end
end

"""
    trimTrailingLevels!(sampler)

Remove all empty `FlexLevel`s from the front and back of `sampler`.
"""
function trimTrailingLevels!(sampler::FlexleSampler)
    first = findfirst(levelIsPopulated, sampler.levels)
    last = findlast(levelIsPopulated, sampler.levels)
    sampler.levels = sampler.levels[first:last]
end


# Sampling methods

"""
    cdfSample(sampler)

Randomly select a `FlexLevel` from `sampler` using inverse transform sampling.

Also returns a "free" random number in [0, 1) for use in subsequent rejection sampling
(see [`rejectionSample`](@ref)).
"""
@inline function cdfSample(sampler::FlexleSampler)
    local chosen_level::FlexLevel
    norm_rand_n = rand() * sampler.sum
    cum_sum = 0.0
    for level in sampler.levels
        cum_sum += level.sum
        chosen_level = level
        (cum_sum > norm_rand_n) && break
    end
    return chosen_level, (norm_rand_n - cum_sum + chosen_level.sum) / chosen_level.sum
end

"""
    rejectionSample(rand_n, level, weights)

Randomly select an index from `level` using rejection sampling given a starting `rand_n`
and a `Vector` of `weights`.

`rand_n` is generated in the course of inverse transform sampling (see [`cdfSample`](@ref))
performed prior to rejection sampling.
"""
@inline function rejectionSample(rand_n::Float64, level::FlexLevel, weights::Vector{Float64})
    # while true
    #     r = rand_n * length(level.indices)
    #     d::Float64, i::Int64 = modf(r)
    #     idx = level.indices[i + 1]   # +1 to offset for 1-indexing
    #     if weights[idx]/level.max > d
    #         return idx
    #     end
    #     rand_n = rand()
    # end
    while true
        i = rand(level.indices)
        if weights[i] > rand_n * level.max
            return i
        end
        rand_n = rand()
    end
end
