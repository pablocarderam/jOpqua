using Printf
using Plots
using StatsPlots
using Random
using BenchmarkTools

# Flexible binary level rejection sampling

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

# Module-internal methods

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
    # level.index_positions[FastHashIndex(i)] = length(level.indices)
    # push!(sampler.index_positions, length(level.indices))
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
    # fhi = FastHashIndex(i)
    idx = sampler.index_positions[i] # level.index_positions[fhi] 
    last = pop!(level.indices)
    # delete!(level.index_positions, i) # fhi)
    if idx != len   # take last index and put it in the place of the one to be removed, unless the last one is itself to be removed
        level.indices[idx] = last
        sampler.index_positions[last] = idx # FastHashIndex(last)] = idx
    end
    sampler.index_positions[i] = 0
    level.sum -= w
    (update_sampler_sum) && (sampler.sum -= w)
    if !levelIsPopulated(level)
        level.max = 0.0
    elseif w == level.max
        level.max = maxLevelWeight(level, sampler)
        # level.max = maximum(i -> sampler.weights[i], level.indices) # potential slowdown: maximum updated by searching the whole vector if max is removed
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

"""
    recalculateFlexleStats!(sampler)

Trim trailing empty levels from `sampler` and recalculate its `sum` from `sampler.weights`.
"""
function recalculateFlexleStats!(sampler::FlexleSampler)
    trimTrailingLevels!(sampler)
    # sampler.max = isempty(sampler.levels) ? 0.0 : sampler.levels[begin].max
    sampler.sum = sum(sampler.weights)
end

"""
    reconstructIndexPositions!(sampler)

Recalculate from scratch the `index_positions` dictionary of each `FlexLevel` in `sampler`.
"""
function reconstructIndexPositions!(sampler::FlexleSampler)
    for level in sampler.levels
        level.index_positions = SparseVector{Int64, Int64}(length(sampler.weights)) # Dict{FastHashIndex, Int64}()
        for i in eachindex(level.indices)
            level.index_positions[level.indices[i]] = i
        end
    end
end

"""
    reconstructIndexPositions!(sampler, i)

Recalculate the `index_positions` dictionary of each `FlexLevel` in `sampler` following the removal of
element `i` from `sampler.weights`.

For performance reasons, takes advantage of that only elements of index at least `i` will require
editing in their appropriate `index_positions` dictionary.

Note that this method of `reconstructIndexPositions!` must be called _after_ `sampler.weights` is modified.
"""
function reconstructIndexPositions!(sampler::FlexleSampler, i::Int64)
    n = length(sampler.weights)
    if i > n
        return
    end
    
    l_idx = getLevel(sampler.weights[i], sampler.levels)
    for idx in i:length(sampler.weights)-1
        l_idx_plus1 = getLevel(sampler.weights[idx+1], sampler.levels)
        if !iszero(sampler.weights[idx])
            if l_idx === l_idx_plus1
                l_idx.index_positions[idx] = l_idx.index_positions[idx+1]
            else
                l_idx.index_positions[idx] = pop!(l_idx.index_positions, idx+1)
            end
        end
        l_idx = l_idx_plus1
    end

    w = sampler.weights[n]
    iszero(w) && return
    d = getLevel(sampler.weights[n], sampler.levels).index_positions
    d[n] = pop!(d, n+1)     # removal of weights[i] means (d[idx+1] => pos) should be updated to (d[idx] => pos) for all idx >= i
end

# Sampling methods

"""
    cdfSample(sampler)

Randomly select a `FlexLevel` from `sampler` using inverse transform sampling.

Also returns a "free" random number in [0, 1) for use in subsequent rejection sampling
(see [`rejectionSample`](@ref)).
"""
function cdfSample(sampler::FlexleSampler)
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
function rejectionSample(rand_n::Float64, level::FlexLevel, weights::Vector{Float64})
    # while true
    #     r = rand_n * length(level.indices)
    #     i = floor(r)
    #     idx = level.indices[Int64(i) + 1]   # +1 to offset for 1-indexing
    #     rand_n = r - i
    #     if weights[idx] > (rand_n * level.max)
    #         return idx
    #     end
    # end
    while true
        i = rand(level.indices)
        if weights[i] > rand_n * level.max
            return i
        end
        rand_n = rand()
    end
end

# function rejectionSample(level::FlexLevel, weights::Vector{Float64})
#     while true
#         i = rand(level.indices)
#         if weights[i] > rand() * level.max
#             return i
#         end
#     end
# end

"""

        --- USER FUNCTIONS ---

"""

"""
    FlexleSampler(weights)

Create a `FlexleSampler` from a `Vector` of `weights`.
"""
function FlexleSampler(weights::AbstractVector{Float64})
    if length(weights) == 0
        throw("Cannot create FlexleSampler from AbstractVector of length 0.")
    end

    w_sum = 0.0
    weights_nonzero = filter(x -> !iszero(x), weights)
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

        logmax = log2(w_max)
        if logmax == ceil(logmax)   # if w_max is a power of 2, can't take its ceiling to get its upper bound; need to add 1 instead
            max_u_log = Int64(logmax) + 1
        else
            max_u_log = Int64(ceil(logmax))
        end

        min_l_log_f = floor(log2(w_min))
        min_l_log = Int64(min_l_log_f)

        l_bound = 2.0^min_l_log_f
        num_levels = max_u_log - min_l_log     # e.g. -2,5 ==> 7 levels [4,3,2,1,0,-1,-2]
    end

    levels = Vector{FlexLevel}(undef, num_levels)   # add check for unreasonable number of levels before allocating space?
    index_positions = zeros(Int64, length(weights))

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
                index_positions[i] = length(l.indices)
                w_sum += w
            end
        end
    end

    return FlexleSampler(levels, Vector(weights), w_sum, index_positions)
end

"""
    updateFlexleSamplerWeight!(sampler, i, w)

Update the weight of element `i` in `sampler` to be `w`, returning the difference between the new and old values of `sampler.weights[i]`.
"""
function updateFlexleSamplerWeight!(sampler::FlexleSampler, i::Int64, w::Float64)
    # remove from current level
    w_old::Float64 = sampler.weights[i]
    levels = sampler.levels
    if !iszero(w_old)
        from = getLevel(w_old, levels)
        removeFromFlexLevel!(i, from, sampler, update_sampler_sum=false)
    end

    # update weights vector - has to be done between removal and addition
    sampler.weights[i] = w

    # move to new level
    if !iszero(w)
        bounds = logBounds(w)
        if !inSampler(bounds, sampler)
            extendLevels!(bounds, levels)
        end
        to = getLevel(bounds, levels)
        addToFlexLevel!(i, to, sampler, update_sampler_sum=false)
    end
    sampler.sum += w - w_old

    # trim excess levels (to save time, only if removed element from a level on the end)
    if !iszero(w_old) && (from === levels[begin] || from === levels[end]) && isempty(from.indices)
        trimTrailingLevels!(sampler)
    end

    return w - w_old
end

"""
    addToFlexleSampler!(sampler, w)

Add a new element with weight `w` to `sampler`, updating all fields accordingly.

Returns the new number of weights in `sampler`, or equivalently, the index corresponding to the new element added.
"""
function addToFlexleSampler!(sampler::FlexleSampler, w::Float64)
    push!(sampler.weights, w)
    if !iszero(w)
        bounds = logBounds(w)
        if !inSampler(bounds, sampler)
            extendLevels!(bounds, sampler.levels)
        end
        to = getLevel(bounds, sampler.levels)
        addToFlexLevel!(length(sampler.weights), to, sampler)
    end
    return length(sampler.weights)
end

"""
    removeFromFlexleSampler!(sampler, i)

Remove element `i` from `sampler` completely, updating all fields accordingly.

Returns the new number of weights in `sampler`.

All elements of index `>i` are updated to account for the removal of element `i`.
"""
function removeFromFlexleSampler!(sampler::FlexleSampler, i::Int64)
    w::Float64 = sampler.weights[i]
    if !iszero(w)
        bounds = logBounds(w)
        from = getLevel(bounds, sampler.levels)
        removeFromFlexLevel!(i, from, sampler)
    end
    deleteat!(sampler.weights, i)
    for level in sampler.levels
        for j in eachindex(level.indices)
            if level.indices[j] > i
                level.indices[j] -= 1
            end
        end
    end
    # reconstructIndexPositions!(sampler)
    # reconstructIndexPositions!(sampler, i)
    deleteat!(sampler.index_positions, i)
    if !iszero(w) && (from === sampler.levels[begin] || from === sampler.levels[end]) && isempty(from.indices)
        trimTrailingLevels!(sampler)
    end
    return length(sampler.weights)
end

"""
    flexleSample(sampler)

Take a single random sample from `sampler`, returning the index of the element sampled.

Samples by inverse transform sampling (see `cdfSample`(@ref)) to select a `FlexLevel` in `sampler`, then rejection
sampling (see `rejectionSample`(@ref)) an index from said `FlexLevel`.
"""
function flexleSample(sampler::FlexleSampler)
    level, rand_n = cdfSample(sampler)
    return rejectionSample(rand_n, level, sampler.weights)
end

"""
    getWeight(sampler, i)

Get the weight of element `i` in `sampler`.
"""
function getWeight(sampler::FlexleSampler, i::Int64)
    return sampler.weights[i]
end

# Testing

"""
    verifyFlexleSampler(sampler, name="")

Performs several checks on `sampler` to ensure that it is internally consistent, returning the total number
of inconsistencies detected.

`name` allows the caller to give an optional label to `sampler` for pretty printing purposes.

Checks are as follows:
- index placement in `sampler.levels`
    - every index of non-zero weight in `sample.weights` appears in some level
    - every index appearing in a level is an index that:
        - exists in `sampler.weights`, and if so,
        - belongs in this level according to its weight
- construction of `index_positions` at each level in `sampler.levels`
    - each index appears in the appropriate level's `index_positions` corresponding to its position in `level.indices`
    - no index appears in an `index_positions` corresponding to a position at which it does not exist
- `sampler` stats
    - each level's `sum` is (approximately) equal to the sum of the weights of the indexes it holds
    - each level's `max` is equal to the maximum of the weights of the indexes it holds
    - `sampler.sum` is (approximately) equal to the sum of `sampler.weights`
"""
function verifyFlexleSampler(sampler::FlexleSampler; name::String="")
    @printf "Verifying FlexleSampler %s...\n" name
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

    # all indices in levels recorded correctedly in index_positions?
    d = sampler.index_positions
    for level in sampler.levels
        for pos in eachindex(level.indices)
            idx = level.indices[pos]
            if d[idx] != pos
                errors += 1
                @printf "Error %i: element %i (level (%f, %f), weight %f) not recorded as position %i in index_positions\n" errors idx level.bounds[1] level.bounds[2] sampler.weights[idx] pos
            end
        end
    end
    for idx in eachindex(d)
        pos = d[idx]
        level = getLevel(sampler.weights[idx], sampler.levels)
        if !iszero(pos) && ((pos > length(level.indices)) || !(level.indices[pos] == idx))
            errors += 1
            @printf "Error %i: element %i (level (%f, %f), weight %f) incorrectly recorded as position %i in index_positions\n" errors idx level.bounds[1] level.bounds[2] sampler.weights[idx] pos
        end
    end
    # show(d)

    # all level sums/maxes correct?
    overall_sum = 0.0
    for level in sampler.levels
        empty = isempty(level.indices)
        expected_sum = empty ? 0.0 : sum(sampler.weights[i] for i in level.indices)
        if !approxeq(expected_sum, level.sum)
            errors += 1
            @printf "Error %i: level (%f, %f) sum incorrect (expected %f, got %f)\n" errors level.bounds[1] level.bounds[2] expected_sum level.sum
        end
        overall_sum += expected_sum

        expected_max = empty ? 0.0 : maximum(sampler.weights[i] for i in level.indices)
        if expected_max != level.max
            errors += 1
            @printf "Error %i: level (%f, %f) max incorrect (expected %f, got %f)\n" errors level.bounds[1] level.bounds[2] expected_max level.max
        end
    end
    if !approxeq(overall_sum, sampler.sum)     # correct for probable floating point error
        errors += 1
        @printf "Error %i: overall sampler sum incorrect (expected %f, got %f)" errors overall_sum sampler.sum
    end

    @printf"Final error count: %i\n" errors
    return errors
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

"""
    testFlexleSampler()

Test the basic internal maintenance of the `FlexleSampler` struct.

Assesses:
- initialization
- updating of weights
- extending of levels
"""
function testFlexleSampler(seed=0)
    Random.seed!(seed)

    # --- INITIALIZATION ---
    println("\nTEST: FlexleSampler initialization from weights vectors\n")
    names = Vector{String}()
    weightss = Vector{Vector{Float64}}()
    samplers = Vector{FlexleSampler}()

    push!(names, "(random weights vector, len 100)")
    push!(weightss, rand(100))
    push!(samplers, FlexleSampler(weightss[end]))
    verifyFlexleSampler(samplers[end], name=names[end])

    push!(names, "(all zero weights vector, len 100)")
    push!(weightss, zeros(100))
    push!(samplers, FlexleSampler(weightss[end]))
    verifyFlexleSampler(samplers[end], name=names[end])

    push!(names, "(some zeros in weights vector, len 100)")
    push!(weightss, zeros(100) .+ [((n >= 0.5) || (n < 0.25) ? n : 0.0) for n in rand(100)])
    push!(samplers, FlexleSampler(weightss[end]))
    verifyFlexleSampler(samplers[end], name=names[end])


    # --- MOVE BETWEEN EXISTING LEVELS ---

    push!(names, "(move between levels, test 1)")
    push!(weightss, [1.5 * n for n in 1.0:100.0])
    push!(samplers, FlexleSampler(weightss[end]))
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(move test, pre-move)")

    # between two levels - element in middle of indices vector
    idx = 24
    w_new = 24.0 # now belongs in (16.0, 32.0) instead of (32.0, 64.0)
    updateFlexleSamplerWeight!(samplers[end], idx, w_new)
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(move test, changed weights[24] from 36.0 to 24.0)")

    # between two levels - element at end of indices vector
    idx = 85
    w_new = 60.0 # now belongs in (32.0, 64.0) instead of (64.0, 128.0)
    updateFlexleSamplerWeight!(samplers[end], idx, w_new)
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(move test, changed weights[85] from 127.5 to 60.0)")

    # # within same level
    # idx = 96
    # w_new = 140.0 # still belongs in (128.0, 256.0)
    # updateFlexleSamplerWeight!(samplers[end], idx, w_new)
    # # printFlexleSampler(samplers[end])
    # verifyFlexleSampler(samplers[end], name="(move test, changed weights[96] from 144.0 to 140.0)")

    # # between two levels - "from" level is now empty
    # idx = 1
    # w_new = 3.2 # now belongs in (2.0, 4.0) instead of (1.0, 2.0)
    # updateFlexleSamplerWeight!(samplers[end], idx, w_new)
    # # printFlexleSampler(samplers[end])
    # verifyFlexleSampler(samplers[end], name="(move test, changed weights[1] from 1.5 to 3.2)")


    # # --- EXTEND LEVELS ---

    # # upwards by 1
    # extendLevels!((256.0, 512.0), samplers[end].levels)
    # # printFlexleSampler(samplers[end])
    # verifyFlexleSampler(samplers[end], name="(extend test, single level upwards)")

    # # upwards by several
    # extendLevels!((4096.0, 8192.0), samplers[end].levels)
    # # printFlexleSampler(samplers[end])
    # verifyFlexleSampler(samplers[end], name="(extend test, several levels upwards)")

    # # downwards by one
    # extendLevels!((0.5, 1.0), samplers[end].levels)
    # # printFlexleSampler(samplers[end])
    # verifyFlexleSampler(samplers[end], name="(extend test, single level downwards)")

    # # downwards by several
    # extendLevels!((0.015625, 0.03125), samplers[end].levels)
    # # printFlexleSampler(samplers[end])
    # verifyFlexleSampler(samplers[end], name="(extend test, several levels downards)")

    # # ERROR: new bound matches existing upper bound
    # try
    #     extendLevels!((4096.0, 8192.0), samplers[end].levels)
    # catch e
    #     println("correctly identified error: ", e)
    # end
    # # printFlexleSampler(samplers[end])
    # verifyFlexleSampler(samplers[end], name="(extend test, attempt to add existing upper bound)")

    # #  ERROR: new bound matches existing lower bound
    # try
    #     extendLevels!((0.015625, 0.03125), samplers[end].levels)
    # catch e
    #     println("correctly identified error: ", e)
    # end
    # # printFlexleSampler(samplers[end])
    # verifyFlexleSampler(samplers[end], name="(extend test, attempt to add existing lower bound)")

    # # ERROR: new bound matches occupied intermediate level
    # try
    #     extendLevels!((2.0, 4.0), samplers[end].levels)
    # catch e
    #     println("correctly identified error: ", e)
    # end
    # # printFlexleSampler(samplers[end])
    # verifyFlexleSampler(samplers[end], name="(extend test, attempt to add existing occupied intermediate level)")

    # # ERROR: new bound matches empty intermediate level
    # try
    #     extendLevels!((1.0, 2.0), samplers[end].levels)
    # catch e
    #     println("correctly identified error: ", e)
    # end
    # # printFlexleSampler(samplers[end])
    # verifyFlexleSampler(samplers[end], name="(extend test, attempt to add existing empty intermediate level)")
end

"""
    levelSamplingExpectedValue(sampler, n=1)

Produce a `Dict` mapping (bounds of) levels in `sampler` to their expected number of selections given `n` samples.

For default value `n=1`, each level is simply assigned the probability of its selection.
"""
function levelSamplingExpectedValue(sampler::FlexleSampler; n::Int64=1)
    norm = Float64(n) / sampler.sum
    d = Dict{Tuple{Float64,Float64},Float64}()
    for level in sampler.levels
        d[level.bounds] = norm * level.sum
    end
    return d
end

"""
    testLevelSampling(sampler, n)

Randomly select a level from `sampler` `n` times and return a `Dict` mapping levels to empirical number of draws.
"""
function testLevelSampling(sampler::FlexleSampler, n::Int64)
    d = Dict{Tuple{Float64,Float64},Int64}()
    for level in sampler.levels
        d[level.bounds] = 0
    end

    for _ in 1:n
        bounds = cdfSample(sampler)[1].bounds
        d[bounds] += 1
    end

    return d
end

"""
    indexSamplingExpectedValue(level, weights, n=1)

Produce a `Dict` mapping indexes in `level` to their expected number of selections (according to `weights`) given `n` samples.

For default value `n=1`, each index is simply assigned the probability of its selection. 
"""
function indexSamplingExpectedValue(level::FlexLevel, weights::Vector{Float64}; n::Int64=1)
    norm = Float64(n) / level.sum
    d = Dict{Int64,Float64}()
    for i in level.indices
        d[i] = norm * weights[i]
    end
    return d
end

"""
    testIndexSampling(level, weights, n)

Randomly select an index from `level` `n` times according to `weights` and return a `Dict` mapping indexes to empirical number of draws.
"""
function testIndexSampling(level::FlexLevel, weights::Vector{Float64}, n::Int64)
    d = Dict{Int64,Int64}()
    for i in level.indices
        d[i] = 0
    end

    for _ in 1:n
        i = rejectionSample(rand(), level, weights)
        d[i] += 1
    end

    return d
end

"""
    flexleSamplingExpectedValue(sampler, n=1)

Produce a `Dict` mapping indexes in `sampler` to their expected number of selections given `n` samples.

For default value `n=1`, each index is simply assigned the probability of its selection. 
"""
function flexleSamplingExpectedValue(sampler::FlexleSampler; n::Int64=1)
    norm = Float64(n) / sampler.sum # (l.sum / sampler.sum) * (Float64(n) / l.sum) -- l terms cancel out
    d = Dict{Int64,Float64}()
    for i in eachindex(sampler.weights)
        d[i] = norm * sampler.weights[i]
    end

    # for l in sampler.levels
    #     for i in l.indices
    #         d[i] = norm * sampler.weights[i]
    #     end
    # end
    return d
end

"""
    testFlexleSampling(sampler, n)

Randomly select an index from `sampler` `n` times and return a `Dict` mapping indexes to empirical number of draws.
"""
function testFlexleSampling(sampler::FlexleSampler, n::Int64)
    d = Dict{Int64,Int64}()
    for i in eachindex(sampler.weights)
        d[i] = 0
    end

    for _ in 1:n
        i = flexleSample(sampler)
        d[i] += 1
    end

    return d
end

"""
    testUserFunctions()

Test the internal maintenance of the `FlexleSampler` struct in response to user commands.

Assesses:
- initialization
- updating of weights
- addition of new weights
- removal of existing weights
"""
function testUserFunctions()
    # test create
    w0 = zeros(1000)
    s0 = FlexleSampler(w0)
    verifyFlexleSampler(s0)

    w1 = w0 .+ [((i < 0.25 || 0.50 < i) ? i : 0.0) for i in rand(1000)]
    s1 = FlexleSampler(w1)
    verifyFlexleSampler(s1)

    w2 = zeros(100) .+ [(i > 0.5 ? i : 0.0) for i in 0.0:99.0]
    s2 = FlexleSampler(w2)
    verifyFlexleSampler(s2)

    # test update
    for i in 2:10:92
        updateFlexleSamplerWeight!(s2, i, 100*rand())
    end
    verifyFlexleSampler(s2)
    
    updateFlexleSamplerWeight!(s2, 87, 0.0)     # make non-zero element 0.0
    verifyFlexleSampler(s2)

    updateFlexleSamplerWeight!(s2, 1, 3000.0)   # make 0.0 element non-zero
    verifyFlexleSampler(s2)

    # test add
    addToFlexleSampler!(s2, 99.5)       # in existing level (occupied)
    addToFlexleSampler!(s2, 480.0)      # in existing level (empty)
    addToFlexleSampler!(s2, 10000.0)    # in non-existing level (above)
    addToFlexleSampler!(s2, 0.001)      # in non-existing level (below)
    addToFlexleSampler!(s2, 0.0)        # zero
    verifyFlexleSampler(s2)

    # test remove
    removeFromFlexleSampler!(s2, 70)      # in middle level, does not empty level
    removeFromFlexleSampler!(s2, 101)     # in middle level, empties level
    removeFromFlexleSampler!(s2, 101)     # empty top level
    removeFromFlexleSampler!(s2, 101)     # empty bottom level
    removeFromFlexleSampler!(s2, 86)      # zero
    removeFromFlexleSampler!(s2, 1)       # first element
    verifyFlexleSampler(s2)
end



"""

        --- EXAMPLES AND GRAPHS ---

"""

function plotCDF(weights::Vector{Float64}, filepath::String, d::Int64; name::String="cdf", xlabel="Element", n::Union{Float64, Nothing}=nothing)
    s = 0.0
    m = length(weights)
    indices = zeros(Int64, m+1)
    cum_sums = zeros(Float64, m+1)
    for i in eachindex(weights)
        indices[i+1] = i
        s += weights[i]
        cum_sums[i+1] = s
    end

    
    p = plot(indices, cum_sums, title="CDF sampling", xlabel=xlabel, ylabel="Cumulative sum", linetype=:steppost, legend=false, size=(d, d), linewidth=2)
    if !isnothing(n)
        k = findfirst(x -> x >= n, cum_sums) - 1
        plot!(p, indices, [n for _ in indices], linewidth=3)
        plot!(p, [k, k], [0, n], linewidth=3)
    end
    png(p, filepath * name)

    return p
end

function plotRejection(weights::Vector{Float64}, filepath::String, d::Int64)
    m = maximum(weights)
    W = hcat(m .- weights, weights)
    p = groupedbar(W, bar_width=1, bar_position=:stack, title="Rejection sampling", xlabel="Element", ylabel="Weight", legend=false, size=(d, d))
    png(p, filepath * "rej")
    return p
end

function makeFlexleExampleGraphs(; seed::Int64=3, filepath="examples/flexlegraphs/")
    d = 800
    d34 = 3*div(d, 4)
    Random.seed!(seed)
    w = 5 * rand(15)
    n = rand() * sum(w)
    plotCDF(w, filepath, d34, n=n)
    plotRejection(w, filepath, d34)
    s = FlexleSampler(w)
    printFlexleSampler(s)
    f = []
    f_lsize = maximum(length(l.indices) for l in s.levels)
    for i in eachindex(s.levels)
        l = s.levels[i] 
        p = [s.weights[j] for j in l.indices]
        push!(f, hcat(maximum(p) .- p, p))
        bound = string(l.bounds)
        groupedbar(f[end], bar_width=1, bar_position=:stack, legend=false, xticks=(1:1:f_lsize, l.indices), xlims=(0, f_lsize+1), ylabel=bound, size=(d/2, d/4))
        png(filepath * "f" * string(i))
    end

    sums = [l.sum for l in s.levels]
    sums_r = reverse(sums)
    plotCDF(sums, filepath, d34, name="f_cdf", xlabel="Level")
    plotCDF(sums_r, filepath, d34, name="f_cdf_r", xlabel="Level")
end

function testSampleRuntime!(w::Vector{Vector{Float64}}, w_sums::Vector{Float64}, sampler::FlexleSampler, r1::Int64, r2::Int64)
    r1 = 1:r1
    r2 = 1:r2
    
    println("CDF sampling:")    
    #randChoose(rand(), weights, w_tot)
    # @benchmark timeRandChoose($r1, $r2, $w, $w_tot)
    display(@benchmark timeRandChoose($r1, $r2, $w, $w_sums))
    # display(@benchmark randChoose(0.5, $weights, $w_tot, regenerate_rand=false))

    println("Flexle sampling:")
    # # flexleSample(sampler)
    # @benchmark timeFlexleSample($r1, $r2, $sampler)
    display(@benchmark timeFlexleSample($r1, $r2, $sampler))
    # display(@benchmark(flexleSample($sampler)))
end

function timeRandChoose(r1::UnitRange{Int64}, r2::UnitRange{Int64}, weights::Vector{Vector{Float64}}, w_sums::Vector{Float64})
    for i in r1
        rand_n = rand()
        for j in r2
            k = i*(j-1) + j
            y, rand_n = randChoose(rand_n, weights[k], w_sums[k], regenerate_rand=true)
        end
    end
end

function timeFlexleSample(r1::UnitRange{Int64}, r2::UnitRange{Int64}, sampler::FlexleSampler)
    for i in r1
        for j in r2
            y = flexleSample(sampler)
        end
    end
end

"""
    testRuntime01(h=10000)

Compare the runtime of performing 1000 samples using CDF (i.e. `randChoose`, inverse transform) versus Flexle sampling on a
`Vector` of `h` weights, all either `0.0` or `1.0`.
"""
function testRuntime01(; h::Int64=10000, seed=0)
    Random.seed!(seed)
    w01::Vector{Float64} = zeros(h)
    for i in eachindex(w01)
        if rand() > 0.5
            w01[i] = 1
        end
    end

    r1, r2 = 200, 5
    w = Vector{Vector{Float64}}()
    w_sums = Vector{Float64}()
    for j in 1:(r1*r2)
        push!(w, zeros(h))
        for i in 1:h
            if rand() > 0.5
                w[j][i] = 1
            end
        end
        push!(w_sums, sum(w[j]))
    end


    testSampleRuntime!(w, w_sums, FlexleSampler(w01), r1, r2)
end

"""
    testRuntimeUniform(h=10000)

Compare the runtime of performing 1000 samples using CDF (i.e. `randChoose`, inverse transform) versus Flexle sampling on a
`Vector` of `h` weights, all randomly selected from a uniform distribution between `0.0` and `1.0`.
"""
function testRuntimeUniform(; h::Int64=10000, seed=0)
    Random.seed!(seed)
    wu::Vector{Float64} = rand(h)

    r1, r2 = 200, 5
    w = Vector{Vector{Float64}}()
    w_sums = Vector{Float64}()
    for j in 1:(r1*r2)
        push!(w, rand(h))
        push!(w_sums, sum(w[j]))
    end

    testSampleRuntime!(w, w_sums, FlexleSampler(wu), r1, r2)
end

function testStandardUpdateWeight!(weights::Vector{Float64}, updates::Tuple{Vector{Int64}, Vector{Float64}}, sum::Float64)
    i = updates[1]
    v = updates[2]
    for idx in eachindex(i)
        sum -= weights[i[idx]]
        weights[i[idx]] = v[idx]
        sum += v[idx]
    end
    return sum
end

function testFlexleUpdateWeight!(sampler::FlexleSampler, updates::Tuple{Vector{Int64}, Vector{Float64}})
    i = updates[1]
    v = updates[2]
    for idx in eachindex(i)
        updateFlexleSamplerWeight!(sampler, i[idx], v[idx])
    end
    # verifyFlexleSampler(sampler)
end

"""
    testUpdateWeight(h=10000, n=1000, seed=0)

Compare the runtime of performing `n` random updates to a collection of `h` weights when storing said weights as
a simple `Vector` versus a `FlexleSampler`.
"""
function testUpdateWeight(; h::Int64=10000, n::Int64=1000, seed=0)
    Random.seed!(seed)
    w = rand(h)
    s = FlexleSampler(w)
    c = sum(w)
    updates = (rand(1:h, n), rand(n))

    println("Vector update weight:")    
    display(@benchmark testStandardUpdateWeight!($w, $updates, $c))

    println("Flexle update weight:")
    display(@benchmark testFlexleUpdateWeight!($s, $updates))
end

function testStandardAddWeight!(weights::Vector{Float64}, new_weights::Vector{Float64}, sum::Float64)
    for i in eachindex(new_weights)
        push!(weights, new_weights[i])
        sum += new_weights[i]
    end
    return sum
end

function testFlexleAddWeight!(sampler::FlexleSampler, new_weights::Vector{Float64})
    for i in eachindex(new_weights)
        addToFlexleSampler!(sampler, new_weights[i])
    end
    # verifyFlexleSampler(sampler)
end

"""
    testAddWeight(h=10000, n=1000, seed=0)

Compare the runtime of adding `n` random values to a starting collection of `h` weights when storing said weights as
a simple `Vector` versus a `FlexleSampler`.
"""
function testAddWeight(; h::Int64=10000, n::Int64=1000, seed=0)
    Random.seed!(seed)
    w = rand(h)
    s = FlexleSampler(w)
    c = sum(w)
    new_weights = rand(n)

    println("Vector add weight:")    
    display(@benchmark testStandardAddWeight!($w, $new_weights, $c))

    println("Flexle add weight:")
    display(@benchmark testFlexleAddWeight!($s, $new_weights))
end

function testStandardRemoveWeight!(weights::Vector{Float64}, indices::Vector{Int64}, sum::Float64)
    for i in eachindex(indices)
        sum -= weights[indices[i]]
        deleteat!(weights, indices[i])
    end
    return sum
end

function testFlexleRemoveWeight!(sampler::FlexleSampler, indices::Vector{Int64})
    for i in eachindex(indices)
        removeFromFlexleSampler!(sampler, indices[i])
    end
    # verifyFlexleSampler(sampler)
end

"""
    testRemoveWeight(h=10000, n=1000, seed=0)

Compare the runtime of removing `n` random elements from a collection of `h` weights when storing said weights as
a simple `Vector` versus a `FlexleSampler`.
"""
function testRemoveWeight(; h::Int64=10000, n::Int64=1000, seed=0)
    Random.seed!(seed)
    w1 = rand(h)
    w2 = copy(w1)
    s = FlexleSampler(w2)
    c = sum(w1)
    indices = Vector{Int64}()
    for i in 1:n
        push!(indices, rand(1:h-i+1))
    end

    println("Vector remove weight:")    
    @time testStandardRemoveWeight!(w1, indices, c)

    println("Flexle remove weight:")
    @time testFlexleRemoveWeight!(s, indices)
end

"""
    flexleTestSuite()

Calls a series of Flexle test function to ensure working status before deployment.
"""
function flexleTestSuite()
    # functionality
    testFlexleSampler()
    testUserFunctions()

    # runtime
    testRuntime01()
    testRuntimeUniform()
    testUpdateWeight()
    testAddWeight()
    testRemoveWeight()
end
