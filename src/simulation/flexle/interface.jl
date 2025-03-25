# User interface containing all exported Flexle functions

using Revise
using Random

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
    getindex(sampler, i)

Get the weight of element `i` in `sampler`.
"""
function Base.getindex(sampler::FlexleSampler, i::Int64)
    return sampler.weights[i]
end

"""
    setindex!(sampler, w, i)

Set the weight of element `i` in `sampler` equal to `w`, returning the difference between the new and old values of `i`.
"""
function Base.setindex!(sampler::FlexleSampler, w::Float64, i::Int64)
    from::Union{Nothing, FlexLevel} = nothing
    to::Union{Nothing, FlexLevel} = nothing
    levels = sampler.levels
    w_old::Float64 = sampler.weights[i]
    delta::Float64 = w - w_old
    nonzero = !iszero(w_old), !iszero(w)
    if nonzero[1]
        from = getLevel(w_old, levels)
    end
    if nonzero[2]
        bounds = logBounds(w)
        if !inSampler(bounds, sampler)
            extendLevels!(bounds, levels)
        end
        to = getLevel(bounds, levels)
        if from == to
            sampler.weights[i] = w
            to.sum += delta
            sampler.sum += delta
            if w > to.max
                to.max = w
            end
            return delta
        end
    end

    if nonzero[1]
        removeFromFlexLevel!(i, from, sampler, update_sampler_sum=false)
    end
    sampler.weights[i] = w  # weight vector update must be between removal and addition
    if nonzero[2]
        addToFlexLevel!(i, to, sampler, update_sampler_sum=false)
    end
    sampler.sum += delta

    # trim excess levels (to save time, only if removed element from a level on the end)
    if nonzero[1] && (from === levels[begin] || from === levels[end]) && !levelIsPopulated(from)
        trimTrailingLevels!(sampler)
    end

    return delta
end

"""
    push!(sampler, w)

Add a new element with weight `w` to `sampler`, updating all fields accordingly.

Returns the new number of weights in `sampler`, or equivalently, the index corresponding to the new element added.
"""
function Base.push!(sampler::FlexleSampler, w::Float64)
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
    deleteat!(sampler, i)

Remove element `i` from `sampler` completely, updating all fields accordingly.

Returns the new number of weights in `sampler`.

All elements of index `>i` are updated to account for the removal of element `i`.
"""
function Base.deleteat!(sampler::FlexleSampler, i::Int64)
    w::Float64 = sampler.weights[i]
    if !iszero(w)
        bounds = logBounds(w)
        from = getLevel(bounds, sampler.levels)
        removeFromFlexLevel!(i, from, sampler)
    end
    deleteat!(sampler.weights, i)
    for level in sampler.levels
        indices = level.indices
        for j in eachindex(indices)
            if indices[j] > i
                indices[j] -= 1
            end
        end
    end
    deleteat!(sampler.index_positions, i)
    if !iszero(w) && (from === sampler.levels[begin] || from === sampler.levels[end]) && isempty(from.indices)
        trimTrailingLevels!(sampler)
    end
    return length(sampler.weights)
end

"""
    sample(sampler)

Take a single random sample from `sampler`, returning the index of the element sampled.

Samples by inverse transform sampling (see `cdfSample`(@ref)) to select a `FlexLevel` in `sampler`, then rejection
sampling (see `rejectionSample`(@ref)) an index from said `FlexLevel`.
"""
@inline function sample(sampler::FlexleSampler)
    level, rand_n = cdfSample(sampler)
    return rejectionSample(rand_n, level, sampler.weights)
end

