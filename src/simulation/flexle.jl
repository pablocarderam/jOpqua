using Printf

# Flexible binary level rejection sampling

# Structs

mutable struct FlexLevel
    bounds::Tuple{Float64,Float64}
    sum::Float64
    max::Float64
    indices::Vector{Int64}
end

mutable struct FlexleSampler
    levels::Vector{FlexLevel}
    weights::AbstractVector{Float64}
    sum::Float64
end

# Initialization

function newFlexleSampler(weights::AbstractVector{Float64})
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

    return FlexleSampler(levels, weights, w_sum)
end

function newFlexLevel(i::Int64, w::Float64)
    return FlexLevel(logBounds(w), w, w, [i])
end

# Maintenance methods

"""
    levelIndex(w, u)

Given a weight `w`, returns the index of the level in some `FlexleSampler.levels` with maximum upper bound `2^u` where `w` belongs.

Returns `0` if `w` is `0.0`, indicating that `w` belongs in no level.

# Examples

`levelIndex(14.2, 6)` ==>  `3`

The value `14.2` belongs in the `(8.0, 16.0)` level, which in a `FlexleSampler` that
starts with a level of bounds `(32.0, 64.0)` (`64` being `2^6`) is at `levels[3]`.

`levelIndex(8.0, 4)` ==> `1`

`levelIndex(8.0, 5)` ==> `2`
"""
function levelIndex(w::Float64, u::Int64)
    return iszero(w) ? 0 : u - Int64(floor(log2(w)))
end

function levelIndex(bounds::Tuple{Float64,Float64}, levels::Vector{FlexLevel})
    for i in eachindex(levels)
        (levels[i].bounds == bounds) && (return i)
    end
    return 0
end

function getLevel(bounds::Tuple{Float64,Float64}, levels::Vector{FlexLevel})
    return levels[levelIndex(bounds, levels)]
end

function getLevel(w::Float64, levels::Vector{FlexLevel})
    return getLevel(logBounds(w), levels)
end

function logDist(a::Float64, b::Float64)
    return Int64(floor(log2(b))) - Int64(floor(log2(a)))
end

function inSampler(i::Int64, sampler::FlexleSampler)
    for level in sampler.levels
        if i in level.indices
            return true
        end
    end
    return false
end

function inSampler(bounds::Tuple{Float64,Float64}, sampler::FlexleSampler)
    return (sampler.levels[begin].bounds[1] >= bounds[1]) && (bounds[1] >= sampler.levels[end].bounds[1])     # bounds between largest and smallest levels' bounds (inclusive)
end

function addToFlexLevel!(i::Int64, w::Float64, level::FlexLevel, sampler::FlexleSampler)
    push!(level.indices, i)
    level.sum += w
    sampler.sum += w
    if w > level.max
        level.max = w
    end
end

function removeFromFlexLevel!(i::Int64, w::Float64, level::FlexLevel, sampler::FlexleSampler)
    len = length(level.indices)
    idx = findfirst(x -> x == i, level.indices)
    last = pop!(level.indices)
    if idx != len   # take last index and put it in the place of the one to be removed, unless the last one is itself to be removed
        level.indices[idx] = last
    end
    level.sum -= w
    sampler.sum -= w
    if isempty(level.indices)
        level.max = 0.0
    elseif w == level.max
        level.max = maximum([sampler.weights[j] for j in level.indices]) # potential slowdown: maximum updated by searching the whole vector if max is removed
    end
end

function moveBetweenFlexLevels!(i::Int64, w_old::Float64, to::FlexLevel, from::FlexLevel, sampler::FlexleSampler)
    removeFromFlexLevel!(i, w_old, from, sampler)
    addToFlexLevel!(i, sampler.weights[i], to, sampler)
end

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

function levelIsPopulated(level::FlexLevel)
    return !isempty(level.indices)
end

function trimTrailingLevels!(sampler::FlexleSampler)
    first = findfirst(levelIsPopulated, sampler.levels)
    last = findlast(levelIsPopulated, sampler.levels)
    if !(first == 1 && last == length(sampler.levels))
        sampler.levels = sampler.levels[first:last]
    end
end

function updateFlexleSamplerWeight!(i::Int64, w_old::Float64, sampler::FlexleSampler)
    levels = sampler.levels
    bounds = logBounds(sampler.weights[i])
    if !inSampler(bounds, sampler)
        extendLevels!(bounds, levels)
    end
    to = getLevel(bounds, levels)
    from = getLevel(w_old, levels)
    moveBetweenFlexLevels!(i, w_old, to, from, sampler)

    # trim excess levels (to save time, only if removed element from a level on the end)
    if (from === levels[begin] && isempty(from.indices)) || (from == levels[end] && isempty(to.indices))
        trimTrailingLevels!(sampler)
    end
end

function recalculateSamplerStats!(sampler::FlexleSampler)
    trimTrailingLevels!(sampler)
    # sampler.max = isempty(sampler.levels) ? 0.0 : sampler.levels[begin].max
    sampler.sum = sum(sampler.weights)
end

# Sampling methods

function cdfSample(sampler::FlexleSampler)
    local chosen_level::FlexLevel
    norm_rand_n = rand() * sampler.sum
    cum_sum = 0.0
    for level in sampler.levels
        cum_sum += level.sum
        chosen_level = level
        (cum_sum > norm_rand_n) && break
    end
    return chosen_level
end

function rejectionSample(level::FlexLevel, weights::Vector{Float64})
    while true
        r = rand() * level.indices
        i = floor(r)
        if weights[Int64(i)] > (r - i) * level.max
            # i = rand(level.indices)
            # if weights[i] > rand() * level.max
            return i
        end
    end
end

function flexleSample(sampler::FlexleSampler)
    return rejectionSample(cdfSample(sampler), sampler.weights)
end

# Testing

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

    # all level sums/maxes correct?
    overall_sum = 0.0
    for level in sampler.levels
        empty = isempty(level.indices)
        expected_sum = empty ? 0.0 : sum(sampler.weights[i] for i in level.indices)
        if expected_sum != level.sum
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
    if abs(overall_sum - sampler.sum) > 1e-6     # correct for probably floating point error
        errors += 1
        @printf "Error %i: overall sampler sum incorrect (expected %f, got %f)" errors overall_sum sampler.sum
    end

    @printf"Final error count: %i\n" errors
    return errors
end

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

function testFlexleSampler()
    # --- INITIALIZATION ---
    println("\nTEST: FlexleSampler initialization from weights vectors\n")
    names = Vector{String}()
    weightss = Vector{Vector{Float64}}()
    samplers = Vector{FlexleSampler}()

    push!(names, "(random weights vector, len 100)")
    push!(weightss, rand(100))
    push!(samplers, newFlexleSampler(weightss[end]))
    verifyFlexleSampler(samplers[end], name=names[end])

    push!(names, "(all zero weights vector, len 100)")
    push!(weightss, zeros(100))
    push!(samplers, newFlexleSampler(weightss[end]))
    verifyFlexleSampler(samplers[end], name=names[end])

    push!(names, "(some zeros in weights vector, len 100)")
    push!(weightss, zeros(100) .+ [((n >= 0.5) || (n < 0.25) ? n : 0.0) for n in rand(100)])
    push!(samplers, newFlexleSampler(weightss[end]))
    verifyFlexleSampler(samplers[end], name=names[end])


    # --- MOVE BETWEEN EXISTING LEVELS ---

    push!(names, "(move between levels, test 1)")
    push!(weightss, [1.5 * n for n in 1.0:100.0])
    push!(samplers, newFlexleSampler(weightss[end]))
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(move test, pre-move)")

    # between two levels - element in middle of indices vector
    idx = 24
    w_old = weightss[end][idx]
    w_new = 24.0
    weightss[end][idx] = w_new    # now belongs in (16.0, 32.0) instead of (32.0, 64.0)
    moveBetweenFlexLevels!(idx, w_old,
        samplers[end].levels[levelIndex(logBounds(w_new), samplers[end].levels)],
        samplers[end].levels[levelIndex(logBounds(w_old), samplers[end].levels)],
        samplers[end])
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(move test, changed weights[24] from 36.0 to 24.0)")

    # between two levels - element at end of indices vector
    idx = 85
    w_old = weightss[end][idx]
    w_new = 60.0
    weightss[end][idx] = w_new    # now belongs in (32.0, 64.0) instead of (64.0, 128.0)
    moveBetweenFlexLevels!(idx, w_old,
        samplers[end].levels[levelIndex(logBounds(w_new), samplers[end].levels)],
        samplers[end].levels[levelIndex(logBounds(w_old), samplers[end].levels)],
        samplers[end])
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(move test, changed weights[85] from 127.5 to 60.0)")

    # within same level
    idx = 96
    w_old = weightss[end][idx]
    w_new = 140.0
    weightss[end][idx] = w_new    # still belongs in (128.0, 256.0)
    moveBetweenFlexLevels!(idx, w_old,
        samplers[end].levels[levelIndex(logBounds(w_new), samplers[end].levels)],
        samplers[end].levels[levelIndex(logBounds(w_old), samplers[end].levels)],
        samplers[end])
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(move test, changed weights[96] from 144.0 to 140.0)")

    # between two levels - "from" level is now empty        TODO: trim empty levels on ends of sampler
    idx = 1
    w_old = weightss[end][idx]
    w_new = 3.2
    weightss[end][idx] = w_new    # now belongs in (2.0, 4.0) instead of (1.0, 2.0)
    moveBetweenFlexLevels!(idx, w_old,
        samplers[end].levels[levelIndex(logBounds(w_new), samplers[end].levels)],
        samplers[end].levels[levelIndex(logBounds(w_old), samplers[end].levels)],
        samplers[end])
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(move test, changed weights[1] from 1.5 to 3.2)")


    # --- EXTEND LEVELS ---

    # upwards by 1
    extendLevels!((256.0, 512.0), samplers[end].levels)
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(extend test, single level upwards)")

    # upwards by several
    extendLevels!((4096.0, 8192.0), samplers[end].levels)
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(extend test, several levels upwards)")

    # downwards by one
    extendLevels!((0.5, 1.0), samplers[end].levels)
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(extend test, single level downwards)")

    # downwards by several
    extendLevels!((0.015625, 0.03125), samplers[end].levels)
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(extend test, several levels downards)")

    # ERROR: new bound matches existing upper bound
    try
        extendLevels!((4096.0, 8192.0), samplers[end].levels)
    catch e
        println("correctly identified error: ", e)
    end
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(extend test, attempt to add existing upper bound)")

    #  ERROR: new bound matches existing lower bound
    try
        extendLevels!((0.015625, 0.03125), samplers[end].levels)
    catch e
        println("correctly identified error: ", e)
    end
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(extend test, attempt to add existing lower bound)")

    # ERROR: new bound matches occupied intermediate level
    try
        extendLevels!((2.0, 4.0), samplers[end].levels)
    catch e
        println("correctly identified error: ", e)
    end
    # printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(extend test, attempt to add existing occupied intermediate level)")

    # ERROR: new bound matches empty intermediate level
    try
        extendLevels!((1.0, 2.0), samplers[end].levels)
    catch e
        println("correctly identified error: ", e)
    end
    printFlexleSampler(samplers[end])
    verifyFlexleSampler(samplers[end], name="(extend test, attempt to add existing empty intermediate level)")
end

function levelSamplingExpectedValue(sampler::FlexleSampler; n::Int64=1)
    norm = Float64(n) / sampler.sum
    d = Dict{Tuple{Float64,Float64},Float64}()
    for level in sampler.levels
        d[level.bounds] = norm * level.sum
    end
    return d
end

function testLevelSampling(sampler::FlexleSampler, n::Int64)
    d = Dict{Tuple{Float64,Float64},Int64}()
    for level in sampler.levels
        d[level.bounds] = 0
    end

    for _ in 1:n
        bounds = cdfSample(sampler).bounds
        d[bounds] += 1
    end

    return d
end

function indexSamplingExpectedValue(level::FlexLevel, weights::Vector{Float64}; n::Int64=1)
    norm = Float64(n) / level.sum
    d = Dict{Int64,Float64}()
    for i in level.indices
        d[i] = norm * weights[i]
    end
    return d
end

function testIndexSampling(level::FlexLevel, weights::Vector{Float64}, n::Int64)
    d = Dict{Int64,Int64}()
    for i in level.indices
        d[i] = 0
    end

    for _ in 1:n
        i = rejectionSample(level, weights)
        d[i] += 1
    end

    return d
end

function flexleSamplingExpectedValue(sampler::FlexleSampler; n::Int64=1)
    norm = Float64(n) / sampler.sum # (l.sum / sampler.sum) * (Float64(n) / l.sum) -- l terms cancel out
    d = Dict{Int64,Float64}()
    for i in eachindex(sampler.weights)
        d[i] = 0.0
    end

    for l in sampler.levels
        for i in l.indices
            d[i] = norm * sampler.weights[i]
        end
    end
    return d
end

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
