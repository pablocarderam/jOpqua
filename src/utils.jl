using PoissonRandom
using Flexle

catCol(a::AbstractMatrix{Float64}, b::AbstractVector{Float64}) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
# Simeon Schaub https://discourse.julialang.org/t/adding-rows-to-a-matrix-dynamically/52984

zeroTruncatedPoisson(rate::Float64) = 1 + pois_rand(rate - 1.0)
# this is a hack of what possibly should be done
# (but also maybe this is correct and I just shouldn't call it a
# zero-truncated Poisson at all, I have to think #TODO:)
# https://en.wikipedia.org/wiki/Zero-truncated_Poisson_distribution

function binomial(n::Integer, p::Real)
    log_q = log(1.0 - p)
    x = 0
    sum = 0.0
    while true
        sum += log(rand()) / (n - x)
        sum < log_q && break
        x += 1
    end
    return x
end
# This is like 25% faster than using Distributions.jl for n=10 and p=0.1 and almost twice as fast for p=0.001, I benchmarked
# Following user anon94023334 (2017) on https://discourse.julialang.org/t/binomial-distribution-without-distributions-jl/3176/8
# who took it from https://stackoverflow.com/questions/23561551/a-efficient-binomial-random-number-generator-code-in-java
# where it was cited as a variant of Luc Devroye's "Second Waiting Time Method" on page 522 of his text "Non-Uniform Random
# Variate Generation."

function approxeq(a::Float64, b::Float64; t::Float64=1e-9)
    return abs(a - b) < t
end

function approxZero(a::Float64; t::Float64=1e-9)
    return a < t
end

"""
    FlexleSamplers(weights, number)

Create a vector with separate `FlexleSampler` built from a `Vector` of `weights`.
"""
function flexleSamplers(weights::AbstractVector{Float64}, number::Int64)
    return [FlexleSampler(weights) for _ in 1:number]
end

"""
    FlexleSamplers(number)

Create a vector with separate, empty `FlexleSampler`.
"""
function flexleSamplers(number::Int64)
    return [FlexleSampler() for _ in 1:number]
end

sample(sampler::FlexleSampler) = Flexle.sample(sampler)
