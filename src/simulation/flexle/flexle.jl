# Flexible binary level rejection sampling

module Flexle

using Revise
using Printf
using Random
using BenchmarkTools

include("sampler.jl")
include("interface.jl")

export FlexleSampler, update!, push!, deleteat!, getindex, setindex!, sample

end
