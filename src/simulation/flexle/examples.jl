using .Flexle
using Plots
using StatsPlots
using Revise

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
    print(s)
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
