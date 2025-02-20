using DataFrames
using Plots

function plotCompartments(
        output::Output, population_ids::Vector{String}, file_name::String;
        xlabel="Time", ylabel="Number",
        size=(800, 800), palette=:tab10, font_scale=2.0, grid=false,
        legend_location=:outerbottom, linewidth=3.0, thickness_scaling=1.0)

    compartment_vars = zeros(Float64, NUM_COMPARTMENTS, length(output.time))
    for id in population_ids
        compartment_vars .+= output.compartment_vars[id]
    end

    Plots.scalefontsizes(font_scale)
    plot(
        output.time, compartment_vars',
        labels=reshape([
            "Uninfected naive", "Uninfected immune",
            "Infected naive", "Infected immune", "Dead"
        ], 1, NUM_COMPARTMENTS),
        xlabel=xlabel, ylabel=ylabel,
        legend=legend_location, palette=palette,
        linewidth=linewidth, thickness_scaling=thickness_scaling,
        grid=grid,
    )
    plot!(size=size)

    savefig(file_name)
    Plots.scalefontsizes(1.0/font_scale)
end

function plotCompartments(
        output::DataFrame, population_ids::Vector{String}, file_name::String;
        xlabel="Time", ylabel="Number",
        size=(800, 800), palette=:tab10, font_scale=2.0, grid=false,
        legend_location=:outerbottom, linewidth=3.0, thickness_scaling=1.0)

    compartment_vars = zeros(Float64, NUM_COMPARTMENTS, length(output[!,"time"]))
    for id in population_ids
        compartment_vars .+= output[(output["Population"].==id),:]
    end

    Plots.scalefontsizes(font_scale)
    plot(
        output[:,"time"], compartment_vars',
        labels=reshape([
            "Uninfected naive", "Uninfected immune",
            "Infected naive", "Infected immune", "Dead"
        ], 1, NUM_COMPARTMENTS),
        xlabel=xlabel, ylabel=ylabel,
        legend=legend_location, palette=palette,
        linewidth=linewidth, thickness_scaling=thickness_scaling,
        grid=grid,
    )
    plot!(size=size)

    savefig(file_name)
    Plots.scalefontsizes(1.0/font_scale)
end
