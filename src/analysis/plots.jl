using DataFrames
using Plots
using NewickTree

function plotCompartments(
    output::DataFrame, population_ids::Vector{String}, file_name::String;
    xlabel="Time", ylabel="Number",
    dimensions=(800, 800), palette=:tab10, font_scale=2.0, grid=false,
    legend_location=:outerbottom, linewidth=3.0, thickness_scaling=1.0)

    compartment_vars = zeros(Float64, length(output[:, "Time"]), NUM_COMPARTMENTS)
    for id in population_ids
        compartment_vars .+= output[(output[:, "Population"].==id), 3:end]
    end

    Plots.scalefontsizes(font_scale)
    plot(
        output[:, "Time"], compartment_vars,
        labels=reshape([
                "Uninfected naive", "Uninfected immune",
                "Infected naive", "Infected immune", "Dead"
            ], 1, NUM_COMPARTMENTS),
        xlabel=xlabel, ylabel=ylabel,
        legend=legend_location, palette=palette,
        linewidth=linewidth, thickness_scaling=thickness_scaling,
        grid=grid,
    )
    plot!(size=dimensions)

    savefig(file_name)
    Plots.scalefontsizes(1.0 / font_scale)
end


function plotComposition(
    output::DataFrame, file_name::String;
    normalized=false,
    xlabel="Time", ylabel="Number", legend=true,
    dimensions=(800, 800), palette=:tab10, font_scale=2.0, grid=false,
    legend_location=:outerbottom, linewidth=0.0, thickness_scaling=1.0)

    dat = DataFrame()
    dat[:, "Time"] = output[:, "Time"]
    if normalized
        for col in names(output)[2:end-1]
            dat[:, col] = output[:, col] ./ max.(output[:, "Total"], 1.0)
        end
        dat[:, "Total"] = output[:, "Total"]
    else
        dat = deepcopy(output)
    end

    if legend
        legend_content = reshape(names(dat)[2:end-1], 1, length(names(dat)) - 2)
    else
        legend_content = false
    end

    Plots.scalefontsizes(font_scale)
    areaplot(
        dat[:, "Time"], Matrix(dat[:, 2:end-1]),
        labels=legend_content,
        xlabel=xlabel, ylabel=ylabel,
        # seriescolor=palette, fillalpha = 1.0,
        legend=legend_location, palette=palette,
        linewidth=linewidth, thickness_scaling=thickness_scaling,
        grid=grid,
    )
    plot!(size=dimensions)

    savefig(file_name)
    Plots.scalefontsizes(1.0 / font_scale)
end

function plotPhylogeny(
    newick_string::String, file_name::String;
    dimensions=(1600, 800), palette=:tab10, font_scale=2.0, grid=false,
    linewidth=2.0, thickness_scaling=1.0)
    Plots.scalefontsizes(font_scale)
    plot(
        readnw(newick_string), transform=true,
        palette=palette,
        linewidth=linewidth, thickness_scaling=thickness_scaling,
        grid=grid,
    )
    plot!(size=dimensions)

    savefig(file_name)
    Plots.scalefontsizes(1.0 / font_scale)
end
