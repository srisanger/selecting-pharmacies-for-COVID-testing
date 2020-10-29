"""
Module containing functions for plotting
"""

using Plots, Statistics, DataStructures, DataFrames

function save_pop_cover_fig(result_dict, pop_tot, fpath)
    x_lab = "Number of ZIP codes with pharmacies"
    plot()
    for chain in sort(collect(keys(result_dict)))
        plot!(result_dict[chain]["RUNS"], result_dict[chain]["POP_COVER"] ./ pop_tot .* 100,
        label=chain,
        xlabel=x_lab, 
        ylabel="% of population with access")
    end
    plot!(ylims=(0, 100), legend=:right)

    savefig(fpath * "pop_cover.pdf")
end

function save_distance_fig(result_dict, n_areas, fpath, percentiles=[0.1, 0.5, 0.9])
    x_lab = "Number of ZIP codes with pharmacies"
    y_lab = "Distance to closest pharmacy [miles]"
    for chain in sort(collect(keys(result_dict)))
        n_runs = length(result_dict[chain]["RUNS"])
        distances = zeros(n_areas, n_runs)
        for (i, key) in enumerate(result_dict[chain]["RUNS"])
            distances[:,i] = collect(values(result_dict[chain]["DISTANCES"][key]))
        end

        plot_data = zeros(length(percentiles), n_runs)
        for p in 1:length(percentiles)
            for i in 1:n_runs
                plot_data[p,i] = quantile(distances[:,i], percentiles[p])
            end
        end

        plot()
        plot!(result_dict[chain]["RUNS"], plot_data[1,:], label="10 percentile")
        plot!(result_dict[chain]["RUNS"], plot_data[2,:], label="Median")
        plot!(result_dict[chain]["RUNS"], plot_data[3,:], label="90 percentile")
        plot!(title=chain,
            xlabel=x_lab,
            ylabel=y_lab,
            ylim=[0,50])
        savefig(fpath * "travel_distances_" * lowercase(chain) * ".pdf")
    end
end

function save_pop_distance_fig(result_dict, input_dict, fpath, percentiles=[0.1, 0.5, 0.9], only_willing=false)
    x_lab = "Number of ZIP codes with pharmacies"
    y_lab = "Distance to closest pharmacy [miles]"
    if !only_willing
        population = input_dict["POP"]
        percentile_pop = percentiles * sum(values(input_dict["POP"]))
    end

    for chain in sort(collect(keys(result_dict)))
        plot_data = zeros(length(percentiles), length(result_dict[chain]["RUNS"]))

        for (i, key) in enumerate(result_dict[chain]["RUNS"])
            if only_willing
                population = result_dict[chain]["TRAVEL_POP"][key]
                percentile_pop = percentiles * sum(values(result_dict[chain]["TRAVEL_POP"][key]))
            end
            sorted_dist = PriorityQueue(result_dict[chain]["DISTANCES"][key])
            sorted_pairs = collect(keys(sorted_dist))
            pop_cumsum = cumsum([population[pair[1]] for pair in sorted_pairs])
            for j in 1:length(percentiles)
                perc_idx = searchsortedfirst(pop_cumsum, percentile_pop[j]) - 1
                perc_idx = perc_idx == 0 ? 1 : perc_idx
                perc_pair = sorted_pairs[perc_idx]
                plot_data[j,i] = sorted_dist[perc_pair]
            end
        end

        plot()
        plot!(result_dict[chain]["RUNS"], plot_data[1,:], label="10 percentile")
        plot!(result_dict[chain]["RUNS"], plot_data[2,:], label="Median")
        plot!(result_dict[chain]["RUNS"], plot_data[3,:], label="90 percentile")
        plot!(title=chain,
            xlabel=x_lab,
            ylabel=y_lab,
            ylim=[0,50])
        fname = "population_travel_distances_" * lowercase(chain)
        if only_willing
            fname = "willing_" * fname
        end
        savefig(fpath * fname * ".pdf")
    end
end