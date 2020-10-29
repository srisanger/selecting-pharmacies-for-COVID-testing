"""
Utility and help functions for optimal distribution of antiviral drugs.

References:
Singh, B., Huang, H.C., Morton, D.P., Johnson, G.P., Gutfraind,A., Galvani, A.P., Clements, B., Meyers, L.A.: 
Optimizing distribution of pandemic influenza antiviral drugs. 
Emerging Infectious Diseases 21(2), 251–258 (2015). DOI 10.3201/eid2102.141024
"""

using JSON, CSV

function prop_willing_travel(dist,
                            dist_treshold=5,
                            K=[1, 1.479],
                            α=[-0.109, -0.4255],
                            β=[1.184, 0.6025])
    """
    Returns proportion of population willing to travel a distance 
    calculated from a piecewise decreasing exponential function.
    Default imputs from Singh et al. (2015, EID).
    """
    param_index = dist < dist_treshold ? 1 : 2
    return K[param_index] * exp(α[param_index] * dist^β[param_index])
end

function haversine_formula(lat1, lon1, lat2, lon2)
    """
    Returns distance in miles between two points in spherical coordinates (lat and lon in degrees).
    Calculated by Haversine formula, see http://www.movable-type.co.uk/scripts/gis-faq-5.1.html.
    NB: Formula assumes spherical Earth, not the actual ellipsoidal.
    """
    EARTH_R = 3956 # miles
    lat1, lon1, lat2, lon2 = (lat1, lon1, lat2, lon2) .* (π / 180)
    a = (sin((lat2 - lat1) / 2))^2 + cos(lat1) * cos(lat2) * (sin((lon2 - lon1) / 2))^2
    c = 2 * asin(min(1, sqrt(a)))
    return EARTH_R * c
end

function make_travel_frac(df_area::DataFrame, df_pharma::DataFrame, included_states; dist_limit=13.0574666)
    """
    Return dictionary of fraction willing to travel from an area to a pharmacy area.
    Assume distance sqrt(area)/2 and travel_frac = 1 for distances within same area.
    included states is an array of states included in the travel_frac.
    Returns area pairs, their distance apart and fraction willing to travel the distance.
    Deafult dist_limit = 13.0574666 correspond to 20% cut-off willingness to travel,
    with Singh et al. (2015, EID) willingness to travel model and parameters.
    """

    unique_pharma_areas = unique(df_pharma, :ZCTA)[:, [:STATE,:ZCTA]]
    pharma_loc = Dict(row[:ZCTA] => [row[:LAT], row[:LON]] for row in eachrow(filter(x -> x[:ZCTA] in unique_pharma_areas[:,:ZCTA], df_area)))
    included_states_areas = Dict(k => filter(x -> x[:STATE] in v, unique_pharma_areas)[:,:ZCTA] for (k, v) in included_states)

    distances = Dict()
    travel_frac = Dict()

    for (i, row) in enumerate(eachrow(df_area))
        lat_area = row[:LAT]
        lon_area = row[:LON]
        for area in included_states_areas[row[:STATE]]
            dist = haversine_formula(lat_area, lon_area, pharma_loc[area][1], pharma_loc[area][2])
            if dist < 0.001
                distances[row[:ZCTA], area] = sqrt(row[:LANDAREA]) / 2
                travel_frac[row[:ZCTA], area] = 1.0
            elseif dist < dist_limit
                travel_frac[row[:ZCTA], area] = prop_willing_travel(dist)
                distances[row[:ZCTA], area] = dist
            end
        end
    end
    return collect(keys(distances)), distances, travel_frac
end

function make_input_dict(df_zcta::DataFrame, df_pharma::DataFrame, pharma_chains, included_states; national=false, chain_combo=false)
    """
    Returns input dictionary from a dataframe
    pharma_chains is array of chain names included in study 
    included_states is an array of states included in travel fraction
    chain_combo binary denote whether run includes combinations of the largest chains: CVS, Walgreens and Walmart
    """
    if national
        input = Dict(
            "AREAS" => df_zcta[:,:ZCTA],
            "PHARMA_AREAS" => Dict("All" => unique(df_pharma[:,:ZCTA])),
            "POP" => Dict(row[:ZCTA] => row[:POP] for row in eachrow(df_zcta)),
        )
    else
        input = Dict(
            "AREAS" => filter(row -> row[:STATE] in keys(included_states), df_zcta)[:,:ZCTA],
            "PHARMA_AREAS" => Dict("All" => unique(filter(row -> row[:STATE] in keys(included_states), df_pharma)[:,:ZCTA])),
            "POP" => Dict(row[:ZCTA] => row[:POP] for row in eachrow(filter(row -> row[:STATE] in keys(included_states), df_zcta))),
        )
    end

    input["TRAVEL_PAIRS"], input["DISTANCE"], input["TRAVEL_FRAC"] = make_travel_frac(df_zcta, df_pharma, included_states)

    if pharma_chains !== nothing
        input["CHAINS"] = pharma_chains
        for chain in pharma_chains
            input["PHARMA_AREAS"][chain] = unique(filter(row -> row[:CHAIN] == chain, df_pharma)[:,:ZCTA])
        end
    end

    if chain_combo
        input["CHAINS"] = ["CVS & Walgreens", "CVS & Walmart", "Walgreens & Walmart", "CVS, Walgreens & Walmart"]
        input["PHARMA_AREAS"]["CVS & Walgreens"] = unique(filter(row -> row[:CHAIN] in ["CVS", "Walgreens"], df_pharma)[:,:ZCTA])
        input["PHARMA_AREAS"]["CVS & Walmart"] = unique(filter(row -> row[:CHAIN] in ["CVS", "Walmart"], df_pharma)[:,:ZCTA])
        input["PHARMA_AREAS"]["Walgreens & Walmart"] = unique(filter(row -> row[:CHAIN] in ["Walgreens", "Walmart"], df_pharma)[:,:ZCTA])
        input["PHARMA_AREAS"]["CVS, Walgreens & Walmart"] = unique(filter(row -> row[:CHAIN] in ["CVS", "Walgreens", "Walmart"], df_pharma)[:,:ZCTA])
    end

    return input
end

function find_largest_chains(all_loc_chains; max_include=5, min_locations=10)
    """Return array of chains with most locations"""
    loc_chains = filter(x -> x != "Other", all_loc_chains)
    loc_count = Dict()
    for chain in unique(loc_chains)
        loc_count[chain] = count(x -> x == chain, loc_chains)
    end
    sorted = sort(collect(loc_count), by=x -> x[2], rev=true)
    n_chains = max_include > length(loc_count) ? length(loc_count) : max_include
    largest_chains = []
    for i in 1:n_chains
        if sorted[i][2] >= min_locations
            largest_chains = [largest_chains; sorted[i][1]]
        end
    end
    return largest_chains
end

function query_results(model, set, param, df_area::DataFrame)
    """
    Finds results for pharmacy and area pairs that is cut off by distance treshold to limit problem size.
    Also referred to as the post-processing stage.
    """

    covered_pop = objective_value(model)
    closest = repeat([(0, 0)], length(set.TRAVEL_PAIRS))
    distances = Dict()

    for (i, pair) in enumerate(set.TRAVEL_PAIRS)
        if value(model[:closest_pharma][pair]) > 0.5
            closest[i] = pair
            distances[pair] = param.DISTANCE[pair]
        end
    end

    closest = filter(x -> x != (0, 0), closest)

    area_selected = map(x -> x[1], closest)
    pharma_selected = map(x -> x[2], closest)

    df_unselected_areas = filter(x -> x[:ZCTA] ∉ area_selected, df_area)
    df_selected_pharma = filter(x -> x[:ZCTA] in pharma_selected, df_area)

    closest_pharma = repeat([(0, 0)], size(df_unselected_areas)[1])

    for (i, row) in enumerate(eachrow(df_unselected_areas))
        dist, idx = findmin(haversine_formula.(row[:LAT], row[:LON], df_selected_pharma[:,:LAT], df_selected_pharma[:,:LON]))
        closest_pharma[i] = (row[:ZCTA], df_selected_pharma[idx,:ZCTA])
        distances[closest_pharma[i]] = dist
        covered_pop += row[:POP] * prop_willing_travel(dist)
    end

    closest = [closest; closest_pharma]

    return closest, covered_pop, distances
end

function get_closest_pharma(model, set)
    """Query model for closest_pharma variables that equals 1."""
    closest = []
    for i in set.AREAS, j in set.PHARMA_AREAS
        if value(model[:closest_pharma][i,j]) > 0.5
            closest = [closest; (i, j)]
        end
    end
    return closest
end

function get_select_pharma(model, pharma_areas)
    """Query model for selected_pharmas variables that equals 1."""
    selected_pharmas = []
    for j in pharma_areas
        if value(model[:select_pharma][j]) > 0.5
            selected_pharmas = [selected_pharmas; j]
        end
    end
    return selected_pharmas
end

function solve_full_models(input_dict)
    """Solves full model without cut-off for willingness to travel."""
    result_dict = Dict()
    for chain in keys(input_dict["PHARMA_AREAS"])
        set = sets(input_dict["AREAS"], input_dict["PHARMA_AREAS"][chain], nothing)
        param = parameters(input_dict["POP"], input_dict["TRAVEL_FRAC"], input_dict["DISTANCE"], 1)

        result_dict[chain] = Dict("N_PHARMA" => length(set.PHARMA_AREAS), "POP_COVER" => zeros(length(set.PHARMA_AREAS)), "CLOSEST" => Dict())

        model = make_full_model(set, param)

        for i in 1:length(set.PHARMA_AREAS)
            println("\n####################\n", chain, " at iteration ", i, "\n####################\n")
            optimize!(model)
            result_dict[chain]["POP_COVER"][i] = objective_value(model)
            result_dict[chain]["CLOSEST"][i] =  get_closest_pharma(model, set)
            set_normalized_rhs(model[:select_pharma_bound], i + 1)
        end
    end
    return result_dict
end

function get_runs(n_pharma)
    """Returns runs when gradually increasing maximum number of selected pharamcy areas"""
    if n_pharma < 500
        runs = 1:n_pharma
    elseif n_pharma < 1000
        runs = 1:5:n_pharma
    elseif n_pharma < 2000
        runs = 1:10:n_pharma
    elseif n_pharma < 4000
        runs = 1:20:n_pharma
    elseif n_pharma < 8000
        runs = 1:40:n_pharma
    elseif n_pharma < 16000
        runs = 1:80:n_pharma
    else
        runs = 1:100:n_pharma
    end
    runs = collect(runs)
    runs[2:end] = runs[2:end] .- 1
    runs = [runs; n_pharma]
    return runs
end

function solve_models(input_dict, df_area)
    """
    Formulates and solve several instances, where instances are different maximum number of selected pharamcy areas decided by get_runs function.
    """
    result_dict = Dict()
    for chain in keys(input_dict["PHARMA_AREAS"])

        set = sets(input_dict["AREAS"], input_dict["PHARMA_AREAS"][chain], filter(x -> x[2] in input_dict["PHARMA_AREAS"][chain], input_dict["TRAVEL_PAIRS"]))
        param = parameters(input_dict["POP"], input_dict["TRAVEL_FRAC"], input_dict["DISTANCE"], 1)

        runs = get_runs(length(set.PHARMA_AREAS))

        result_dict[chain] = Dict("N_PHARMA" => length(set.PHARMA_AREAS), "POP_COVER" => zeros(length(runs)), "CLOSEST" => Dict(), "DISTANCES" => Dict(), "RUNS" => runs)

        model = make_reduced_model(set, param)

        for (i, run) in enumerate(runs)
            println("\n####################\n", chain, " at iteration ", run, "\n####################\n")
            set_normalized_rhs(model[:select_pharma_bound], run)
            optimize!(model)
            result_dict[chain]["CLOSEST"][run], result_dict[chain]["POP_COVER"][i], result_dict[chain]["DISTANCES"][run] = query_results(model, set, param, df_area)
        end
    end
    return result_dict
end

function make_single_reduced_model(input_dict, chain, max_pharma)
    """Returns reduced optimization model from input_dict"""
    set = sets(input_dict["AREAS"], input_dict["PHARMA_AREAS"][chain], filter(x -> x[2] in input_dict["PHARMA_AREAS"][chain], input_dict["TRAVEL_PAIRS"]))
    param = parameters(input_dict["POP"], input_dict["TRAVEL_FRAC"], input_dict["DISTANCE"], max_pharma)
    model = make_reduced_model(set, param)
    return model
end

function include_distance_and_travel_prop(results, input_dict)
    """Includes distances and estimated population that travels in the result dictionary."""
    for chain in keys(results)
        results[chain]["DISTANCES"] = Dict()
        results[chain]["TRAVEL_POP"] = Dict()
        for i in 1:results[chain]["N_PHARMA"]
            results[chain]["DISTANCES"][i] = Dict()
            results[chain]["TRAVEL_POP"][i] = Dict()
            for pair in results[chain]["CLOSEST"][i]
                results[chain]["DISTANCES"][i][pair] = input_dict["DISTANCE"][pair]
                results[chain]["TRAVEL_POP"][i][pair[1]] = input_dict["POP"][pair[1]] * input_dict["TRAVEL_FRAC"][pair]
            end
        end
    end
    return results
end

function save_results(result_dict, fpath)
    """Saves full result dictionary as a json file"""
    open(fpath * "results.json", "w") do f
        JSON.print(f, result_dict)
    end
end

function save_pop_cover_results(result_dict, pop_tot, fpath)
    """Save estimated population covered as a json file"""
    pop_cover_results = Dict()
    pop_cover_results["POP_TOT"] = pop_tot
    for chain in keys(result_dict)
        pop_cover_results[chain] = result_dict[chain]["POP_COVER"]
    end
    open(fpath * "pop_cover.json", "w") do f
        JSON.print(f, pop_cover_results)
    end
end

function save_distance_results(result_dict, n_areas, fpath, percentiles=[0, 0.1, 0.5, 0.9, 1])
    """Save different percentiles of distances the full population must travel to their closest pharmacy."""
    for chain in keys(result_dict)
        n_runs = length(result_dict[chain]["RUNS"])
        distances = zeros(n_areas, n_runs)
        for (i, key) in enumerate(result_dict[chain]["RUNS"])
            distances[:,i] = collect(values(result_dict[chain]["DISTANCES"][key]))
        end

        df_results = DataFrame()
        for p in 1:length(percentiles)
            percs = zeros(n_runs)
            for i in 1:n_runs
                percs[i] = quantile(distances[:,i], percentiles[p])
            end
            df_results[!,Symbol("p" * string(percentiles[p]))] = percs
        end
        CSV.write(fpath * "travel_distances_" * lowercase(chain) * ".csv", df_results)
    end
end

function save_pop_distance_results(result_dict, input_dict, fpath, percentiles=[0, 0.1, 0.5, 0.9, 1], only_willing=false)
    """Save different percentiles of distances the population that decides to travel must travel to their closest pharmacy."""

    if !only_willing
        population = input_dict["POP"]
        percentile_pop = percentiles * sum(values(input_dict["POP"]))
    end

    for chain in keys(result_dict)
        plot_data = zeros(length(percentiles), length(result_dict[chain]["RUNS"]))

        for (i, key) in enumerate(result_dict[chain]["RUNS"])
            if only_willing
                population = result_dict[chain]["TRAVEL_POP"][key]
                percentile_pop = percentiles * sum(values(result_dict[chain]["TRAVEL_POP"][key]))
            end
            sorted_dist = PriorityQueue(result_dict[chain]["DISTANCES"][key])
            sorted_pairs = collect(keys(sorted_dist))
            pop_cumsum = cumsum([population[pair[1]] for pair in sorted_pairs])
            for p in 1:length(percentiles)
                perc_idx = searchsortedfirst(pop_cumsum, percentile_pop[p]) - 1
                perc_idx = perc_idx == 0 ? 1 : perc_idx
                perc_pair = sorted_pairs[perc_idx]
                plot_data[p,i] = sorted_dist[perc_pair]
            end
        end

        df_results = DataFrame()
        for (idx, p) in enumerate(percentiles)
            df_results[!, Symbol("p" * string(p))] = plot_data[idx,:]
        end

        fname = "population_travel_distances_" * lowercase(chain)

        if only_willing
            fname = "willing_" * fname
        end

        CSV.write(fpath * fname * ".csv", df_results)
    end
end

function query_results_from_selected(pharma_selected, df_area)
    """Query results from selected pharmas"""
    covered_pop = 0.0
    distances = Dict()
    df_selected_pharma = filter(x -> x[:ZCTA] in pharma_selected, df_area)

    closest_pharma = repeat([(0, 0)], size(df_area)[1])

    for (i, row) in enumerate(eachrow(df_area))
        dist, idx = findmin(haversine_formula.(row[:LAT], row[:LON], df_selected_pharma[:,:LAT], df_selected_pharma[:,:LON]))
        closest_pharma[i] = (row[:ZCTA], df_selected_pharma[idx,:ZCTA])

        if dist < 0.001
            distances[closest_pharma[i]] = sqrt(row[:LANDAREA]) / 2
            covered_pop += row[:POP]
        else
            distances[closest_pharma[i]] = dist
            covered_pop += row[:POP] * prop_willing_travel(dist)
        end
    end

    return closest_pharma, covered_pop, distances
end

