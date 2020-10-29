"""
Module to run case study of pharmacy cover of the 50 states
"""

using CSV, DataFrames, JSON
include("./model/utils.jl");
include("./model/model.jl");
include("./model/plotting.jl");

FPATH_IN = "./data/"
FPATH_OUT = "./results/states/"
df_areas = CSV.read(FPATH_IN * "areas.csv");
df_pharma = CSV.read(FPATH_IN * "pharmacies.csv");

STATES = Dict("State One" => "ONE", "State Two" => "TWO", "State Three" => "THREE");

for state in values(STATES)
    println("\n####################\n", state, "\n####################\n")

    df_state_areas = filter(row -> row[:STATE] == state, df_areas);
    df_state_pharma = filter(row -> row[:STATE] == state, df_pharma);
    state_pop = sum(df_state_areas[:,:POP])
    chains = find_largest_chains(df_state_pharma[:,:CHAIN])
    input_dict = make_input_dict(df_state_areas, df_state_pharma, chains, Dict(state => [state]))
    results = solve_models(input_dict, df_state_areas)

    fpath_results = FPATH_OUT * state * "/"
    save_pop_cover_results(results, state_pop, fpath_results)
    save_distance_results(results, length(input_dict["AREAS"]), fpath_results)
    save_pop_distance_results(results, input_dict, fpath_results)
    save_pop_cover_fig(results, state_pop, fpath_results)
    save_distance_fig(results, length(input_dict["AREAS"]), fpath_results)
    save_pop_distance_fig(results, input_dict, fpath_results)
end
