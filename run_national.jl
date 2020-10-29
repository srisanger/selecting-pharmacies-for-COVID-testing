"""
Module to run case study of pharmacy cover nationally
"""

using CSV, DataFrames, JSON
include("./model/utils.jl");
include("./model/model.jl");
include("./model/plotting.jl");

# To reduce computation, in particular for e.g. national runs of the whole US, the 
# make_input_dict fuction only calculates travel distances to states in INCLUDED_STATES, 
# which is usually the state being calculated and its neighbors. 
# Here, State ONE adjoins TWO, state TWO adjoins ONE and THREE, while state THREE adjoins TWO.
INCLUDED_STATES = Dict(
"ONE" => ["ONE","TWO"],
"TWO" => ["ONE", "TWO", "THREE"],
"THREE" => ["TWO", "THREE"]);

FPATH_IN = "./data/"
FPATH_OUT = "./results/national/"
df_areas = CSV.read(FPATH_IN * "/areas.csv");
df_areas = filter(x -> x[:STATE] in keys(INCLUDED_STATES), df_areas);
df_pharma = CSV.read(FPATH_IN * "/pharmacies.csv");
df_pharma = filter(x -> x[:STATE] in keys(INCLUDED_STATES), df_pharma);

population = sum(df_areas[:,:POP])
chains = find_largest_chains(df_pharma[:,:CHAIN])


input_dict = make_input_dict(df_areas, df_pharma, chains, INCLUDED_STATES, national=true)
results = solve_models(input_dict, df_areas)

save_pop_cover_results(results, population, FPATH_OUT)
save_distance_results(results, length(input_dict["AREAS"]), FPATH_OUT)
save_pop_distance_results(results, input_dict, FPATH_OUT)

save_pop_cover_fig(results, population, FPATH_OUT)
save_distance_fig(results, length(input_dict["AREAS"]), FPATH_OUT)
save_pop_distance_fig(results, input_dict, FPATH_OUT)
