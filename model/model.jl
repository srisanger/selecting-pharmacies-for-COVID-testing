"""
Model from Singh et al. (2015, EID) to optimize distribution of antiviral drugs.
"""

using JuMP, GLPK

struct sets
    AREAS # i - area
    PHARMA_AREAS # j - area with pharmacy
    TRAVEL_PAIRS # (i,j)
end

struct parameters
    POP # Population at area i
    TRAVEL_FRAC # Fraction of population in area i willing to travel to area j
    DISTANCE # Distance between pairs (i,j)
    DEFAULT_MAX_PHARMAS # Hyperparameter that bound selected areas
end

function init_sets(input_dict)
    AREAS = input_dict["AREAS"]
    PHARMA_AREAS = input_dict["PHARMA_AREAS"]
    TRAVEL_PAIRS = input_dict["TRAVEL_PAIRS"]
    return sets(AREAS, PHARMA_AREAS, TRAVEL_PAIRS)
end

function init_parameters(input_dict, max_pharmas::Int64)
    POP = input_dict["POP"]
    TRAVEL_FRAC = input_dict["TRAVEL_FRAC"]
    DISTANCE = input_dict["DISTANCE"]
    DEFAULT_MAX_PHARMAS = max_pharmas
    return parameters(POP, TRAVEL_FRAC, DISTANCE, DEFAULT_MAX_PHARMAS)
end

function make_full_model(sets, param)
    model = Model(with_optimizer(GLPK.Optimizer))

    @variable(model, select_pharma[sets.PHARMA_AREAS], Bin)
    @variable(model, closest_pharma[sets.AREAS, sets.PHARMA_AREAS], Bin)

    @objective(model, Max, sum(param.POP[i] * param.TRAVEL_FRAC[i,j] * closest_pharma[i,j] for i in sets.AREAS, j in sets.PHARMA_AREAS))

    @constraint(model, select_pharma_bound, sum(select_pharma[j] for j in sets.PHARMA_AREAS) <= param.DEFAULT_MAX_PHARMAS)
    @constraint(model, find_closest_pharma[i in sets.AREAS], sum(closest_pharma[i,j] for j in sets.PHARMA_AREAS) == 1)
    @constraint(model, check_selected_pharma[i in sets.AREAS, j in sets.PHARMA_AREAS], closest_pharma[i,j] <= select_pharma[j])

    return model
end

function make_reduced_model(sets, param)
    """Solve model where travel fractions less than a certain threshold is omitted"""
    model = Model(with_optimizer(GLPK.Optimizer))
    
    @variable(model, select_pharma[sets.PHARMA_AREAS], Bin)
    @variable(model, closest_pharma[sets.TRAVEL_PAIRS], Bin)

    @objective(model, Max, sum(param.POP[i] * param.TRAVEL_FRAC[i,j] * closest_pharma[(i, j)] for (i, j) in sets.TRAVEL_PAIRS))

    @constraint(model, select_pharma_bound, sum(select_pharma[j] for j in sets.PHARMA_AREAS) <= param.DEFAULT_MAX_PHARMAS)

    travel_pair_areas = unique(map(x -> x[1], sets.TRAVEL_PAIRS))
    travel_pair_pharmas = Dict(i => Int64[] for i in travel_pair_areas)
    for (i, j) in sets.TRAVEL_PAIRS
        travel_pair_pharmas[i] = [travel_pair_pharmas[i]; j]
    end

    @constraint(model, find_closest_pharma[i in travel_pair_areas], sum(closest_pharma[(i, j)] for j in travel_pair_pharmas[i]) <= 1)
    @constraint(model, check_selected_pharma[(i, j) in sets.TRAVEL_PAIRS], closest_pharma[(i, j)] <= select_pharma[j])

    return model
end

