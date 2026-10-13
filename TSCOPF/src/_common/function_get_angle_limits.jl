# =====================================================================
# ----------------------------- Auxiliar ------------------------------
# =====================================================================

# Function to standardize the min and max angular deviation between adjacent buses
function get_angle_limits(
                          DCIR::DataFrame, 
                          active_branches::Vector{Int64},
    )

    # Initialize dictionary with count, min_ang, max_ang
    pair_info = Dict{Tuple{Int, Int}, NamedTuple{(:count, :min_ang, :max_ang), Tuple{Int, Float64, Float64}}}()
    pair_circ_map = Dict{Tuple{Int, Int}, Int}()

    # This loop maps the circuits that have parallel branches, count the number of parallel lines/transformers,
    # and save the minimum and maximum angular difference for the buses related to these circuits
    for branch in active_branches
        a, b = DCIR.from_bus[branch], DCIR.to_bus[branch]  # Bus from and bus to
        pair = (min(a, b), max(a, b))                      # Sort the buses id
        circ_num = DCIR.id[branch]                         # Get the number of the branch


        if !haskey(pair_circ_map, pair)     # Check if these buses were already addded
            pair_circ_map[pair] = circ_num  # Get only the number of the first branch ON connecting these buses
        else
            pair_circ_map[pair] = min(pair_circ_map[pair], circ_num) # Get the number of the first branch ON connecting these buses
        end
        
        if haskey(pair_info, pair) # Check if this pair of buses were already added
            # Update count
            old       = pair_info[pair]
            new_count = old.count + 1                          # Count the number of branches ON connecting these buses
            new_min   = max(old.min_ang, DCIR.ang_min[branch]) # Get the minimum angle defined for the branches connecting these buses
            new_max   = min(old.max_ang, DCIR.ang_max[branch]) # Get the maximum angle defined for the branches connecting these buses
            pair_info[pair] = (new_count, new_min, new_max)
        else
            pair_info[pair] = (1, DCIR.ang_min[branch], DCIR.ang_max[branch]) # Add this pair of buses for the first time
        end
    end
    sorted_pair_info  = sort(collect(pair_info), by = x -> pair_circ_map[x[1]]) # Sort the data inside pair_info Dict
    sorted_pair_circ_map = sort(collect(pair_circ_map), by = x -> x.first) # Sort the data inside pair_circ_map Dict

    return sorted_pair_info, sorted_pair_circ_map
end