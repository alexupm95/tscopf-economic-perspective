# Auxiliar functions to retrieve the Jacobian, Hessian matrices of the solved optimization model

# Only performing some tests trying to recover the Jacobian and Hessian matrices from the optimization problem
# More details in the JuMP documentation at: https://jump.dev/JuMP.jl/dev/tutorials/nonlinear/querying_hessians/#Jacobians

# Function to calculate and save the matrices built from the optimization model.
# Everything is written as sparse COO triplets (row;col;value) — only the structural
# nonzeros — into RESULTS/.../Optim_Matrices/, with index→name legends so the numeric
# row/column indices stay interpretable. This is gated by save_optim_matrices, which is
# forced false for TSC runs (see resolve_save_optim_matrices) because these matrices have
# thousands of rows/columns there.
function Compute_and_Save_Optimal_Matrices(model::Model, path_names::OrderedDict{Symbol, String}; return_matrices::Bool = false)

    out_dir = joinpath(path_names[:pf_results_date], "Optim_Matrices")
    mkpath(out_dir)

    # Decision-variable order — the column index of the Jacobian/Hessian and the row
    # index of the Hessian/gradient all refer to this ordering.
    save_index_legend_csv(joinpath(out_dir, "Variable_Index.csv"),
                          name.(all_variables(model)), :variable)

    # Jacobian of the constraints (sparse), plus its row legend (constraint expressions,
    # in the SAME loop order the Jacobian rows are built).
    jacobian_matrix = compute_optimal_jacobian(model)
    save_sparse_coo_csv(jacobian_matrix, joinpath(out_dir, "Jacobian_COO.csv"))
    save_index_legend_csv(joinpath(out_dir, "Jacobian_Constraint_Index.csv"),
                          constraint_labels(model), :constraint)

    # Hessian of the Lagrangian (sparse, symmetric; rows and cols indexed by Variable_Index).
    hessian_matrix = compute_optimal_hessian(model)
    save_sparse_coo_csv(hessian_matrix, joinpath(out_dir, "Hessian_COO.csv"))

    # Gradient of the Lagrangian (dense vector, indexed by Variable_Index).
    grad_lagrange = compute_lagrangian_gradient(model)
    save_gradient_to_csv(model, grad_lagrange, path_names)

    println("Optimal matrices (sparse COO) saved to: ", out_dir)

    if return_matrices
        return jacobian_matrix, hessian_matrix, grad_lagrange
    else
        return nothing
    end

end

"""
    build_constraint_evaluator(model; features=[:Jac]) -> (evaluator, rows, x)

Nonlinear evaluator for every constraint of `model` that is not a plain variable
bound, together with the two orderings the derivative matrices are indexed by:

- `rows::Vector{ConstraintRef}` — Jacobian **row** order,
- `x::Vector{VariableRef}` — Jacobian **column** (and Hessian row/column) order.

Returning the orderings is the point. The row order is an artefact of iterating
`list_of_constraint_types` and then `all_constraints`, so without it an exported
Jacobian can only be matched back to the model by parsing the constraint strings in
`Jacobian_Constraint_Index.csv`. With it, `Dict(rows .=> eachindex(rows))` maps any
`ConstraintRef` straight to its row, which is what any analysis of the KKT system
needs.

`features` is forwarded to `MOI.initialize`: `[:Jac]`, `[:Hess]`, `[:Grad, :Jac]`.
"""
function build_constraint_evaluator(model::Model; features::Vector{Symbol} = [:Jac])
    rows = ConstraintRef[]
    nlp = MOI.Nonlinear.Model()
    for (F, S) in list_of_constraint_types(model)
        F <: VariableRef && continue      # variable bounds carry their duals separately
        for ci in all_constraints(model, F, S)
            push!(rows, ci)
            object = constraint_object(ci)
            MOI.Nonlinear.add_constraint(nlp, object.func, object.set)
        end
    end
    MOI.Nonlinear.set_objective(nlp, objective_function(model))
    x = all_variables(model)
    evaluator = MOI.Nonlinear.Evaluator(nlp, MOI.Nonlinear.SparseReverseMode(), index.(x))
    MOI.initialize(evaluator, features)
    return evaluator, rows, x
end

# Function to compute the Jacobian at the optimal solution
function compute_optimal_jacobian(model::Model)
    evaluator, rows, x = build_constraint_evaluator(model; features = [:Jac])
    # Query the Jacobian structure
    sparsity = MOI.jacobian_structure(evaluator)
    I, J, V = first.(sparsity), last.(sparsity), zeros(length(sparsity))
    # Query the Jacobian values
    MOI.eval_constraint_jacobian(evaluator, V, value.(x))

    return SparseArrays.sparse(I, J, V, length(rows), length(x))
end

# Function to fill off diagonal elements of a sparse matrix
function fill_off_diagonal(H)
    ret = H + H'
    row_vals = SparseArrays.rowvals(ret)
    non_zeros = SparseArrays.nonzeros(ret)
    for col in 1:size(ret, 2)
        for i in SparseArrays.nzrange(ret, col)
            if col == row_vals[i]
                non_zeros[i] /= 2
            end
        end
    end
    return ret
end

# Function to compute the Hessian at the optimal solution
function compute_optimal_hessian(model::Model)
    evaluator, rows, x = build_constraint_evaluator(model; features = [:Hess])
    hessian_sparsity = MOI.hessian_lagrangian_structure(evaluator)
    I = [i for (i, _) in hessian_sparsity]
    J = [j for (_, j) in hessian_sparsity]
    V = zeros(length(hessian_sparsity))
    MOI.eval_hessian_lagrangian(evaluator, V, value.(x), 1.0, dual.(rows))
    H = SparseArrays.sparse(I, J, V, length(x), length(x))

    return fill_off_diagonal(H) # sparse, symmetric (off-diagonals mirrored); saved as COO
end

# Function to compute the Gradient at the optimal solution
function compute_lagrangian_gradient(model::Model)
    # 1. Setup the NLP evaluator, initialized for both gradient and Jacobian
    evaluator, rows, x_vars = build_constraint_evaluator(model; features = [:Grad, :Jac])
    x_val = value.(x_vars)

    # 3. Compute Objective Gradient: ∇f(x)
    grad_f = zeros(length(x_vars))
    MOI.eval_objective_gradient(evaluator, grad_f, x_val)
    
    # 4. Compute Jacobian: J(x)
    sparsity = MOI.jacobian_structure(evaluator)
    I, J, V = first.(sparsity), last.(sparsity), zeros(length(sparsity))
    MOI.eval_constraint_jacobian(evaluator, V, x_val)
    Jac = SparseArrays.sparse(I, J, V, length(rows), length(x_vars))
    
    # 5. Get Dual Variables (Multipliers): λ
    # Note: Check sign conventions; JuMP duals for inequalities might need a sign flip 
    # depending on the formulation (f + λg vs f - λg).
    λ = dual.(rows)
    
    # 6. Calculate ∇L = ∇f - J' * λ
    grad_L = grad_f - Jac' * λ
    
    return grad_L
end

# Function to show all the constraints in the terminal
function show_constraints(model::Model)
    # List all types of constraints in the model
    for (F, S) in list_of_constraint_types(model) 
        println("Constraints of type $F in $S:")
        for con in all_constraints(model, F, S)
            println("  ", constraint_object(con))
        end
    end
end

# Write a sparse matrix as COO triplets (row;col;value) — only the structural nonzeros.
# Scales with nnz, not with n², so it stays small for the sparse Jacobian/Hessian and
# loads in one line elsewhere (Julia readdlm, Python scipy.io / pandas, MATLAB).
function save_sparse_coo_csv(M::SparseMatrixCSC, filepath::String)
    I, J, V = SparseArrays.findnz(M)
    CSV.write(filepath, DataFrame(row = I, col = J, value = V); delim=";")
    return filepath
end

# Write an index→name legend so the numeric COO row/column indices are interpretable.
function save_index_legend_csv(filepath::String, labels::Vector{String}, colname::Symbol)
    CSV.write(filepath, DataFrame(:index => collect(1:length(labels)), colname => labels); delim=";")
    return filepath
end

# Constraint expression strings in the SAME loop order used to build the Jacobian rows
# (compute_optimal_jacobian iterates list_of_constraint_types skipping VariableRef bounds).
function constraint_labels(model::Model)
    labels = String[]
    for (F, S) in list_of_constraint_types(model)
        for ci in all_constraints(model, F, S)
            if !(F <: VariableRef)
                push!(labels, string(constraint_object(ci)))
            end
        end
    end
    return labels
end

# Function to save the Gradient of the Lagrangian (dense vector) next to the COO matrices.
function save_gradient_to_csv(model::Model, grad_vector::Vector{Float64}, path_names::OrderedDict{Symbol, String})
    out_dir = joinpath(path_names[:pf_results_date], "Optim_Matrices")
    mkpath(out_dir)
    df = DataFrame(Variable = name.(all_variables(model)), Gradient_Value = grad_vector)
    CSV.write(joinpath(out_dir, "Gradient_Lagrangian.csv"), df; delim=";")
    println("Gradient of Lagrangian saved to: ", out_dir)
    return nothing
end

using DataFrames, CSV, SparseArrays

"""
Saves a report of constraints, their current values, and their duals.
Also decomposes the gradient of the Lagrangian for a specific variable.
"""
function save_constraint_diagnostics(model::Model, filename::String = "constraint_duals.csv")
    # 1. Collect all non-variable constraints
    all_cons = Any[]
    for (F, S) in list_of_constraint_types(model)
        # Skip simple variable bounds (Pg >= Pmin) as they are handled by dual.(VariableRef)
        if !(F <: VariableRef)
            append!(all_cons, all_constraints(model, F, S))
        end
    end

    # 2. Extract Data
    df = DataFrame(
        Constraint_Name = String[],
        Constraint_Type = String[],
        Dual_Value      = Float64[],
        Slack           = Float64[],
        Is_Active       = Bool[]
    )

    for ci in all_cons
        name  = name(ci) == "" ? string(constraint_object(ci).func) : name(ci)
        c_dual = dual(ci)
        # Slack tells you how far from binding the constraint is
        c_slack = JuMP.slack(ci) 
        
        push!(df, (
            name, 
            string(typeof(ci)), 
            c_dual, 
            c_slack, 
            abs(c_slack) < 1e-6 # True if binding
        ))
    end

    # 3. Save to CSV
    CSV.write(filename, df)
    println("Constraint diagnostics saved to $filename")
    
    return df
end

"""
Decomposes the Stationarity of a specific variable (e.g., Pg[2]).
This shows you exactly which constraint 'dual' is balancing the cost.
"""
function analyze_variable_stationarity(model::Model, target_var::VariableRef)
    # Get objective gradient component for this variable (e.g., c_g)
    # For linear objective sum(c_g * Pg), it's just the coefficient
    obj_grad = derivative(objective_function(model), target_var)
    
    println("\n--- Stationarity Analysis for $(name(target_var)) ---")
    println("Objective Marginal Cost (df/dx): $obj_grad")
    
    # Check all constraints that involve this variable
    total_constraint_contribution = 0.0
    
    for (F, S) in list_of_constraint_types(model)
        for ci in all_constraints(model, F, S)
            # Get the function part of the constraint
            func = constraint_object(ci).func
            
            # Check if our variable is in this constraint
            # This works for Linear/Quadratic. For Nonlinear, use the Evaluator approach.
            coeff = 0.0
            try
                coeff = derivative(func, target_var)
            catch
                continue # Skip if variable not present
            end
            
            if coeff != 0.0
                c_dual = dual(ci)
                contribution = c_dual * coeff
                total_constraint_contribution += contribution
                
                if abs(contribution) > 1e-5
                    println("Constraint: $(name(ci)) | Dual: $c_dual | Coeff: $coeff | Contrib: $contribution")
                end
            end
        end
    end
    
    # Don't forget Variable Bounds!
    lower_dual = has_lower_bound(target_var) ? reduced_cost(target_var) : 0.0
    
    println("Variable Bound Dual (Reduced Cost): $lower_dual")
    
    residual = obj_grad + total_constraint_contribution + lower_dual
    println("--------------------------------------------------")
    println("LAGRANGIAN GRADIENT RESIDUAL: $residual")
end


"""
Iterates through all constraints and variable bounds, 
saving the symbolic expression and the dual value to a CSV.
"""
function save_constraints_to_csv(model::Model, path_names::OrderedDict{Symbol, String}; filename::String = "model_constraints_duals.csv")
    out_path = joinpath(path_names[:pf_results_date], filename)
    mkpath(path_names[:pf_results_date])
    # Initialize a DataFrame to store our data
    df = DataFrame(Expression = String[], Dual_Value = Float64[])

    # 1. Capture Linear, Quadratic, and Nonlinear Constraints
    for (F, S) in list_of_constraint_types(model)
        # Skip VariableRef here, we handle bounds in step 2
        if !(F <: VariableRef)
            for ci in all_constraints(model, F, S)
                # Get the symbolic function (e.g., P_g - P_elec)
                expr_str = string(constraint_object(ci).func)
                # Get the dual multiplier (λ)
                dual_val = dual(ci)
                
                push!(df, (expr_str, dual_val))
            end
        end
    end

    # 2. Capture Variable Bounds (Pg >= Pmin, etc.)
    for v in all_variables(model)
        if has_lower_bound(v)
            push!(df, ("Lower Bound: $(name(v))", extract_var_side_dual(v, :lower)))
        end
        if has_upper_bound(v)
            push!(df, ("Upper Bound: $(name(v))", extract_var_side_dual(v, :upper)))
        end
    end

    # 3. Write to file
    CSV.write(out_path, df; delim=";")
    println("Successfully saved $(nrow(df)) constraints/bounds to ", out_path)
end

# Usage:
# save_constraints_to_csv(model, "Lagrangian_Components.csv")