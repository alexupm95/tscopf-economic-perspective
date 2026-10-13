# =====================================================================================
#  Parametric Jacobian ∂c/∂p for parameters that are NOT decision variables
# =====================================================================================
#
#  Demand is data in this codebase, never a JuMP object: `DBUS.p_d` is inlined as a
#  Float64 (FULL_BUS ZIP right-hand side) or as a ComplexF64 matrix entry (the Kron
#  constant-impedance shunt, which then propagates through `Reduce_Matrix` into every
#  entry of `Y_red`). `build_constraint_evaluator` differentiates with respect to
#  `all_variables(model)`, so there is no `P_d` column and MOI cannot supply ∂c/∂P_d.
#
#  This file supplies it WITHOUT hand-differentiating anything: rebuild the model with
#  the load perturbed, evaluate the perturbed model's own constraint functions at the
#  BASE solution x*, and difference.
#
#      ∂c_r/∂p ≈ [ c_r(x*; p+ε) − c_r(x*; p−ε) ] / (2ε),     c_r(x) := g_r(x) − rhs_r
#
#  `Y_red`, the faulted-bus diagonal overwrite and the GL scaling chain rule all come
#  along for free, because the perturbed build runs the same builders the real model
#  does. Structural discovery is a side effect: any row with |∂c_r/∂p| > tol IS an
#  entry point, so no family-name list exists to go stale.
#
#  TRAP — differentiate `g − rhs`, never `g` alone. Load can sit in the set constant
#  as easily as in the function; on the ACOPF balance it is inside the function
#  (`Pg_sum - p_d - Σterms == 0`), but a builder rewrite could move it. Dropping the
#  rhs term would silently zero the ENERGY channel — the largest one — and report the
#  price as 100 % stability.
#
#  TRAP — row rescaling. `functions_4_TS_fullbus_eqconst.jl:112,122` writes a DIFFERENT
#  algebraic form at buses with |p_d| < 1e-12: `terms_Pb == inj_P` instead of
#  `terms_Pb*V² == inj_P*V² − p_d*ZIP(V)`. Perturbing such a bus off zero multiplies the
#  row by V[i]², and the dual scales by 1/V[i]², so a naive difference is wrong by that
#  factor. On the 9-bus case this is SIX of nine buses, not an edge case. `rownorm`
#  below detects it generically (any uniform rescale, not just V²); `param_column`
#  divides it out.
#
#  Everything here runs inside `RunConfig.post_solve_hook`. `release_solver_backend!`
#  empties the model immediately afterwards, so no ConstraintRef survives the return.
# =====================================================================================

using JuMP
using CSV, DataFrames, DataStructures, SparseArrays, LinearAlgebra
using TSCOPF
import MathOptInterface as MOI
import TSCOPF: build_constraint_evaluator

# `label_jacobian`, `_set_rhs` and `_walk_store` live in kkt_sensitivity.jl. The driver
# includes that file first (every existing example already does), so they are in scope.

"""
    Save_Param_Residuals!(path_names, opf_dict, dyn_model_dict, model;
                          tag, x_ref_path = nothing) -> String

Evaluate every constraint of THIS model at a reference primal point and dump the
residuals, one row per Jacobian row, to `Sensitivity/resid_<tag>.csv`.

`tag` names the leg: `"base"` for the unperturbed run, `"bus07_P_plus"` and friends for
the perturbed ones.

`x_ref_path === nothing` marks the base leg: the model's own solution is the reference,
and `Sensitivity/x_ref.csv` is written alongside so the perturbed legs can find it.
Otherwise `x_ref_path` points at that file and its values are mapped onto this model's
columns.

The mapping is positional and asserted, not assumed. A load perturbation changes numeric
coefficients only, so the variable vector of the perturbed model must match the base one
element for element; if it does not, the builder took a different structural branch and
every difference computed downstream would be meaningless.

Columns written:

    row;source;family;kind;key1;key2;rhs;dual;resid;rownorm

`resid` is `g_r(x_ref) − rhs_r`. `rownorm` is `‖Jac_r(x_ref)‖₁`, the rescale detector.
"""
function Save_Param_Residuals!(path_names, opf_dict, dyn_model_dict, model;
                               tag::String, x_ref_path = nothing)::String

    evaluator, rows, xvars = build_constraint_evaluator(model; features = [:Jac])
    row_of = Dict(rows[i] => i for i in eachindex(rows))
    col_of = Dict(xvars[j] => j for j in eachindex(xvars))
    rows_df, cols_df = label_jacobian(opf_dict, dyn_model_dict, rows, xvars, row_of, col_of)

    out_dir = joinpath(path_names[:pf_results_date], "Sensitivity")
    mkpath(out_dir)

    # --- the reference primal point ------------------------------------------------
    if x_ref_path === nothing
        x_ref = JuMP.value.(xvars)
        CSV.write(joinpath(out_dir, "x_ref.csv"),
                  select(cols_df, :col, :source, :family, :key1, :key2, :name, :value);
                  delim = ';')
    else
        ref = CSV.read(x_ref_path, DataFrame; delim = ';')
        nrow(ref) == nrow(cols_df) || error(
            "x_ref has $(nrow(ref)) columns but this model has $(nrow(cols_df)) — " *
            "the perturbed build is not structurally identical to the base.")
        # Positional identity, asserted on the labels rather than trusted.
        mismatch = findfirst(i -> (ref.family[i], ref.key1[i], ref.key2[i]) !=
                                  (cols_df.family[i], cols_df.key1[i], cols_df.key2[i]),
                             1:nrow(ref))
        mismatch === nothing || error(
            "variable ordering diverged at column $mismatch: base " *
            "$(ref.family[mismatch])[$(ref.key1[mismatch]),$(ref.key2[mismatch])] vs " *
            "$(cols_df.family[mismatch])[$(cols_df.key1[mismatch]),$(cols_df.key2[mismatch])]")
        x_ref = Vector{Float64}(ref.value)
    end

    # --- constraint values at that point --------------------------------------------
    # `g − rhs`, not `g`: see the header trap.
    g = zeros(length(rows))
    MOI.eval_constraint(evaluator, g, x_ref)
    resid = g .- rows_df.rhs

    # --- row norms, the rescale detector ---------------------------------------------
    # Accumulate ‖Jac_r‖₁ straight from the COO triplets, in ONE pass.
    #
    # Do NOT assemble a SparseMatrixCSC and slice rows out of it: CSC stores columns,
    # so `view(Jac, r, :)` rescans the whole structure for every row and the loop
    # becomes O(n_rows · nnz) — ~1e10 operations on a 30k-row DQ Jacobian, which
    # hangs long enough to look like a solver problem rather than a bug.
    sparsity = MOI.jacobian_structure(evaluator)
    V = zeros(length(sparsity))
    MOI.eval_constraint_jacobian(evaluator, V, x_ref)
    rownorm = zeros(length(rows))
    for (k, (r, _c)) in enumerate(sparsity)
        rownorm[r] += abs(V[k])
    end

    out = DataFrame(
        row = rows_df.row, source = rows_df.source, family = rows_df.family,
        kind = rows_df.kind, key1 = rows_df.key1, key2 = rows_df.key2,
        rhs = rows_df.rhs, dual = rows_df.dual, resid = resid, rownorm = rownorm,
    )
    file_path = joinpath(out_dir, "resid_$(tag).csv")
    CSV.write(file_path, out; delim = ';')
    return file_path
end

"""
    param_column(base_dir, plus_dir, minus_dir; ε_pu, tol = 1e-12,
                 base_tag = "base", plus_tag = "plus", minus_tag = "minus") -> DataFrame

Central-difference one parametric column from three exported residual dumps and return
the rows that actually depend on the parameter. The three `*_tag` arguments name the
`resid_<tag>.csv` files written by `Save_Param_Residuals!`; the three directories may all
be the same run folder if the legs were dumped side by side.

`ε_pu` is the HALF-perturbation in per unit — the same units the constraints are stated
in, so the returned derivative needs no rescaling to land in the €/MW convention the
generator-side Tier 1 already uses.

Returned columns:

    source;family;kind;key1;key2;dual;dcdp;scale;scale_plus;scale_minus

`dcdp` is `∂c_r/∂p`, divided by the symmetrised `scale` so a row the builder rewrote at
a different scale contributes its true coefficient. `scale ≈ 1` on every ordinary row —
a free assertion that nothing else got renormalised — and `≈ V[i]²` exactly on the rows
of a bus whose base load is zero. `scale_plus`/`scale_minus` are kept unsymmetrised for
diagnosis: their spread is the O(ε) coefficient drift, and it should shrink with ε while
`scale` does not.

Rows with `|dcdp| ≤ tol` are dropped: what survives IS the census of where the parameter
enters the built model.
"""
function param_column(base_dir::String, plus_dir::String, minus_dir::String;
                      ε_pu::Float64, tol::Float64 = 1e-12,
                      base_tag::String = "base", plus_tag::String = "plus",
                      minus_tag::String = "minus")::DataFrame

    rd(dir, tag) = CSV.read(joinpath(dir, "Sensitivity", "resid_$(tag).csv"),
                            DataFrame; delim = ';')
    b = rd(base_dir, base_tag)
    p = rd(plus_dir, plus_tag)
    m = rd(minus_dir, minus_tag)

    key(df, i) = (df.source[i], df.family[i], df.kind[i], df.key1[i], df.key2[i])
    nrow(b) == nrow(p) == nrow(m) || error(
        "row counts differ: base $(nrow(b)), plus $(nrow(p)), minus $(nrow(m)) — " *
        "a perturbed build took a different structural branch.")
    bad = findfirst(i -> !(key(b, i) == key(p, i) == key(m, i)), 1:nrow(b))
    bad === nothing || error("row ordering diverged at row $bad ($(key(b, bad)))")

    # A row whose base norm is zero carries no information and cannot be rescaled.
    safe(x) = iszero(x) ? 1.0 : x
    scale_p = p.rownorm ./ safe.(b.rownorm)
    scale_m = m.rownorm ./ safe.(b.rownorm)

    # SYMMETRISE, do not use either leg's ratio on its own.
    #
    # Two different things move a row's norm between base and perturbed. (1) The row's
    # own coefficients depend on p_d, so every row drifts by O(ε) — antisymmetric,
    # s₊ ≈ 1+kε and s₋ ≈ 1−kε. (2) At a bus with |p_d| < 1e-12 the builder switches
    # algebraic branch and multiplies the whole row by V[i]² — symmetric, present in
    # both legs. Only (2) is a real rescale that the dual compensates for; dividing by
    # (1) injects an O(ε/p_d) error, ~0.3 % at 0.5 MW on a 187.5 MW bus.
    #
    # The mean kills the antisymmetric part and keeps the symmetric one. Magnitude alone
    # could never separate them: V² ≈ 1 in per unit, so a genuine branch change and pure
    # drift look identical on one leg.
    scale = (scale_p .+ scale_m) ./ 2

    dcdp = ((p.resid .- m.resid) ./ safe.(scale)) ./ (2 * ε_pu)

    out = DataFrame(
        source = b.source, family = b.family, kind = b.kind,
        key1 = b.key1, key2 = b.key2, dual = b.dual,
        dcdp = dcdp, scale = scale, scale_plus = scale_p, scale_minus = scale_m,
    )
    return out[abs.(out.dcdp) .> tol, :]
end
