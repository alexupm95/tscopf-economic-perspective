#=
================================================================================
 kkt_sensitivity.jl — Tier 2: S = ∂h/∂P from the KKT matrices
================================================================================
 Not part of the TSCOPF module; loaded by `runs/run_39bus.jl` (J_KKT = 2) and
 reached through `RunConfig.post_solve_hook`.

 With h_g^t := δ_g^t − δ_COI^t (what the corridor bounds), the implicit function
 theorem on the dynamic equality block F(x, P) = 0 gives

     dx/dP = −F_x⁻¹ F_P ,     S_P = h_P + h_x·dx/dP

 Stationarity in x is NOT just F_xᵀα + h_xᵀγ̃ = 0. The speed corridor h_ω, and any
 other active inequality on a state variable — an optional `bound_P_m` / `bound_E`
 box, the FULL_BUS voltage floor — enter too:

     F_xᵀα + h_xᵀγ̃ + h_ω,xᵀγ̃_ω + Σ_b b_xᵀμ_b = 0        (γ̃ = Γ⁺ − Γ⁻)

 which turns the causer-pays identity into a SUM over channels:

     γ̃ᵀS_P + γ̃_ωᵀS_ω,P + Σ_b μ_bᵀS_b,P
         = γ̃ᵀh_P + γ̃_ωᵀh_ω,P + αᵀF_P + Σ_b μ_bᵀb_P

 δ and Δω are separate CHANNELS, not one h: h_ω bounds a speed deviation rather than
 an angle, so S_ω is a different matrix with different units. γ̃ is in EUR/rad and γ̃_ω
 in EUR/(pu speed), but both charges land in EUR/MW — so the split into angle-corridor
 rent and speed-corridor rent is a genuine decomposition in one unit. = −dyn

 with S_b,P = b_P + b_x·dx/dP, the same construction as S_P against the same
 factorisation. Dropping the Σ_b term does not produce a visible error — it produces a
 plausible number that is short by exactly μ_bᵀS_b,P (measured: a 179-bus run with
 `bound_P_m = true` came up short by that machine's own P_m multiplier, on 15 of 29
 machines, while every other diagnostic looked clean).

 So the blocks are discovered STRUCTURALLY — every row that touches a state column and
 is not already F or h — rather than from a list of family names. A name list is the
 thing that goes stale when a builder gains a box, which is how the P_m leak survived.

 `dyn` is exactly the Tier-1 dynamic contribution to ∂L/∂P_g. h_P = 0 whenever the
 corridor is a plain box on δ; it is nonzero on `:swing_propagated`, which states the
 corridor in substituted form and so carries P_g directly.

 Every derivative here comes from `build_constraint_evaluator`, i.e. from MOI's
 automatic differentiation of the actual model. Nothing is hand-derived, so this
 file does not need editing when a builder changes, and it extends to FULL_BUS /
 DQ / GFM by adding rows to the registry below — not by adding calculus.
================================================================================
=#

using JuMP
using CSV, DataFrames, DataStructures, SparseArrays, LinearAlgebra, TOML
using TSCOPF
import MathOptInterface as MOI
import TSCOPF: build_constraint_evaluator

# --- which rows are F, which columns are x -----------------------------------------
#
# Mirrors `dual_spec_catalog` (DynDualRegistry.jl): the path is decided by
# meta[:gen_order] / meta[:network_form], and each path names the families that make
# up the square system. Adding a path = adding one entry here.

struct DynStateSpec
    label::String
    F::Vector{Symbol}      # equality families defining the trajectory
    x::Vector{Symbol}      # state variable families
    h::Vector{Symbol}      # δ corridor rows, UPPER side (the lower side is −h)
    h_lower::Vector{Symbol}
    h_ω::Vector{Symbol}    # Δω corridor rows, same convention
    h_ω_lower::Vector{Symbol}
end

const KRON_CLASSICAL_STATE = DynStateSpec(
    "KRON_REDUCED/CLASSICAL_2ND",
    # `eq_const_ΔωCOI_*` and `ΔωCOI_*` must be added or omitted TOGETHER — one row and
    # one column per time step each, so listing both keeps n_F == n_x. They exist only
    # under `constrain_Δω` with `bound_style_Δω = :coi_box`; absent families contribute
    # nothing, which is what makes listing them unconditionally safe (the same mechanism
    # that already covers δCOI when it degrades to an expression).
    [:eq_const_P_init, :eq_const_Q_init, :eq_const_Pm_init,
     :eq_const_Pe_tf, :eq_const_δ_tf, :eq_const_Δω_tf, :eq_const_δCOI_tf,
     :eq_const_Pe_tpf, :eq_const_δ_tpf, :eq_const_Δω_tpf, :eq_const_δCOI_tpf,
     :eq_const_ΔωCOI_tf, :eq_const_ΔωCOI_tpf],
    [:E, :P_m, :δ,
     :Pe_tf, :δ_tf, :Δω_tf, :δCOI_tf,
     :Pe_tpf, :δ_tpf, :Δω_tpf, :δCOI_tpf,
     :ΔωCOI_tf, :ΔωCOI_tpf],
    # Both δ-corridor flavours. `:coi_box` / `:swing_propagated` store under the
    # δ_COI keys; `:highest_H` / `:ref_gen` store under δ_ref. Only one flavour
    # exists in a given run, and absent families contribute no rows, so listing
    # both means the machine-referenced corridor works without a second spec.
    #
    [:ineq_const_δ_COI_tf_upper, :ineq_const_δ_COI_tpf_upper,
     :ineq_const_δ_ref_tf_upper, :ineq_const_δ_ref_tpf_upper],
    [:ineq_const_δ_COI_tf_lower, :ineq_const_δ_COI_tpf_lower,
     :ineq_const_δ_ref_tf_lower, :ineq_const_δ_ref_tpf_lower],
    # The speed corridor is its OWN channel, not more rows of the δ one: h_ω bounds a
    # speed deviation rather than an angle, so S_ω is a different matrix with different
    # units. Both flavours listed — `:coi_box` stores under Δω_COI, `:abs` under Δω_abs,
    # and only one exists per run.
    #
    # Left as a discovered x-block it would still close the identity (adjoint and charge
    # use the same slice) but report the wrong functional under `:coi_box`, because the
    # ΔωCOI column would be dropped and S_ω would measure ABSOLUTE speed.
    [:ineq_const_Δω_COI_tf_upper, :ineq_const_Δω_COI_tpf_upper,
     :ineq_const_Δω_abs_tf_upper, :ineq_const_Δω_abs_tpf_upper],
    [:ineq_const_Δω_COI_tf_lower, :ineq_const_Δω_COI_tpf_lower,
     :ineq_const_Δω_abs_tf_lower, :ineq_const_Δω_abs_tpf_lower],
)

const FULLBUS_CLASSICAL_STATE = DynStateSpec(
    "FULL_BUS/CLASSICAL_2ND",
    # Same three pre-fault links as Kron, then per window: the machine's own P/Q
    # injection, the NODAL balances (which is what FULL_BUS adds — Kron folds the
    # network into Y_red and has no bus equations), and the two swing rows.
    [:eq_const_P_init, :eq_const_Q_init, :eq_const_Pm_init,
     :eq_const_Pe_tf, :eq_const_Qe_tf,
     :eq_const_Pbalance_tf, :eq_const_Qbalance_tf,
     :eq_const_δ_tf, :eq_const_Δω_tf,
     :eq_const_δCOI_tf, :eq_const_ΔωCOI_tf,
     :eq_const_Pe_tpf, :eq_const_Qe_tpf,
     :eq_const_Pbalance_tpf, :eq_const_Qbalance_tpf,
     :eq_const_δ_tpf, :eq_const_Δω_tpf,
     :eq_const_δCOI_tpf, :eq_const_ΔωCOI_tpf],
    # The matching columns. `V_tf` / `θ_tf` are the DYNAMIC bus voltages, one per bus
    # per step, paired with the nodal balances; the steady-state `V` / `θ` live in
    # opf_dict under different names and stay outside x, as on Kron.
    [:E, :P_m, :δ,
     :V_tf, :θ_tf, :Pe_tf, :Qe_tf, :δ_tf, :Δω_tf, :δCOI_tf, :ΔωCOI_tf,
     :V_tpf, :θ_tpf, :Pe_tpf, :Qe_tpf, :δ_tpf, :Δω_tpf, :δCOI_tpf, :ΔωCOI_tpf],
    # Corridors are network-form agnostic: `_add_δ_bounds_*!` and `_add_Δω_bounds!` are
    # shared by every path, so the family names are identical to Kron.
    [:ineq_const_δ_COI_tf_upper, :ineq_const_δ_COI_tpf_upper,
     :ineq_const_δ_ref_tf_upper, :ineq_const_δ_ref_tpf_upper],
    [:ineq_const_δ_COI_tf_lower, :ineq_const_δ_COI_tpf_lower,
     :ineq_const_δ_ref_tf_lower, :ineq_const_δ_ref_tpf_lower],
    [:ineq_const_Δω_COI_tf_upper, :ineq_const_Δω_COI_tpf_upper,
     :ineq_const_Δω_abs_tf_upper, :ineq_const_Δω_abs_tpf_upper],
    [:ineq_const_Δω_COI_tf_lower, :ineq_const_Δω_COI_tpf_lower,
     :ineq_const_Δω_abs_tf_lower, :ineq_const_Δω_abs_tpf_lower],
)

const DQ_FULLBUS_STATE = DynStateSpec(
    "FULL_BUS/DQ_4TH",
    # Pre-fault is 7 rows per machine, matching the 7 columns below. `Vd_init` and
    # `Vq_init` are EXPRESSIONS, not variables (functions_2_build_TS_model_w_DqFullBus.jl
    # :396-398), so their equalities pin `Id`/`Iq` rather than adding columns — the
    # counting only balances once you notice that.
    [:eq_const_P_init, :eq_const_Q_init, :eq_const_Pm_init,
     :eq_const_Ed_init, :eq_const_Eq_init, :eq_const_Vd_init, :eq_const_Vq_init,
     # Per window: the four dq machine rows, the torque, the machine's P/Q injection,
     # the nodal balances, and the two swing rows. Same shape in both windows.
     :eq_const_Ed_tf, :eq_const_Eq_tf, :eq_const_Vd_tf, :eq_const_Vq_tf,
     :eq_const_Te_tf, :eq_const_Pe_tf, :eq_const_Qe_tf,
     :eq_const_Pbalance_tf, :eq_const_Qbalance_tf,
     :eq_const_δ_tf, :eq_const_Δω_tf, :eq_const_δCOI_tf, :eq_const_ΔωCOI_tf,
     :eq_const_Ed_tpf, :eq_const_Eq_tpf, :eq_const_Vd_tpf, :eq_const_Vq_tpf,
     :eq_const_Te_tpf, :eq_const_Pe_tpf, :eq_const_Qe_tpf,
     :eq_const_Pbalance_tpf, :eq_const_Qbalance_tpf,
     :eq_const_δ_tpf, :eq_const_Δω_tpf, :eq_const_δCOI_tpf, :eq_const_ΔωCOI_tpf,
     # AVR (include_avr): V_ref per machine, then the lead-lag, the unsaturated field
     # voltage and its smooth saturation per machine per step.
     :eq_const_Efd_init,
     :eq_const_avr_leadlag_tf, :eq_const_E_fd_unlim_tf, :eq_const_E_fd_tf,
     :eq_const_avr_leadlag_tpf, :eq_const_E_fd_unlim_tpf, :eq_const_E_fd_tpf,
     # Governor (include_governor): P_ref per machine, then the valve and mechanical
     # rows. Under GOV_SMOOTH the valve ODE state and its clamped output are separate
     # variables; under GOV_NO_LIMIT / GOV_HARD_BOUND they are the SAME container, so
     # one of the two column families simply comes back empty.
     :eq_const_Pref_init,
     :eq_const_gov_valve_tf, :eq_const_gov_mech_tf,
     :eq_const_gov_valve_tpf, :eq_const_gov_mech_tpf,
     # GOV_SMOOTH only: Pv = smoothclip(Pv_raw), the row that pins the clamped output
     # against the ODE state. Stored under a dynamically built key
     # (functions_4_TS_governor.jl:522), so it does not appear in a grep for literal
     # `[:eq_const][:...]` keys — omitting it left n_x exceeding n_F by exactly the
     # Pv_raw column count. Absent under GOV_NO_LIMIT / GOV_HARD_BOUND, where Pv and
     # Pv_raw are the same container and no clamp row exists.
     :eq_const_gov_valve_limit_tf, :eq_const_gov_valve_limit_tpf],
    # `E_fd` is the scalar field voltage and is built on EVERY dq run, not only with an
    # AVR (:354). The AVR instead adds the time-indexed `E_fd_tf/tpf` and their
    # saturation equalities, which is why controls need their own listing.
    [:E_fd, :δ, :Ed, :Eq, :Id, :Iq, :P_m,
     :V_tf, :θ_tf, :Ed_tf, :Eq_tf, :Id_tf, :Iq_tf, :Te_tf,
     :Pe_tf, :Qe_tf, :δ_tf, :Δω_tf, :δCOI_tf, :ΔωCOI_tf,
     :V_tpf, :θ_tpf, :Ed_tpf, :Eq_tpf, :Id_tpf, :Iq_tpf, :Te_tpf,
     :Pe_tpf, :Qe_tpf, :δ_tpf, :Δω_tpf, :δCOI_tpf, :ΔωCOI_tpf,
     # AVR columns matching the rows above. `E_fd` is NOT here — it is already listed
     # as a pre-fault dq state, because it exists with or without an AVR.
     :V_ref, :E_LL_tf, :E_fd_unlim_tf, :E_fd_tf,
     :E_LL_tpf, :E_fd_unlim_tpf, :E_fd_tpf,
     # Governor columns. `P_m` is already listed for the same reason as `E_fd`.
     :P_ref, :Pv_raw_tf, :Pv_tf, :Pm_tf,
     :Pv_raw_tpf, :Pv_tpf, :Pm_tpf],
    [:ineq_const_δ_COI_tf_upper, :ineq_const_δ_COI_tpf_upper,
     :ineq_const_δ_ref_tf_upper, :ineq_const_δ_ref_tpf_upper],
    [:ineq_const_δ_COI_tf_lower, :ineq_const_δ_COI_tpf_lower,
     :ineq_const_δ_ref_tf_lower, :ineq_const_δ_ref_tpf_lower],
    [:ineq_const_Δω_COI_tf_upper, :ineq_const_Δω_COI_tpf_upper,
     :ineq_const_Δω_abs_tf_upper, :ineq_const_Δω_abs_tpf_upper],
    [:ineq_const_Δω_COI_tf_lower, :ineq_const_Δω_COI_tpf_lower,
     :ineq_const_Δω_abs_tf_lower, :ineq_const_Δω_abs_tpf_lower],
)

"""
State-space spec for this run, or `nothing` when the path has no entry yet.

Neither the Kron builder nor the classical FULL_BUS one writes `meta[:gen_order]` —
only the DQ builder does — so an absent key means CLASSICAL_2ND. The Kron builder
writes no `:network_form` either, hence the default. Both defaults are load-bearing;
do not "tidy" them away without checking what each builder actually stores.

`ineq_const_V_tf_lower`, the always-on FULL_BUS voltage floor, needs no entry here: it
is found by `x_block_rows` and charged as a block. Enabling the governor or AVR DOES
need work — they add state families (`Pv_*`, `Pm_*`, `E_fd_*`) and their equalities,
so `square` will go false until those are listed.
"""
function state_spec(dyn_model_dict)
    meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
    gen_order = get(meta, :gen_order, "CLASSICAL_2ND")
    network_form = get(meta, :network_form, "KRON_REDUCED")
    if gen_order == "CLASSICAL_2ND"
        network_form == "KRON_REDUCED" && return KRON_CLASSICAL_STATE
        network_form == "FULL_BUS" && return FULLBUS_CLASSICAL_STATE
    elseif gen_order == "DQ_4TH"
        # DQ exists on FULL_BUS only; the enum permits nothing else.
        network_form == "FULL_BUS" && return DQ_FULLBUS_STATE
    end
    return nothing
end

# --- labelling the Jacobian --------------------------------------------------------

"""
Walk a `[:vars]` / `[:eq_const]` / `[:ineq_const]` store and yield one row per
JuMP reference, with the container keys attached.

Handles both shapes used across this codebase: `id → ref` and `id → t → ref`. The
first key is a generator id in most families, a bus id in the nodal ones and a time
index in `δCOI_*`; it is recorded verbatim and interpreted by whoever reads it.
"""
function _walk_store(store, source::String, kind::String, sink::Function)
    for (family, container) in store
        container isa AbstractDict || continue
        for (k1, v1) in container
            if v1 isa AbstractDict
                for (k2, v2) in v1
                    sink(String(family), source, kind, Int(k1), Int(k2), v2)
                end
            else
                sink(String(family), source, kind, Int(k1), -1, v1)
            end
        end
    end
    return nothing
end

"""
    label_jacobian(opf_dict, dyn_model_dict, row_of, col_of) -> (rows_df, cols_df)

One table per Jacobian axis, labelled with the model's own family names. References
the walker cannot place keep `family = "UNMAPPED"` so a missing family is visible in
the export rather than silently absent from the slices.
"""
function label_jacobian(opf_dict, dyn_model_dict, rows, xvars, row_of, col_of)
    row_family = fill("UNMAPPED", length(rows))
    row_source = fill("", length(rows))
    row_kind = fill("", length(rows))
    row_k1 = fill(-1, length(rows))
    row_k2 = fill(-1, length(rows))

    function claim_row(family, source, kind, k1, k2, ref)
        i = get(row_of, ref, nothing)
        i === nothing && return nothing         # variable-bound rows are not in the Jacobian
        row_family[i] = family; row_source[i] = source
        row_kind[i] = kind; row_k1[i] = k1; row_k2[i] = k2
        return nothing
    end

    for (dict, source) in ((opf_dict, "dispatch"), (dyn_model_dict, "dynamic"))
        dict === nothing && continue
        haskey(dict, :eq_const) && _walk_store(dict[:eq_const], source, "eq", claim_row)
        haskey(dict, :ineq_const) && _walk_store(dict[:ineq_const], source, "ineq", claim_row)
    end

    rows_df = DataFrame(
        row = collect(1:length(rows)),
        source = row_source, family = row_family, kind = row_kind,
        key1 = row_k1, key2 = row_k2,
        dual = [JuMP.dual(c) for c in rows],
        rhs = [_set_rhs(constraint_object(c).set) for c in rows],
    )

    col_family = fill("UNMAPPED", length(xvars))
    col_source = fill("", length(xvars))
    col_k1 = fill(-1, length(xvars))
    col_k2 = fill(-1, length(xvars))

    function claim_col(family, source, _kind, k1, k2, ref)
        j = get(col_of, ref, nothing)
        j === nothing && return nothing
        col_family[j] = family; col_source[j] = source
        col_k1[j] = k1; col_k2[j] = k2
        return nothing
    end

    for (dict, source) in ((opf_dict, "dispatch"), (dyn_model_dict, "dynamic"))
        dict === nothing && continue
        haskey(dict, :vars) && _walk_store(dict[:vars], source, "var", claim_col)
    end

    cols_df = DataFrame(
        col = collect(1:length(xvars)),
        source = col_source, family = col_family,
        key1 = col_k1, key2 = col_k2,
        name = JuMP.name.(xvars),
        value = JuMP.value.(xvars),
    )

    return rows_df, cols_df
end

"""
Make a value writable by `TOML.print`.

TOML has no `nan` or `inf` literal, and `TOML.print` throws on them rather than
degrading — which would lose the entire meta file over one unset diagnostic. Any
non-finite float becomes its string form instead.
"""
_toml_safe(x::AbstractFloat) = isfinite(x) ? Float64(x) : string(x)
_toml_safe(x::AbstractDict) = OrderedDict{String, Any}(string(k) => _toml_safe(v) for (k, v) in x)
_toml_safe(x::AbstractVector) = [_toml_safe(v) for v in x]
_toml_safe(x) = x

_set_rhs(s::MOI.EqualTo) = s.value
_set_rhs(s::MOI.LessThan) = s.upper
_set_rhs(s::MOI.GreaterThan) = s.lower
_set_rhs(::Any) = NaN

"""Row (or column) indices of `families`, in Jacobian index order."""
function _select(df::DataFrame, families::Vector{Symbol}, idcol::Symbol)
    wanted = Set(String.(families))
    return [row[idcol] for row in eachrow(df) if row.family in wanted]
end

"""
    corridor_row_pairs(rows_df, upper_families, lower_families) -> (upper_rows, lower_rows)

Upper corridor rows and the lower row of the SAME `(generator, step)`, aligned
element by element.

Takes the family lists rather than the spec so it serves both corridors: the δ box and
the Δω box use the same `id → t` container shape, so one implementation covers them.

Selecting each side independently and trusting the order would be wrong: JuMP groups
constraints by function/set type, so the upper and lower families need not be
interleaved in the same way, and a mismatch would silently pair generator 1's upper
row with generator 2's lower one — producing a `γ̃` that is nonsense but finite.
Pairing is therefore done on the container keys.
"""
function corridor_row_pairs(rows_df::DataFrame, upper_families::Vector{Symbol},
                            lower_families::Vector{Symbol})
    by_key = Dict{Tuple{String, Int, Int}, Int}()
    for r in eachrow(rows_df)
        by_key[(r.family, r.key1, r.key2)] = r.row
    end
    upper, lower = Int[], Int[]
    for (fam_up, fam_lo) in zip(upper_families, lower_families)
        keys_up = [(r.key1, r.key2) for r in eachrow(rows_df) if r.family == String(fam_up)]
        sort!(keys_up)
        for (k1, k2) in keys_up
            lo = get(by_key, (String(fam_lo), k1, k2), nothing)
            lo === nothing && throw(ArgumentError(
                "corridor upper row $(fam_up)[$k1][$k2] has no matching $(fam_lo) row"))
            push!(upper, by_key[(String(fam_up), k1, k2)])
            push!(lower, lo)
        end
    end
    return upper, lower
end

# --- extra constraint blocks on the state variables ---------------------------------

"""
    x_block_rows(Jac, claimed, x_cols) -> Vector{Int}

Every Jacobian row that touches at least one state column and is not already claimed by
`F` or by the corridor pair.

Structural discovery, deliberately not a list of family names. The identity needs *every*
active constraint on `x` in the ledger, and a name list is exactly the thing that goes
stale the moment a builder gains a box — which is how the `bound_P_m` leak survived
unnoticed through a full validation pass. Asking the Jacobian "what else touches x?"
cannot go stale, and it makes new model paths (FULL_BUS, DQ, GFM) need no block
registration at all: only the F / x / h lists.

Rows with no state entry are simply absent, so dispatch-side families (power balance, P_g
and Q_g boxes, line limits) drop out on their own rather than by being excluded by name.
"""
function x_block_rows(Jac::SparseMatrixCSC, claimed::Set{Int}, x_cols::Vector{Int})
    touched = falses(size(Jac, 1))
    for j in x_cols
        for k in nzrange(Jac, j)
            iszero(Jac.nzval[k]) && continue
            touched[Jac.rowval[k]] = true
        end
    end
    return [i for i in 1:size(Jac, 1) if touched[i] && !(i in claimed)]
end

"""
    block_charges(Jac, rows_df, block_rows, x_cols, P_cols, dxdP)
        -> (per_family::DataFrame, w_total::Vector, charge_total::Vector, bP_total::Vector)

Charge `μ_bᵀS_b,P` for each discovered block, grouped by constraint family, plus the
adjoint contribution `Σ_b b_xᵀμ_b` those blocks make and the direct term `Σ_b μ_bᵀb_P`
that appears on the other side of the identity.

`bP_total` is zero for every block that is a plain box on a state — which is all of them
today. It is computed rather than assumed so that a block which ever does carry `P_g`
shows up as a contribution instead of as an identity residual.

`S_b` is never materialised. The accounting only needs `μ_bᵀ(b_P + b_x·dx/dP)`, so the
`b_x` product is taken first and the result is a single row vector — necessary rather
than merely tidy, because `ineq_const_V_tf_lower` is one row per bus per time step
(≈250k on a 179-bus case) and the materialised `S_b` would be 250k×n_P.
"""
function block_charges(Jac::SparseMatrixCSC, rows_df::DataFrame, block_rows::Vector{Int},
                       x_cols::Vector{Int}, P_cols::Vector{Int}, dxdP::Matrix{Float64})
    n_P = length(P_cols)
    per_family = DataFrame(family = String[], kind = String[], n_rows = Int[],
                           max_abs_dual = Float64[], charge = Vector{Float64}[])
    w_total = zeros(length(x_cols))
    charge_total = zeros(n_P)
    bP_total = zeros(n_P)
    isempty(block_rows) && return per_family, w_total, charge_total, bP_total

    for fam in unique(rows_df.family[block_rows])
        rows_f = [r for r in block_rows if rows_df.family[r] == fam]
        μ_f = rows_df.dual[rows_f]
        B_x = Jac[rows_f, x_cols]
        B_P = Jac[rows_f, P_cols]
        # Both products stay sparse-times-vector. Densifying B_P would allocate
        # n_rows × n_P, which for the FULL_BUS voltage floor is ~250k × n_gen.
        w_f = B_x' * μ_f                       # n_x — this block's adjoint contribution
        bP_f = B_P' * μ_f                      # n_P — direct P_g content, normally zero
        charge_f = dxdP' * w_f .+ bP_f         # μ_fᵀS_f,P
        w_total .+= w_f
        charge_total .+= charge_f
        bP_total .+= bP_f
        push!(per_family, (fam, rows_df.kind[first(rows_f)], length(rows_f),
                           maximum(abs, μ_f), collect(charge_f)))
    end
    return per_family, w_total, charge_total, bP_total
end

# --- the export --------------------------------------------------------------------

"""
    Save_KKT_Sensitivity!(path_names, cfg, sys, opf_dict, dyn_model_dict, model) -> String

Write `<run>/Sensitivity/` and return the folder path. Must run while the model is
alive (the Jacobian is evaluated at the solution).
"""
function Save_KKT_Sensitivity!(path_names, cfg, sys, opf_dict, dyn_model_dict, model)::String
    spec = state_spec(dyn_model_dict)
    if spec === nothing
        meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
        @warn "No DynStateSpec for this model path — Tier 2 export skipped." *
              " Add one to kkt_sensitivity.jl." gen_order=get(meta, :gen_order, "?") network_form=get(meta, :network_form, "?")
        return ""
    end

    out_dir = joinpath(path_names[:pf_results_date], "Sensitivity")
    mkpath(out_dir)

    # --- the exact Jacobian at the solution, plus its two orderings ----------------
    evaluator, rows, xvars = build_constraint_evaluator(model; features = [:Jac])
    row_of = Dict(rows[i] => i for i in eachindex(rows))
    col_of = Dict(xvars[j] => j for j in eachindex(xvars))

    sparsity = MOI.jacobian_structure(evaluator)
    I, J, V = first.(sparsity), last.(sparsity), zeros(length(sparsity))
    x_val = JuMP.value.(xvars)
    MOI.eval_constraint_jacobian(evaluator, V, x_val)
    Jac = sparse(I, J, V, length(rows), length(xvars))

    rows_df, cols_df = label_jacobian(opf_dict, dyn_model_dict, rows, xvars, row_of, col_of)

    # --- slice out the blocks -----------------------------------------------------
    F_rows = _select(rows_df, spec.F, :row)
    h_rows, hlo_rows = corridor_row_pairs(rows_df, spec.h, spec.h_lower)
    hω_rows, hωlo_rows = corridor_row_pairs(rows_df, spec.h_ω, spec.h_ω_lower)
    x_cols = _select(cols_df, spec.x, :col)
    P_cols = _select(cols_df, [:P_g], :col)
    Q_cols = _select(cols_df, [:Q_g], :col)

    # Everything else that touches a state column. Discovered from the Jacobian, not
    # from a name list — see `x_block_rows`. Both corridors are claimed first so a
    # named channel is never also counted as a block.
    claimed = Set{Int}(vcat(F_rows, h_rows, hlo_rows, hω_rows, hωlo_rows))
    block_rows = x_block_rows(Jac, claimed, x_cols)

    # Stamp the partition onto the tables themselves, so a reader slices on `block`
    # instead of re-declaring the registry on the other side of the CSV.
    block_set = Set(block_rows)
    rows_df.block = [r.family in Set(String.(spec.F)) ? "F" :
                     r.family in Set(String.(spec.h)) ? "h" :
                     r.family in Set(String.(spec.h_lower)) ? "h_lower" :
                     r.family in Set(String.(spec.h_ω)) ? "h_omega" :
                     r.family in Set(String.(spec.h_ω_lower)) ? "h_omega_lower" :
                     r.row in block_set ? "x_block" : "other"
                     for r in eachrow(rows_df)]

    # An equality outside F that touches x is NOT a block — it means the state-space spec
    # is incomplete, and `dx/dP` is then wrong even though the identity below will still
    # close (the term is carried on both sides). Loud, because a closing identity would
    # otherwise make this invisible.
    stray_eq = unique(rows_df.family[[r for r in block_rows if rows_df.kind[r] == "eq"]])
    isempty(stray_eq) || @warn "Equality families touch the state block but are not in " *
        "spec.F — dx/dP is computed from an incomplete system. Add them to the spec." stray_eq

    # Under VARIABLE encoding simple boxes never become Jacobian rows, so a binding one
    # cannot be discovered here at all. The leak still happens; only the attribution is lost.
    #
    # Compared through `string`: meta holds the `BoundEncoding` ENUM, not a Symbol or a
    # String, so a direct `== :CONSTRAINT` silently never matches and warns on every run.
    let enc = get(get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}()),
                  :bound_encoding, "CONSTRAINT")
        string(enc) == "CONSTRAINT" ||
            @warn "bound_encoding is not CONSTRAINT — state boxes are variable bounds and " *
                  "are invisible to block discovery. Re-run with CONSTRAINT before " *
                  "trusting the charges." bound_encoding=string(enc)
    end
    cols_df.block = [c.family in Set(String.(spec.x)) ? "x" :
                     c.family == "P_g" ? "P" :
                     c.family == "Q_g" ? "Q" : "other"
                     for c in eachrow(cols_df)]

    n_F, n_x = length(F_rows), length(x_cols)
    square = n_F == n_x
    if !square
        counts = combine(groupby(filter(r -> r.row in Set(F_rows), rows_df), :family), nrow)
        @warn "Dynamic block is not square — S cannot be computed." n_F n_x
        println(counts)
    end
    length(h_rows) == length(hlo_rows) || @warn "Corridor upper/lower row counts differ" nu=length(h_rows) nl=length(hlo_rows)
    # An empty h block is the one failure mode that produces a plausible answer:
    # S becomes 0×n, every charge comes out exactly zero, and nothing complains.
    isempty(h_rows) && @warn "No δ corridor rows found — S will be empty and every " *
        "causer charge exactly zero. Either the run built no δ corridor (constrain_δ " *
        "= false, or both ineq_δ_COI_* toggles off), or it uses a corridor family " *
        "this spec does not list." families_h=spec.h

    meta_out = OrderedDict{String, Any}(
        "path" => spec.label,
        "families_F" => String.(spec.F),
        "families_x" => String.(spec.x),
        "families_h" => String.(spec.h),
        "families_h_lower" => String.(spec.h_lower),
        "n_F_rows" => n_F, "n_x_cols" => n_x, "square" => square,
        "families_h_omega" => String.(spec.h_ω),
        "families_h_omega_lower" => String.(spec.h_ω_lower),
        "n_h_rows" => length(h_rows), "n_P_cols" => length(P_cols),
        "n_h_omega_rows" => length(hω_rows),
        "n_x_block_rows" => length(block_rows),
        "x_block_families" => sort(unique(rows_df.family[block_rows])),
        "n_unmapped_rows" => count(==("UNMAPPED"), rows_df.family),
        "n_unmapped_cols" => count(==("UNMAPPED"), cols_df.family),
        "jacobian_nnz" => nnz(Jac),
    )

    CSV.write(joinpath(out_dir, "jacobian_rows.csv"), rows_df; delim = ';')
    CSV.write(joinpath(out_dir, "jacobian_cols.csv"), cols_df; delim = ';')
    let (Ic, Jc, Vc) = findnz(Jac)
        CSV.write(joinpath(out_dir, "jacobian_coo.csv"),
                  DataFrame(row = Ic, col = Jc, value = Vc); delim = ';')
    end

    if !square
        open(io -> TOML.print(io, _toml_safe(meta_out); sorted = true),
             joinpath(out_dir, "sensitivity_meta.toml"), "w")
        return out_dir
    end

    F_x = Jac[F_rows, x_cols]
    F_P = Jac[F_rows, P_cols]
    F_Q = isempty(Q_cols) ? spzeros(n_F, 0) : Jac[F_rows, Q_cols]
    h_x = Jac[h_rows, x_cols]
    h_P = Jac[h_rows, P_cols]
    hω_x = Jac[hω_rows, x_cols]
    hω_P = Jac[hω_rows, P_cols]

    # A corridor row whose nonzeros are not all inside x ∪ P is being silently
    # truncated: the slice then measures a DIFFERENT functional from the one the model
    # bounds, at plausible magnitude and with the identity still closing. This is
    # exactly what a missing `ΔωCOI` column does to the `:coi_box` speed corridor —
    # S_ω would come out as the sensitivity of the absolute speed.
    for (name, rws, hx, hp) in (("δ", h_rows, h_x, h_P), ("Δω", hω_rows, hω_x, hω_P))
        isempty(rws) && continue
        total, kept = nnz(Jac[rws, :]), nnz(hx) + nnz(hp)
        total == kept || @warn "$(name) corridor rows touch columns outside x ∪ P — the " *
            "sensitivity is being truncated to a different functional. Add the missing " *
            "variable family to spec.x (and its defining equality to spec.F)." total kept
    end

    # --- S = h_P + h_x(−F_x⁻¹F_P) -------------------------------------------------
    F_lu = lu(F_x)
    dxdP = -(F_lu \ Matrix(F_P))
    S_P = Matrix(h_P) + h_x * dxdP
    S_Q = isempty(Q_cols) ? zeros(length(h_rows), 0) : h_x * (-(F_lu \ Matrix(F_Q)))

    # The speed corridor: same construction, same factorisation, different functional.
    S_ω = Matrix(hω_P) + hω_x * dxdP

    # --- the discovered blocks: charges and adjoint contribution -------------------
    α = rows_df.dual[F_rows]
    γ̃ = rows_df.dual[h_rows] .- rows_df.dual[hlo_rows]
    γ̃ω = rows_df.dual[hω_rows] .- rows_df.dual[hωlo_rows]
    blocks_df, w_blocks, charge_blocks, bP_blocks =
        block_charges(Jac, rows_df, block_rows, x_cols, P_cols, dxdP)

    # Stationarity in x assembled from EVERY Jacobian row against its own dual, with no
    # reference to the F / h / block classification. That independence is the whole point:
    # `adj .+ w_blocks` below is built only from rows this file managed to classify, so a
    # row belonging to no channel is missing from it — but present here. The difference
    # between the two is exactly the unclassified charge (see the split further down).
    # The objective depends only on P_g, which is a P column and never an x column, so
    # ∇ₓf = 0 and the constraint Jacobian alone gives stationarity in x.
    r_true = Vector(Jac[:, x_cols]' * rows_df.dual)

    # --- validation 1: adjoint identity ------------------------------------------
    # Two versions, and the difference between them is the diagnostic.
    #
    #   `adjoint_residual_by_state`  — named channels only (both corridors), no blocks. A
    #                                  family above noise here names a state carrying an
    #                                  extra active block. Keep it: this is what pointed at
    #                                  `P_m` on the 179-bus run.
    #   `adjoint_residual_full_max`  — with the blocks folded in. THIS is the one that
    #                                  must be ~0; if it is not, a block was missed.
    adj = F_x' * α + h_x' * γ̃ + (isempty(hω_rows) ? zeros(n_x) : hω_x' * γ̃ω)
    adj_by_family = OrderedDict{String, Any}()
    for fam in unique(cols_df.family[x_cols])
        idx = findall(==(fam), cols_df.family[x_cols])
        adj_by_family[fam] = maximum(abs, adj[idx])
    end
    meta_out["adjoint_residual_max"] = maximum(abs, adj)
    meta_out["adjoint_residual_by_state"] = adj_by_family
    meta_out["adjoint_residual_full_max"] = maximum(abs, adj .+ w_blocks)

    # --- validation 2: the identity γ̃ᵀS_P = −dyn ----------------------------------
    # γ̃ᵀS_P = γ̃ᵀh_P + αᵀF_P. The h_P term is NOT always zero: it vanishes on
    # `bound_style_δ = :coi_box`, where the corridor is a plain box on the δ variable,
    # but on `:swing_propagated` the corridor is stated in substituted form and
    # carries P_g itself — so the corridor prices active power directly. Dropping the
    # term silently turns that whole channel into an "identity residual".
    charge = vec(γ̃' * S_P)                        # per P column
    gamma_h_P = vec(γ̃' * Matrix(h_P))
    alpha_F_P = vec(α' * Matrix(F_P))
    n_P = length(P_cols)
    charge_ω = isempty(hω_rows) ? zeros(n_P) : vec(γ̃ω' * S_ω)
    gamma_hω_P = isempty(hω_rows) ? zeros(n_P) : vec(γ̃ω' * Matrix(hω_P))
    meta_out["gamma_S_P"] = collect(charge)
    meta_out["gamma_h_P"] = collect(gamma_h_P)
    meta_out["alpha_F_P"] = collect(alpha_F_P)
    meta_out["gamma_binds"] = isempty(γ̃) ? 0.0 : maximum(abs, γ̃)
    meta_out["h_P_is_zero"] = iszero(nnz(h_P))
    meta_out["gamma_S_omega_P"] = collect(charge_ω)
    meta_out["gamma_omega_binds"] = isempty(γ̃ω) ? 0.0 : maximum(abs, γ̃ω)
    # Expected true: both speed styles are plain boxes on the Δω variable, and there is
    # no `:swing_propagated` analogue for speed. Recorded so a builder that ever pushes
    # P_g into the speed corridor shows up as a contribution, not as a residual.
    meta_out["h_omega_P_is_zero"] = iszero(nnz(hω_P))
    meta_out["block_charge_total"] = collect(charge_blocks)
    meta_out["block_b_P"] = collect(bP_blocks)

    # Corridor-only residual. Kept deliberately: it is not noise. Rearranging the joint
    # identity gives  residual = Σμᵀb_P − Σμᵀ S_b,P,  so with the usual b_P = 0 it is
    # exactly MINUS the summed block charge. That makes it an independent check on the
    # new machinery, and a one-number measure of how much of the dynamic price runs
    # outside the corridor. On a run with no active block it stays at solver noise.
    resid_corridor = charge .- (gamma_h_P .+ alpha_F_P)
    meta_out["identity_residual"] = maximum(abs, resid_corridor)
    meta_out["block_charge_vs_residual"] =
        maximum(abs, charge_blocks .+ charge_ω .+ resid_corridor .- bP_blocks .- gamma_hω_P)

    # The identity over every channel. Read `identity_residual_unexplained` below, NOT this
    # one: substituting S_P = h_P + h_x·dxdP, S_ω = hω_P + hω_x·dxdP and S_b,P = b_P + b_x·dxdP
    # — all built from the single `dxdP` above — and using dxdP = −F_x⁻¹F_P collapses the
    # line below to `(F_xᵀα + h_xᵀγ̃ + hω_xᵀγ̃ω + Σb_xᵀμ_b)ᵀ·dxdP` exactly. That leading factor
    # is the solver's own dual infeasibility (`adjoint_residual_full_max`, which reproduces
    # Ipopt's reported value digit for digit), so this quantity is ‖solver tolerance‖ scaled
    # by ‖dx/dP‖ and GROWS with the state count: ~1e-11 at n_x ≈ 9e3, ~1e-4 at n_x ≈ 2e4 over
    # a 5 s horizon. Bounding it by a fixed tolerance therefore fails on large correct runs.
    joint = charge .+ charge_ω .+ charge_blocks
    joint_vec = joint .- (gamma_h_P .+ gamma_hω_P .+ alpha_F_P .+ bP_blocks)
    meta_out["identity_residual_joint"] = maximum(abs, joint_vec)

    # The split that makes the residual a test again.
    #
    #   explained    — the same product formed with the CLASSIFICATION-FREE `r_true`. This is
    #                  the part attributable to the solver stopping at a finite dual
    #                  tolerance; it is expected to be nonzero and to grow with n_x.
    #   unexplained  — joint_vec − explained = ((adj + w_blocks) − r_true)ᵀ·dxdP, i.e. the
    #                  charge carried by rows that reached no channel. THIS is the number
    #                  that must be ~0. It is insensitive to solver tolerance, because both
    #                  terms inherit the same dual infeasibility and it cancels.
    explained = vec(r_true' * dxdP)
    meta_out["stationarity_residual_max"] = maximum(abs, r_true)
    meta_out["identity_residual_explained"] = collect(explained)
    meta_out["identity_residual_unexplained"] = maximum(abs, joint_vec .- explained)

    # --- validation 3: finite difference on the dynamic block ---------------------
    # Primal-side and dual-free, so it validates each S independently of every
    # multiplier. The helper only needs the corridor's own slices, so the speed
    # corridor reuses it unchanged.
    # Each check is skipped when its corridor has no rows. Running it on an empty
    # corridor divides by a zero scale and yields NaN, which `_toml_safe` then writes
    # as the STRING "NaN" (TOML has no nan literal) — and a reader formatting that with
    # a float spec dies. An absent key is honest; a stringified NaN is a trap.
    if !isempty(h_rows)
        meta_out["fd"] = _finite_difference_check(evaluator, rows_df, x_val, F_rows,
                                                  x_cols, P_cols, h_rows, h_x, S_P, F_lu)
    end
    if !isempty(hω_rows)
        meta_out["fd_omega"] = _finite_difference_check(evaluator, rows_df, x_val, F_rows,
                                                        x_cols, P_cols, hω_rows, hω_x, S_ω, F_lu)
    end

    # --- write ---------------------------------------------------------------------
    hkeys = rows_df[h_rows, :]
    S_long = DataFrame(t_index = Int[], gen_h = Int[], p_col = Int[],
                       gen_P = Int[], value = Float64[])
    P_gen = cols_df.key1[P_cols]
    for (r, hrow) in enumerate(eachrow(hkeys)), (c, pc) in enumerate(P_cols)
        push!(S_long, (hrow.key2, hrow.key1, pc, P_gen[c], S_P[r, c]))
    end
    CSV.write(joinpath(out_dir, "S_P.csv"), S_long; delim = ';')
    if !isempty(Q_cols)
        SQ_long = DataFrame(t_index = Int[], gen_h = Int[], gen_Q = Int[], value = Float64[])
        Q_gen = cols_df.key1[Q_cols]
        for (r, hrow) in enumerate(eachrow(hkeys)), (c, _) in enumerate(Q_cols)
            push!(SQ_long, (hrow.key2, hrow.key1, Q_gen[c], S_Q[r, c]))
        end
        CSV.write(joinpath(out_dir, "S_Q.csv"), SQ_long; delim = ';')
    end

    # The speed-corridor sensitivity, same long format as S_P.csv. Written only when the
    # corridor exists, so its absence in a run folder means "no speed corridor" rather
    # than "computed and empty".
    if !isempty(hω_rows)
        ωkeys = rows_df[hω_rows, :]
        Sω_long = DataFrame(t_index = Int[], gen_h = Int[], gen_P = Int[], value = Float64[])
        for (r, hrow) in enumerate(eachrow(ωkeys)), (c, _) in enumerate(P_cols)
            push!(Sω_long, (hrow.key2, hrow.key1, P_gen[c], S_ω[r, c]))
        end
        CSV.write(joinpath(out_dir, "S_omega_P.csv"), Sω_long; delim = ';')
    end

    # --- the paper table: one row per generator, one column per priced channel ------
    channels = DataFrame(gen_P = P_gen, charge_delta = charge, charge_omega = charge_ω)
    for r in eachrow(blocks_df)
        # Family names carry Greek and are not safe as bare column symbols; prefix and
        # keep them verbatim so the column maps back to the model family by inspection.
        channels[!, Symbol("charge_", r.family)] = r.charge
    end
    channels.charge_blocks_total = charge_blocks
    channels.charge_total = charge .+ charge_ω .+ charge_blocks
    channels.minus_dyn = gamma_h_P .+ gamma_hω_P .+ alpha_F_P .+ bP_blocks
    channels.residual = channels.charge_total .- channels.minus_dyn
    # Per machine, the same split as the two scalar meta keys: how much of this row's gap is
    # the solver's dual tolerance, and how much is a channel nobody claimed.
    channels.residual_explained = explained
    channels.residual_unexplained = channels.residual .- explained
    CSV.write(joinpath(out_dir, "charges_by_channel.csv"), channels; delim = ';')

    if !isempty(blocks_df)
        CSV.write(joinpath(out_dir, "x_blocks.csv"),
                  select(blocks_df, :family, :kind, :n_rows, :max_abs_dual); delim = ';')
    end

    open(io -> TOML.print(io, _toml_safe(meta_out); sorted = true),
         joinpath(out_dir, "sensitivity_meta.toml"), "w")

    println("Tier-2 sensitivity saved in: ", out_dir)
    println("  adjoint residual (channels only) = ", meta_out["adjoint_residual_max"])
    println("  adjoint residual (with blocks)   = ", meta_out["adjoint_residual_full_max"])
    println("  identity residual (δ only)       = ", meta_out["identity_residual"])
    println("  identity residual (joint)        = ", meta_out["identity_residual_joint"],
            "   (expected ≈ solver dual tol × ‖dx/dP‖, grows with n_x — not an error)")
    println("    of which explained by solver   = ", maximum(abs, explained))
    println("    UNEXPLAINED (must be ~0)       = ",
            meta_out["identity_residual_unexplained"])
    println("  stationarity residual max |r|    = ", meta_out["stationarity_residual_max"],
            "   (compare with Ipopt's \"Dual infeasibility\" in solver_log.txt)")
    haskey(meta_out, "fd") &&
        println("  FD rel err  δ                    = ", meta_out["fd"]["rel_err"])
    haskey(meta_out, "fd_omega") &&
        println("  FD rel err  Δω                   = ", meta_out["fd_omega"]["rel_err"])
    println("  corridor rows  δ=", length(h_rows), " (max|γ̃|=", meta_out["gamma_binds"],
            ")   Δω=", length(hω_rows), " (max|γ̃|=", meta_out["gamma_omega_binds"], ")")
    println("  max|charge|    δ=", isempty(charge) ? 0.0 : maximum(abs, charge),
            "   Δω=", isempty(charge_ω) ? 0.0 : maximum(abs, charge_ω))
    if !isempty(blocks_df)
        println("  extra blocks on x (", nrow(blocks_df), " families):")
        for r in eachrow(blocks_df)
            println("    ", rpad(r.family, 34), " rows=", lpad(r.n_rows, 7),
                    "  max|μ|=", r.max_abs_dual, "  max|charge|=", maximum(abs, r.charge))
        end
    end
    return out_dir
end

"""
Primal-side check, independent of every dual: perturb one `P_g` column by ε, solve
F(x, P+ε) = 0 by frozen-Jacobian Newton, and compare the resulting change in h
against the column of S_P.

`MOI.eval_constraint` on the evaluator that produced the Jacobian evaluates F, so no
part of the dynamics is re-implemented here.
"""
function _finite_difference_check(evaluator, rows_df, x_val, F_rows, x_cols, P_cols,
                                  h_rows, h_x, S_P, F_lu; ε::Float64 = 1e-6, which::Int = 1)
    n_all = length(rows_df.row)
    g_all = zeros(n_all)

    # Reference the residual against the evaluation at the SOLUTION rather than
    # against the stored right-hand sides: F(x*, P*) = 0 by construction, so
    # g(x) − g(x*) is the residual whatever convention MOI.Nonlinear uses for the
    # set constant. One less assumption to be wrong about.
    MOI.eval_constraint(evaluator, g_all, x_val)
    g_star = g_all[F_rows]

    xv = copy(x_val)
    xv[P_cols[which]] += ε
    resid = Float64[]
    for _ in 1:10
        MOI.eval_constraint(evaluator, g_all, xv)
        resid = g_all[F_rows] .- g_star
        maximum(abs, resid) < 1e-12 && break
        xv[x_cols] .-= F_lu \ resid
    end

    # h is affine in x on every path here (δ_g − δ_COI), so Δh = h_x·Δx exactly.
    Δh = h_x * (xv[x_cols] .- x_val[x_cols])
    S_fd = Δh ./ ε
    scale = maximum(abs, S_P[:, which])
    rel = scale > 0 ? maximum(abs, S_fd .- S_P[:, which]) / scale : NaN
    return OrderedDict{String, Any}(
        "column" => which, "epsilon" => ε,
        "newton_residual" => isempty(resid) ? NaN : maximum(abs, resid),
        "rel_err" => rel,
    )
end
