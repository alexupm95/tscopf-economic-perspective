#=
================================================================================
 runtests_kkt_sensitivity_blocks.jl — offline check of the x-block ledger
================================================================================
 Replays the Tier-2 accounting against a run that has ALREADY been exported, using
 `Sensitivity/jacobian_{rows,cols,coo}.csv`. No solve, no live model.

 That is possible because those three files carry everything the identity needs: the
 exact Jacobian at the solution, every dual, and the family label of every row and
 column. So this exercises the real `x_block_rows` / `block_charges` on real data
 rather than on a constructed fixture.

 What it proves:

   1. adjoint  F_xᵀα + h_xᵀγ̃ + Σ_b b_xᵀμ_b = 0            (blocks complete the adjoint)
   2. identity γ̃ᵀS_P + Σ_b μ_bᵀS_b,P = γ̃ᵀh_P + αᵀF_P + Σ_b μ_bᵀb_P
   3. the corridor-only residual equals −Σ_b μ_bᵀS_b,P, i.e. the leak we measured
      before the blocks existed IS the block charge — an independent cross-check that
      does not reuse the same arithmetic path

 Run by hand (this file is deliberately not in test/runtests.jl, which ships):

     julia --project=TSCOPF kkt/runtests_kkt_sensitivity_blocks.jl <run_dir>

 With no argument it uses the run folder named by KKT_TEST_RUN.
================================================================================
=#

using CSV, DataFrames, SparseArrays, LinearAlgebra, DataStructures, TOML, Test

include(joinpath(@__DIR__, "kkt_sensitivity.jl"))


"""
    replay_run(run_dir) -> Bool

Replay one exported run. A FUNCTION rather than a top-level script so several runs
can be checked in ONE Julia process: the cold compile dominates the runtime, so a
process per directory both wastes it and pushes a multi-path sweep past the point
where a background task gets killed. Re-`include`ing this file instead would fail
anyway, on the redefinition of `struct DynStateSpec`.

Returns whether the run's testset passed, so a sweep can report one verdict.
"""
function replay_run(run_dir::AbstractString)::Bool
    # Tolerate a stray trailing CR: paths often arrive from a file written on Windows,
    # and the resulting directory name is invalid in a way that reads as "missing export".
    run_dir = rstrip(run_dir, ['\r', '\n', ' '])
    sens_dir = joinpath(run_dir, "Sensitivity")
    if !isdir(sens_dir)
        # Report and keep going: in a sweep, one unusable folder must not hide the
        # verdicts of every run after it.
        @warn "No Sensitivity/ folder — skipping" run_dir
        return false
    end

    @info "Replaying Tier-2 accounting" run_dir
    rows_df = CSV.read(joinpath(sens_dir, "jacobian_rows.csv"), DataFrame; delim = ';')
    cols_df = CSV.read(joinpath(sens_dir, "jacobian_cols.csv"), DataFrame; delim = ';')
    coo     = CSV.read(joinpath(sens_dir, "jacobian_coo.csv"), DataFrame; delim = ';')

    # Runs exported before `kind` was added to the row table still carry the information in
    # the family name, so old exports stay replayable.
    if !hasproperty(rows_df, :kind)
        rows_df.kind = [startswith(f, "eq_const") ? "eq" :
                        startswith(f, "ineq_const") ? "ineq" : "" for f in rows_df.family]
    end

    Jac = sparse(coo.row, coo.col, coo.value, nrow(rows_df), nrow(cols_df))

    # Pick the spec from the run's own record rather than assuming Kron, so a FULL_BUS
    # export replays through the same checks. `path` is written by the exporter from
    # `spec.label`, so the two cannot drift apart.
    meta_path = joinpath(sens_dir, "sensitivity_meta.toml")
    run_path = isfile(meta_path) ? get(TOML.parsefile(meta_path), "path", "") : ""
    spec = run_path == "FULL_BUS/CLASSICAL_2ND" ? FULLBUS_CLASSICAL_STATE :
           run_path == "FULL_BUS/DQ_4TH" ? DQ_FULLBUS_STATE :
           run_path in ("", "KRON_REDUCED/CLASSICAL_2ND") ? KRON_CLASSICAL_STATE :
           error("no spec for exported path $(run_path)")
    @info "replaying" path=(isempty(run_path) ? "(assumed Kron)" : run_path) spec=spec.label

    F_rows = _select(rows_df, spec.F, :row)
    h_rows, hlo_rows = corridor_row_pairs(rows_df, spec.h, spec.h_lower)
    hω_rows, hωlo_rows = corridor_row_pairs(rows_df, spec.h_ω, spec.h_ω_lower)
    x_cols = _select(cols_df, spec.x, :col)
    P_cols = _select(cols_df, [:P_g], :col)

    # Preconditions, NOT assertions: a bare `@test` outside a testset throws on failure,
    # which in a sweep aborts every later run. And an empty δ corridor is legitimate —
    # a speed-corridor-only run has none — so the requirement is that SOME corridor
    # exists, not that the angle one does.
    if length(F_rows) != length(x_cols)
        @warn "Dynamic block is not square — S is undefined, skipping" n_F=length(F_rows) n_x=length(x_cols)
        return false
    end
    if isempty(h_rows) && isempty(hω_rows)
        @warn "No corridor of either kind — nothing to price, skipping" run_dir
        return false
    end

    F_x = Jac[F_rows, x_cols];  F_P = Jac[F_rows, P_cols]
    h_x = Jac[h_rows, x_cols];  h_P = Jac[h_rows, P_cols]
    hω_x = Jac[hω_rows, x_cols]; hω_P = Jac[hω_rows, P_cols]

    F_lu = lu(F_x)
    dxdP = -(F_lu \ Matrix(F_P))
    S_P  = Matrix(h_P) + h_x * dxdP
    S_ω  = Matrix(hω_P) + hω_x * dxdP

    α = rows_df.dual[F_rows]
    γ̃ = rows_df.dual[h_rows] .- rows_df.dual[hlo_rows]
    γ̃ω = rows_df.dual[hω_rows] .- rows_df.dual[hωlo_rows]

    claimed    = Set{Int}(vcat(F_rows, h_rows, hlo_rows, hω_rows, hωlo_rows))
    block_rows = x_block_rows(Jac, claimed, x_cols)
    blocks_df, w_blocks, charge_blocks, bP_blocks =
        block_charges(Jac, rows_df, block_rows, x_cols, P_cols, dxdP)

    charge      = vec(γ̃' * S_P)
    charge_ω    = isempty(hω_rows) ? zeros(length(P_cols)) : vec(γ̃ω' * S_ω)
    gamma_h_P   = vec(γ̃' * Matrix(h_P))
    gamma_hω_P  = isempty(hω_rows) ? zeros(length(P_cols)) : vec(γ̃ω' * Matrix(hω_P))
    alpha_F_P   = vec(α' * Matrix(F_P))
    resid_corr  = charge .- (gamma_h_P .+ alpha_F_P)
    resid_joint = (charge .+ charge_ω .+ charge_blocks) .-
                  (gamma_h_P .+ gamma_hω_P .+ alpha_F_P .+ bP_blocks)

    adj_corr = F_x' * α + h_x' * γ̃
    adj_full = adj_corr .+ (isempty(hω_rows) ? zero(adj_corr) : hω_x' * γ̃ω) .+ w_blocks

    # Stationarity in x from every row and its own dual, independent of the F/h/block
    # classification — mirrors the exporter. `resid_joint` is algebraically adj_fullᵀ·dxdP,
    # so it carries the solver's dual tolerance amplified by ‖dx/dP‖ and grows with n_x.
    # Subtracting the classification-free version cancels that and leaves only the charge
    # of rows nobody claimed, which is what these assertions are actually about.
    r_true    = Vector(Jac[:, x_cols]' * rows_df.dual)
    explained = vec(r_true' * dxdP)

    # Values here run to ~45 EUR/MW, so 1e-6 absolute is ~2e-8 relative — tight enough to
    # catch a wrong sign or a dropped family, loose enough to survive the CSV round trip.
    TOL = 1e-6

    # Gross-breakage bound for the raw joint residual, which cannot use TOL: on a 23k-state
    # 5 s run it legitimately reaches 1e-4–1e-3 at Ipopt's default dual tolerance. Scaled to
    # the charges so it stays meaningful across case sizes.
    JOINT_SANITY = 1e-3 * max(maximum(abs, charge; init = 0.0),
                              maximum(abs, charge_ω; init = 0.0), 1.0)

    println("\n--- blocks discovered on x ---")
    if isempty(blocks_df)
        println("  (none — every constraint on x is F or the corridor)")
    else
        for r in eachrow(blocks_df)
            println("  ", rpad(r.family, 32), " kind=", rpad(r.kind, 5),
                    " rows=", lpad(r.n_rows, 7), "  max|μ|=", round(r.max_abs_dual, digits = 6),
                    "  max|charge|=", round(maximum(abs, r.charge), digits = 6))
        end
    end

    println("\n--- channels ---")
    println("  delta corridor rows    : ", length(h_rows), "  max|gamma|=", isempty(γ̃) ? 0.0 : maximum(abs, γ̃))
    println("  speed corridor rows    : ", length(hω_rows), "  max|gamma|=", isempty(γ̃ω) ? 0.0 : maximum(abs, γ̃ω))

    println("\n--- residuals ---")
    println("  adjoint, corridor only : ", maximum(abs, adj_corr))
    println("  adjoint, with blocks   : ", maximum(abs, adj_full))
    println("  identity, corridor only: ", maximum(abs, resid_corr))
    println("  identity, joint        : ", maximum(abs, resid_joint),
            "  (solver tol × ‖dx/dP‖, grows with n_x)")
    println("    explained by solver  : ", maximum(abs, explained))
    println("    UNEXPLAINED          : ", maximum(abs, resid_joint .- explained))
    println("  stationarity max |r|   : ", maximum(abs, r_true))

    # A failing `@testset` THROWS at its end, which in a sweep would abort every
    # remaining run and hide their results. Catch it, report false, keep going.
    ok = true
    try
    @testset "x-block ledger" begin
        # 1. the blocks complete the adjoint
        @test maximum(abs, adj_full) < TOL

        # 2. the joint identity closes once the solver's own dual infeasibility is removed.
        #    Bounding `resid_joint` itself would fail on any correct run large enough for
        #    ‖dx/dP‖ to amplify a 1e-9 dual residual past 1e-6 — which is every FULL_BUS
        #    DQ run at a 5 s horizon. The unexplained part is the classification test.
        @test maximum(abs, resid_joint .- explained) < TOL

        # 2b. …and the raw residual is still bounded, loosely, so a genuinely diverged run
        #     cannot pass by having its breakage land inside `explained`.
        @test maximum(abs, resid_joint) < JOINT_SANITY

        # 3. cross-check: the old corridor-only residual IS the block charge. Reached by a
        #    different route (S_P and γ̃) than charge_blocks (μᵀb_x·dx/dP), so agreement is
        #    evidence about the maths rather than about one code path. Same tolerance
        #    correction as 2 — this expression differs from `resid_joint` only by the speed
        #    channel, so it inherits the same amplification.
        @test maximum(abs, charge_blocks .+ charge_ω .+ resid_corr .- bP_blocks
                           .- gamma_hω_P .- explained) < TOL

        # 4. state boxes carry no P_g, so their direct term must vanish. If this ever fires,
        #    a block genuinely prices active power and the charge column is the place to look.
        @test maximum(abs, bP_blocks) < TOL

        # 5. The speed corridor must be a NAMED channel, never a discovered block. If any
        #    Delta-omega family shows up in the block list, it was not claimed by h_omega and
        #    its charge is being reported under the wrong heading.
        @test isempty(filter(r -> occursin("Δω", r.family) || occursin("Domega", r.family),
                             blocks_df))

        # 6. Same for the COI reference of the speed corridor. `eq_const_ΔωCOI_*` touching x
        #    while absent from spec.F is the exact bug this step fixes: it would leave S_omega
        #    measuring absolute speed instead of the COI-relative quantity the corridor bounds.
        @test isempty(filter(r -> occursin("ΔωCOI", r.family), blocks_df))

        # 7. Both corridors must be fully captured by the x and P slices. A nonzero in any
        #    other column means the functional is silently truncated.
        for (name, rows, hx, hp) in (("delta", h_rows, h_x, h_P), ("omega", hω_rows, hω_x, hω_P))
            isempty(rows) && continue
            @test nnz(Jac[rows, :]) == nnz(hx) + nnz(hp)
        end

        # 5. On the run that motivated all this, the leak was the P_m box and nothing else.
        #    Only asserted when that family is present, so the fixed run passes too.
        pm = filter(r -> occursin("P_m", r.family), blocks_df)
        if !isempty(pm)
            pm_charge = reduce(.+, pm.charge)
            @test maximum(abs, pm_charge .+ resid_corr) < TOL
            @test isapprox(maximum(abs, pm_charge), 44.37863245887598; atol = 1e-6)
            println("\n  P_m block reproduces the measured leak: max|charge| = ",
                    maximum(abs, pm_charge))
        end
    end   # @testset
    catch e
        e isa Test.TestSetException || rethrow()
        ok = false
    end
    return ok
end

runs = isempty(ARGS) ? [ENV["KKT_TEST_RUN"]] : ARGS
results = Pair{String, Bool}[]
for (i, r) in enumerate(runs)
    println("\n", "="^78)
    println("[", i, "/", length(runs), "] ", r)
    println("="^78)
    push!(results, string(r) => replay_run(r))
end

println("\n", "="^78)
for (r, ok) in results
    println(ok ? "  PASS  " : "  FAIL  ", r)
end
all_ok = all(last, results)
println(all_ok ? "ALL RUNS PASSED" : "SOME RUNS FAILED")
exit(all_ok ? 0 : 1)
