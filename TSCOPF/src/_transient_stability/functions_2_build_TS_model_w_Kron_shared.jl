#=
================================================================================
 functions_2_build_TS_model_w_Kron_shared.jl — shared TS helpers
================================================================================
=#

# ===================================================================================
# Phase 3 — classical init warm-start (E′, δ from solved or started OPF point)
# ===================================================================================

"""
    classical_E_delta_warmstart(v, θ, P_g, Q_g, Xd_tr) -> (E_mag, δ_angle)

Internal voltage behind transient reactance X′_d from the steady-state bus
voltage and generator dispatch.  Used by FULL_BUS ACOPF warm-start extraction.
"""
function classical_E_delta_warmstart(
    v_val::Float64,
    th_val::Float64,
    pg_val::Float64,
    qg_val::Float64,
    Xd_tr::Float64,
)::Tuple{Float64, Float64}
    V_c = v_val * exp(1im * th_val)
    I_c = (pg_val - 1im * qg_val) / conj(V_c)
    E_c = V_c + 1im * Xd_tr * I_c
    return abs(E_c), angle(E_c)
end
