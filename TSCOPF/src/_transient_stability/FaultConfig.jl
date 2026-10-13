#=
================================================================================
 FaultConfig.jl  —  transient-stability disturbance specification
================================================================================
 Decouples fault configuration from steady-state settings.
 SC contingencies are read from `contingencies.csv`; GL (generator / load trip)
 and OB (open branch) parameters are set explicitly on `FaultConfig` inside
 `DynModelConfig`.

 `build_fault_details` returns the legacy `OrderedDict` consumed by Kron and
 FULL_BUS builders so existing save/export code keeps working.
================================================================================
=#

"""
    FaultType

Disturbance class for a transient-stability run:

- `SC` — three-phase bus short-circuit (from `contingencies.csv`)
- `GL` — generator or load disconnection (explicit ids on `FaultConfig`)
- `OB` — open branch: line/transformer trip with no short-circuit stage
"""
@enum FaultType SC GL OB

"""
    FaultConfig

Disturbance specification for a transient-stability run
(`RunConfig.transient.dyn_model.fault`).

`fault_type = SC` reads the faulted bus and branch to disconnect from
`contingencies.csv` via `contingency_id`. `fault_type = GL` trips a
generator (`gl_gen_ids`) or scales demand at specific buses
(`gl_load_bus_ids`, `gl_percent_power`); a run may specify one or the
other, not both. Load scaling: `p_d_new = p_d * (1 + α)` (same for `q_d`),
so `α = -1` is a full trip and `α = -0.5` sheds half the pre-fault demand.

`fault_type = OB` opens one or more in-service branches (`ob_branch_ids`,
DCIR row indices) with the same single-window timeline as GL (no fault-on
short-circuit stage). Opening a set that islands the network or isolates
the slack is rejected at validation.

`fault_type = SC` stamps `sc_fault_admittance_pu` on the faulted-bus diagonal
during the fault window. The default is the inductive solid short
`0 - im*FAULT_BUS_SHUNT` with `:add` (residual angle near the pre-fault phasor,
no active draw). The resistive solid short remains available as
`FAULT_BUS_SHUNT + 0im` with `:overwrite` (residual voltage near −90°, and a
small active-power sink at the faulted bus). These two fields are read only
for SC. In Ybus terms `+jB` is a capacitor, so a positive imaginary part is
warned about at validation.

See the user guide / parameter reference, section 5.
"""
Base.@kwdef struct FaultConfig
    fault_type::FaultType = SC
    contingency_id::Int = 2

    # SC — faulted-bus shunt. Default is the inductive solid short, added to the
    # network diagonal. The resistive option is `FAULT_BUS_SHUNT + 0im` with
    # `:overwrite`. Only read when fault_type == SC.
    # +jB is capacitive; an inductive short is a negative imaginary part.
    sc_fault_admittance_pu::ComplexF64 = -im * FAULT_BUS_SHUNT
    sc_fault_stamp::Symbol = :add

    # GL — generator disconnection (DGEN row indices)
    gl_gen_ids::Vector{Int} = Int[]

    # GL — load modification (parallel vectors; bus index ↔ scaling α)
    # Demand after disturbance: p_d_new = p_d * (1 + α), q_d_new = q_d * (1 + α).
    # Examples: α = -1 → full trip; α = -0.5 → 50% shed; α = 0 → unchanged; α = +1 → double.
    gl_load_bus_ids::Vector{Int} = Int[]
    gl_percent_power::Vector{Float64} = Float64[]

    # OB — open branch (DCIR row indices; lines and transformers)
    ob_branch_ids::Vector{Int} = Int[]
end

"""Fraction of pre-fault demand retained after GL load scaling with factor `α`."""
gl_load_demand_retained(α::Real) = 1 + α

"""
    apply_gl_load_scaling!(DBUS, bus_ids, percent_power)

Scale active/reactive demand at `bus_ids` for a GL load disturbance.

``p_{d,k}^{\\mathrm{new}} = p_{d,k}\\,(1 + \\alpha_k)`` (same for ``q_d``).
"""
function apply_gl_load_scaling!(
    DBUS::DataFrame,
    bus_ids::AbstractVector{<:Integer},
    percent_power::AbstractVector{<:Real},
)
    length(bus_ids) == length(percent_power) ||
        throw(ArgumentError("bus_ids and percent_power must have the same length."))
    for α in percent_power
        α < -1.0 && throw(ArgumentError(
            "gl_percent_power entries must be ≥ -1 (got α=$α); values below -1 imply negative demand."))
    end
    DBUS.p_d[bus_ids] = DBUS.p_d[bus_ids] .* (1 .+ percent_power)
    DBUS.q_d[bus_ids] = DBUS.q_d[bus_ids] .* (1 .+ percent_power)
    return DBUS
end

"""Deep-copy a `FaultConfig` (for `reconfigure` / `copy_dyn_model_config`)."""
function copy_fault_config(fault::FaultConfig)::FaultConfig
    return FaultConfig(;
        fault_type=fault.fault_type,
        contingency_id=fault.contingency_id,
        sc_fault_admittance_pu=fault.sc_fault_admittance_pu,
        sc_fault_stamp=fault.sc_fault_stamp,
        gl_gen_ids=copy(fault.gl_gen_ids),
        gl_load_bus_ids=copy(fault.gl_load_bus_ids),
        gl_percent_power=copy(fault.gl_percent_power),
        ob_branch_ids=copy(fault.ob_branch_ids),
    )
end

"""
Return `true` when every bus in `DBUS` remains reachable from the slack after
opening the circuits in `open_branch_ids` (and ignoring already-off branches).
"""
function _network_connected_from_slack(
    DBUS::DataFrame,
    DCIR::DataFrame,
    open_branch_ids::AbstractVector{<:Integer},
)::Bool
    nbus = nrow(DBUS)
    nbus == 0 && return true

    slack_rows = findall(==(3), DBUS.type)
    length(slack_rows) == 1 || return false
    slack_bus = Int(DBUS.bus[slack_rows[1]])

    open_set = Set{Int}(Int(b) for b in open_branch_ids)
    adj = [Int[] for _ in 1:nbus]
    for br in 1:nrow(DCIR)
        br in open_set && continue
        DCIR.l_status[br] != 1 && continue
        i = Int(DCIR.from_bus[br])
        k = Int(DCIR.to_bus[br])
        (1 <= i <= nbus && 1 <= k <= nbus) || continue
        push!(adj[i], k)
        push!(adj[k], i)
    end

    visited = falses(nbus)
    stack = Int[slack_bus]
    visited[slack_bus] = true
    while !isempty(stack)
        u = pop!(stack)
        for v in adj[u]
            if !visited[v]
                visited[v] = true
                push!(stack, v)
            end
        end
    end
    return all(visited)
end

"""
    _validate_sc_fault_admittance!(fault)

Fail fast on an unusable SC shunt. Called only for `fault_type == SC`.

A positive imaginary part is a capacitor (`+jB` in the Ybus), not an inductive
short, so it warns and is still accepted — the caller may have meant it.
"""
function _validate_sc_fault_admittance!(fault::FaultConfig)
    y = fault.sc_fault_admittance_pu
    (isfinite(real(y)) && isfinite(imag(y))) ||
        throw(ArgumentError("sc_fault_admittance_pu must be finite (got $y)."))
    abs(y) > 0 ||
        throw(ArgumentError("sc_fault_admittance_pu must be nonzero (got $y)."))
    real(y) >= 0 ||
        throw(ArgumentError(
            "sc_fault_admittance_pu must have real(y) ≥ 0 (got $y); " *
            "a negative conductance is not a fault resistance."))
    fault.sc_fault_stamp in (:overwrite, :add) ||
        throw(ArgumentError(
            "sc_fault_stamp must be :overwrite or :add (got $(fault.sc_fault_stamp))."))
    if imag(y) > 0
        @warn "sc_fault_admittance_pu has a positive imaginary part ($y). " *
              "In the Ybus, +jB is a capacitor, not an inductive short; " *
              "an inductive short is −jB."
    end
    return nothing
end

"""
    validate_fault_config!(fault, DGEN, DBUS, DCIR)

Fail fast on incoherent disturbance settings before building the TS model.
`DCIR` is required for `OB` (open-branch) checks and is accepted for all types.
"""
function validate_fault_config!(
    fault::FaultConfig,
    DGEN::DataFrame,
    DBUS::DataFrame,
    DCIR::DataFrame,
)
    has_ob = !isempty(fault.ob_branch_ids)

    if fault.fault_type == SC
        fault.contingency_id > 0 ||
            throw(ArgumentError("SC fault requires contingency_id > 0."))
        has_ob && throw(ArgumentError(
            "SC fault must leave ob_branch_ids empty (got $(fault.ob_branch_ids))."))
        _validate_sc_fault_admittance!(fault)
        return nothing
    end

    if fault.fault_type == GL
        has_ob && throw(ArgumentError(
            "GL fault must leave ob_branch_ids empty (got $(fault.ob_branch_ids))."))

        has_gen = !isempty(fault.gl_gen_ids)
        has_load = !isempty(fault.gl_load_bus_ids)
        if has_gen && has_load
            throw(ArgumentError(
                "GL fault cannot specify both gl_gen_ids and gl_load_bus_ids in one run."))
        end
        if !has_gen && !has_load
            throw(ArgumentError(
                "GL fault requires gl_gen_ids or gl_load_bus_ids to be non-empty."))
        end

        if has_gen
            for gen_id in fault.gl_gen_ids
                1 <= gen_id <= nrow(DGEN) ||
                    throw(ArgumentError("gl_gen_ids contains invalid generator id $gen_id."))
                DGEN.g_status[gen_id] == 1 ||
                    throw(ArgumentError("Generator $gen_id is not in service (g_status != 1)."))
                _is_slack_generator(gen_id, DGEN, DBUS) &&
                    throw(ArgumentError(
                        "Tripping the slack generator (gen $gen_id) is not supported."))
            end
            n_in_service = count(==(1), DGEN.g_status)
            length(fault.gl_gen_ids) >= n_in_service &&
                throw(ArgumentError(
                    "At least one generator must remain in service after GL gen trip."))
        end

        if has_load
            length(fault.gl_load_bus_ids) == length(fault.gl_percent_power) ||
                throw(ArgumentError(
                    "gl_load_bus_ids and gl_percent_power must have the same length."))
            for bus_id in fault.gl_load_bus_ids
                bus_id in DBUS.bus ||
                    throw(ArgumentError("gl_load_bus_ids contains invalid bus id $bus_id."))
            end
            for α in fault.gl_percent_power
                α < -1.0 && throw(ArgumentError(
                    "gl_percent_power entries must be ≥ -1 (got α=$α); values below -1 imply negative demand."))
            end
        end

        return nothing
    end

    fault.fault_type == OB || throw(ArgumentError("Unknown fault_type: $(fault.fault_type)."))

    (!isempty(fault.gl_gen_ids) || !isempty(fault.gl_load_bus_ids)) &&
        throw(ArgumentError(
            "OB fault must leave gl_gen_ids and gl_load_bus_ids empty."))
    has_ob || throw(ArgumentError(
        "OB fault requires ob_branch_ids to be non-empty."))

    ncir = nrow(DCIR)
    for br in fault.ob_branch_ids
        1 <= br <= ncir ||
            throw(ArgumentError("ob_branch_ids contains invalid circuit id $br."))
        DCIR.l_status[br] == 1 ||
            throw(ArgumentError(
                "Branch $br is not in service (l_status != 1); cannot open it."))
    end
    length(unique(fault.ob_branch_ids)) == length(fault.ob_branch_ids) ||
        throw(ArgumentError("ob_branch_ids must not contain duplicates."))

    _network_connected_from_slack(DBUS, DCIR, fault.ob_branch_ids) ||
        throw(ArgumentError(
            "OB open of branches $(fault.ob_branch_ids) would island the network " *
            "or isolate the slack bus."))

    return nothing
end

"""Return `true` when `gen_id` is the slack generator on `DBUS`."""
function _is_slack_generator(gen_id::Int, DGEN::DataFrame, DBUS::DataFrame)::Bool
    slack_rows = findall(==(3), DBUS.type)
    length(slack_rows) == 1 || return false
    return DGEN.bus[gen_id] == DBUS.bus[slack_rows[1]]
end

"""
    build_fault_details(fault, path_input_files; contingency_id) -> OrderedDict

Build the legacy `ts_fault_details` container used by TS model builders.
For SC, `contingency_id` keyword overrides `fault.contingency_id` when provided
(threaded from `RunConfig.contingency_id` for backward compatibility).
"""
function build_fault_details(
    fault::FaultConfig,
    path_input_files::String;
    contingency_id::Union{Nothing, Int}=nothing,
)::OrderedDict{Symbol, Any}
    ts_fault_details = OrderedDict{Symbol, Any}()
    ts_fault_details[:reduced_model] = true

    if fault.fault_type == SC
        ts_fault_details[:fault_type] = "SC"
        ts_fault_details[:sc] = OrderedDict{Symbol, Any}()
        ts_fault_details[:sc][:fault_location] = "bus"

        cid = something(contingency_id, fault.contingency_id)
        df_cont = CSV.read(joinpath(path_input_files, "contingencies.csv"), DataFrame; delim=';')
        row = df_cont[df_cont.contingency_nr .== cid, :]
        nrow(row) == 1 || throw(ArgumentError(
            "contingency_id=$cid not found in contingencies.csv."))

        ts_fault_details[:sc][:bus] = OrderedDict{Symbol, Any}(
            :bus_id => Int(row[1, :faulted_bus]),
            :disconnect_branch => true,
            :branch_id_2_disconnect => [Int(row[1, :circuit_id])],
            # Shunt stamped on the faulted-bus diagonal for the fault window only.
            :admittance_pu => fault.sc_fault_admittance_pu,
            :stamp => fault.sc_fault_stamp,
        )
        return ts_fault_details
    end

    if fault.fault_type == GL
        ts_fault_details[:fault_type] = "GL"
        ts_fault_details[:gl] = OrderedDict{Symbol, Any}()

        if !isempty(fault.gl_gen_ids)
            ts_fault_details[:gl][:element_2_disconnect] = "gen"
            ts_fault_details[:gl][:gen] = OrderedDict{Symbol, Any}(
                :gen_id => copy(fault.gl_gen_ids),
                :bus_id => Int[],  # filled by builder after reading DGEN copy
            )
        else
            ts_fault_details[:gl][:element_2_disconnect] = "load"
            ts_fault_details[:gl][:load] = OrderedDict{Symbol, Any}(
                :bus_id => copy(fault.gl_load_bus_ids),
                :percent_power => copy(fault.gl_percent_power),
            )
        end
        return ts_fault_details
    end

    fault.fault_type == OB || throw(ArgumentError("Unknown fault_type: $(fault.fault_type)."))
    ts_fault_details[:fault_type] = "OB"
    branch_ids = copy(fault.ob_branch_ids)
    ts_fault_details[:ob] = OrderedDict{Symbol, Any}(
        :branch_id => branch_ids,
    )
    # Annotate endpoints for logging when line_data.csv is available.
    line_path = joinpath(path_input_files, "line_data.csv")
    if isfile(line_path)
        df_line = CSV.read(line_path, DataFrame; delim=';')
        from_col = hasproperty(df_line, :fbus) ? :fbus :
                   hasproperty(df_line, :from_bus) ? :from_bus : nothing
        to_col = hasproperty(df_line, :tbus) ? :tbus :
                 hasproperty(df_line, :to_bus) ? :to_bus : nothing
        if from_col !== nothing && to_col !== nothing
            ts_fault_details[:ob][:from_bus] = [Int(df_line[br, from_col]) for br in branch_ids]
            ts_fault_details[:ob][:to_bus]   = [Int(df_line[br, to_col]) for br in branch_ids]
        end
    end
    return ts_fault_details
end
