# ===================================================================================
# Generator dynamic data — header-based parser (minimal 5-col or full 19-col CSV)
# ===================================================================================

const _GEN_DYN_OPTIONAL_COLS = (
    :Xq_tr, :Xd, :Xq, :Td, :Tq, :Ra, :T_exc, :K_exc, :Ta_exc, :Tb_exc,
    :R, :T1, :T2, :T3,
)

"""Normalise a CSV header string for column matching."""
_normalize_gen_dyn_header(name) = lowercase(strip(string(name)))

# CSV.jl types any column holding an empty cell as `Union{Missing, Float64}`. Map those
# cells to NaN so a blank parameter is indistinguishable from an absent column downstream:
# `validate_dyn_data!` then reports both spellings with the same message, naming the column
# and the flag that requires it. Plain `Float64.(col)` instead dies here with a bare
# `MethodError: no method matching Float64(::Missing)` that names neither.
_gen_dyn_float(v)::Float64 = ismissing(v) ? NaN : Float64(v)

"""
    Parse_Gen_Dynamic_DataFrame(df_raw::DataFrame) -> DataFrame

Parse `gen_dynamic_data.csv` (minimal) or `gen_dynamic_data_full.csv` (dq/AVR/GVNR).

Minimal format (`bus;Xd;H;D`): `Xd` is mapped to `Xd_tr` (legacy convention).
Legacy files may still include `Eg`; it is ignored (E is solved from init constraints).
Full format keeps separate `Xd_tr` (transient) and `Xd` (synchronous) columns.
Optional machine / control columns absent from the file are filled with `NaN`, and so
are blank cells in a column that is present — a partly-filled parameter and a missing
one are the same failure, reported together by `validate_dyn_data!`.
"""
function Parse_Gen_Dynamic_DataFrame(df_raw::DataFrame)::DataFrame
    raw_headers = [_normalize_gen_dyn_header(n) for n in names(df_raw)]
    has_xd_tr = "xd_tr" in raw_headers

    col_index = Dict{String, Int}()
    for (j, h) in enumerate(raw_headers)
        if h == "xd"
            if has_xd_tr
                col_index["xd_sync"] = j   # full file: synchronous reactance
            else
                col_index["xd_tr"] = j     # minimal file: Xd means Xd_tr
            end
        elseif h == "xd_tr"
            col_index["xd_tr"] = j
        elseif h == "bus"
            col_index["bus"] = j
        elseif h == "eg"
            col_index["eg"] = j
        elseif h == "xq_tr"
            col_index["xq_tr"] = j
        elseif h == "xq"
            col_index["xq"] = j
        elseif h == "td"
            col_index["td"] = j
        elseif h == "tq"
            col_index["tq"] = j
        elseif h == "h"
            col_index["h"] = j
        elseif h == "d"
            col_index["d"] = j
        elseif h == "ra"
            col_index["ra"] = j
        elseif h == "t_exc"
            col_index["t_exc"] = j
        elseif h == "k_exc"
            col_index["k_exc"] = j
        elseif h == "ta_exc"
            col_index["ta_exc"] = j
        elseif h == "tb_exc"
            col_index["tb_exc"] = j
        elseif h == "r"
            col_index["r"] = j
        elseif h == "t1"
            col_index["t1"] = j
        elseif h == "t2"
            col_index["t2"] = j
        elseif h == "t3"
            col_index["t3"] = j
        end
    end

    for req in ("bus", "xd_tr", "h", "d")
        haskey(col_index, req) ||
            throw(ArgumentError("gen_dynamic_data missing required column \"$req\" (headers: $(names(df_raw)))."))
    end

    n = nrow(df_raw)
    out = DataFrame(
        id     = collect(1:n),
        bus    = Vector{Int64}(undef, n),
        Xd_tr  = Vector{Float64}(undef, n),
        H      = Vector{Float64}(undef, n),
        D      = Vector{Float64}(undef, n),
    )
    for col in _GEN_DYN_OPTIONAL_COLS
        out[!, col] = fill(NaN, n)
    end

    for i in 1:n
        # `bus` stays a strict Int64 conversion — there is no NaN sentinel for an index,
        # so a blank bus is a malformed file with nothing sensible to defer to validation.
        out.bus[i]    = Int64(df_raw[i, col_index["bus"]])
        out.Xd_tr[i]  = _gen_dyn_float(df_raw[i, col_index["xd_tr"]])
        out.H[i]      = _gen_dyn_float(df_raw[i, col_index["h"]])
        out.D[i]      = _gen_dyn_float(df_raw[i, col_index["d"]])
    end

    optional_map = Dict{String, Symbol}(
        "xq_tr" => :Xq_tr, "xd_sync" => :Xd, "xq" => :Xq,
        "td" => :Td, "tq" => :Tq, "ra" => :Ra,
        "t_exc" => :T_exc, "k_exc" => :K_exc,
        "ta_exc" => :Ta_exc, "tb_exc" => :Tb_exc,
        "r" => :R, "t1" => :T1, "t2" => :T2, "t3" => :T3,
    )
    for (key, sym) in optional_map
        if haskey(col_index, key)
            out[!, sym] = _gen_dyn_float.(df_raw[!, col_index[key]])
        end
    end

    return out
end

"""Read and parse a generator-dynamic CSV from `folder_path` / `filename`."""
function Read_Gen_Dynamic_Data(folder_path::String; filename::String="gen_dynamic_data.csv")::DataFrame
    path = joinpath(folder_path, filename)
    isfile(path) || throw(ArgumentError("Generator dynamic data file not found: $path"))
    df_raw = CSV.read(path, DataFrame; delim=';')
    return Parse_Gen_Dynamic_DataFrame(df_raw)
end

"""`true` when every dq 4th-order machine column is present and finite."""
function has_full_machine_data(DGEN_DYN::DataFrame)::Bool
    dq_cols = (:Xq_tr, :Xd, :Xq, :Td, :Tq, :Ra)
    return all(c -> c in propertynames(DGEN_DYN) && all(isfinite, DGEN_DYN[!, c]), dq_cols)
end

# ===================================================================================
# GFM dynamic data — separate CSV from SG gen_dynamic_data (header-based)
# ===================================================================================

const _GFM_DYN_REQUIRED_COLS = (
    :Xl, :mq, :Kpv, :Kiv, :Emax, :Emin, :mp, :Pmax, :Pmin,
    :Qmax, :Qmin, :Kpqmax, :Kiqmax, :Kppmax, :Kipmax,
    :TPf, :TQf, :TVf, :Imax,
)

"""Measurement-filter time constants, in the order the P/Q/V filters consume them."""
const _GFM_FILTER_COLS = ("tpf", "tqf", "tvf")

"""
    Parse_GFM_Dynamic_DataFrame(df_raw; sys_base_MVA) -> DataFrame

Parse `gfm_dynamic_data.csv`. Required headers: `id`, `bus`, plus GFM controls
(`Xl`, `mq`, `mp`, `Kpv`, `Kiv`, `Emax`, `Emin`, `Pmax`, `Pmin`, `Qmax`, `Qmin`,
`Kppmax`, `Kipmax`, `Kpqmax`, `Kiqmax`, `Imax`) and the three measurement-filter
time constants `TPf`, `TQf`, `TVf` (REGFM_A1 eq. 4–6). `Kiqmax` is also accepted
as `kiqmax`. Optional `mach_base_MVA` (alias `InvBase`); defaults to `sys_base_MVA`.

There is no legacy path for a missing rating or gain: a GFM file that omits
`Qmax`, `Qmin`, `Kpqmax`, `Kiqmax`, `Kppmax` or `Kipmax` fails here, whether or
not `gfm_pq_overload` is on. Archived `RESULTS/*/Inputs/` copies and external
case folders have to gain those columns before they will parse.

A legacy single-constant `Tf` column is still accepted in place of the three and is
broadcast into all of them, so case folders and archived `RESULTS/*/Inputs/` copies
written before the split keep parsing unchanged. When a file carries both spellings the
three split columns win and `Tf` is ignored; carrying only some of the three and no `Tf`
throws.

Applies the reference machine→system base conversion on GFM rows after parse.
"""
function Parse_GFM_Dynamic_DataFrame(
    df_raw::DataFrame;
    sys_base_MVA::Float64,
)::DataFrame
    raw_headers = [_normalize_gen_dyn_header(n) for n in names(df_raw)]
    col_index = Dict{String, Int}()
    # `tpf`/`tqf`/`tvf` here are the lowercased CSV headers `TPf`/`TQf`/`TVf` — the
    # REGFM_A1 Table 1 symbols for the P/Q/V measurement lags. They are NOT the `_tpf`
    # post-fault window suffix the builders use; `col_index` is a private header
    # namespace and the two never meet. Keep the spec spelling so the CSV diffs against
    # PNNL-32278 without translation.
    for (j, h) in enumerate(raw_headers)
        if h in ("id", "bus", "xl", "mq", "kpv", "kiv", "emax", "emin", "mp",
                 "pmax", "pmin", "qmax", "qmin", "kpqmax", "kiqmax",
                 "kppmax", "kipmax", "tf", "tpf", "tqf", "tvf", "imax")
            col_index[h] = j
        elseif h in ("mach_base_mva", "invbase", "inv_base")
            col_index["mach_base_mva"] = j
        end
    end

    for req in ("id", "bus", "xl", "mq", "kpv", "kiv", "emax", "emin", "mp",
                "pmax", "pmin", "qmax", "qmin", "kpqmax", "kiqmax",
                "kppmax", "kipmax", "imax")
        haskey(col_index, req) ||
            throw(ArgumentError(
                "gfm_dynamic_data missing required column \"$req\" (headers: $(names(df_raw)))."))
    end

    # Either the three split constants or the legacy shared one. A file carrying only
    # some of the three and no `Tf` is a typo, not a legacy file, so it must fail loudly.
    split_filters = all(h -> haskey(col_index, h), _GFM_FILTER_COLS)
    legacy_filter = haskey(col_index, "tf")
    if !split_filters && !legacy_filter
        throw(ArgumentError(
            "gfm_dynamic_data needs the measurement-filter time constants: either all " *
            "of \"TPf\", \"TQf\", \"TVf\" (REGFM_A1 eq. 4–6) or the legacy single " *
            "\"Tf\" column, which is broadcast into all three " *
            "(headers: $(names(df_raw)))."))
    end

    n = nrow(df_raw)
    out = DataFrame(
        id   = Vector{Int64}(undef, n),
        bus  = Vector{Int64}(undef, n),
        Xl   = Vector{Float64}(undef, n),
        mq   = Vector{Float64}(undef, n),
        Kpv  = Vector{Float64}(undef, n),
        Kiv  = Vector{Float64}(undef, n),
        Emax = Vector{Float64}(undef, n),
        Emin = Vector{Float64}(undef, n),
        mp   = Vector{Float64}(undef, n),
        Pmax = Vector{Float64}(undef, n),
        Pmin = Vector{Float64}(undef, n),
        Qmax = Vector{Float64}(undef, n),
        Qmin = Vector{Float64}(undef, n),
        Kpqmax = Vector{Float64}(undef, n),
        Kiqmax = Vector{Float64}(undef, n),
        Kppmax = Vector{Float64}(undef, n),
        Kipmax = Vector{Float64}(undef, n),
        TPf  = Vector{Float64}(undef, n),
        TQf  = Vector{Float64}(undef, n),
        TVf  = Vector{Float64}(undef, n),
        Imax = Vector{Float64}(undef, n),
        mach_base_MVA = fill(sys_base_MVA, n),
    )

    for i in 1:n
        out.id[i]   = Int64(df_raw[i, col_index["id"]])
        out.bus[i]  = Int64(df_raw[i, col_index["bus"]])
        out.Xl[i]   = Float64(df_raw[i, col_index["xl"]])
        out.mq[i]   = Float64(df_raw[i, col_index["mq"]])
        out.Kpv[i]  = Float64(df_raw[i, col_index["kpv"]])
        out.Kiv[i]  = Float64(df_raw[i, col_index["kiv"]])
        out.Emax[i] = Float64(df_raw[i, col_index["emax"]])
        out.Emin[i] = Float64(df_raw[i, col_index["emin"]])
        out.mp[i]   = Float64(df_raw[i, col_index["mp"]])
        out.Pmax[i] = Float64(df_raw[i, col_index["pmax"]])
        out.Pmin[i] = Float64(df_raw[i, col_index["pmin"]])
        out.Qmax[i] = Float64(df_raw[i, col_index["qmax"]])
        out.Qmin[i] = Float64(df_raw[i, col_index["qmin"]])
        out.Kpqmax[i] = Float64(df_raw[i, col_index["kpqmax"]])
        out.Kiqmax[i] = Float64(df_raw[i, col_index["kiqmax"]])
        out.Kppmax[i] = Float64(df_raw[i, col_index["kppmax"]])
        out.Kipmax[i] = Float64(df_raw[i, col_index["kipmax"]])
        for col in ("pmax", "pmin", "qmax", "qmin", "kpqmax", "kiqmax", "kppmax", "kipmax")
            isfinite(Float64(df_raw[i, col_index[col]])) || throw(ArgumentError(
                "gfm_dynamic_data row $i column \"$col\" must be finite " *
                "(got $(df_raw[i, col_index[col]])). A missing rating is not a " *
                "silent \"no limit\"."))
        end
        # Split constants win when present; otherwise the legacy `Tf` fans out to all
        # three, which is exactly what the shared-constant model did before the split.
        if split_filters
            out.TPf[i] = Float64(df_raw[i, col_index["tpf"]])
            out.TQf[i] = Float64(df_raw[i, col_index["tqf"]])
            out.TVf[i] = Float64(df_raw[i, col_index["tvf"]])
        else
            Tf_legacy = Float64(df_raw[i, col_index["tf"]])
            out.TPf[i] = Tf_legacy
            out.TQf[i] = Tf_legacy
            out.TVf[i] = Tf_legacy
        end
        out.Imax[i] = Float64(df_raw[i, col_index["imax"]])
    end
    if haskey(col_index, "mach_base_mva")
        out.mach_base_MVA = Float64.(df_raw[!, col_index["mach_base_mva"]])
        out.mach_base_MVA[out.mach_base_MVA .== 0.0] .= sys_base_MVA
    end

    apply_gfm_base_conversion!(out, sys_base_MVA)
    return out
end

"""Scale GFM impedances/droops/limits from machine base to system base (in place)."""
function apply_gfm_base_conversion!(DGFM::DataFrame, sys_base_MVA::Float64)
    base_ratio = sys_base_MVA ./ DGFM.mach_base_MVA
    DGFM.Xl   .*= base_ratio
    DGFM.mq   .*= base_ratio
    DGFM.mp   .*= base_ratio
    # Power error → per-unit Δω, same normalisation as `mp`.
    DGFM.Kppmax .*= base_ratio
    DGFM.Kipmax .*= base_ratio
    # Power error → voltage command, same normalisation as `mq`.
    DGFM.Kpqmax .*= base_ratio
    DGFM.Kiqmax .*= base_ratio
    DGFM.Pmax ./= base_ratio
    DGFM.Pmin ./= base_ratio
    DGFM.Qmax ./= base_ratio
    DGFM.Qmin ./= base_ratio
    DGFM.Imax ./= base_ratio
    return DGFM
end

"""Read and parse `gfm_dynamic_data.csv` from `folder_path`."""
function Read_GFM_Dynamic_Data(
    folder_path::String;
    filename::String = "gfm_dynamic_data.csv",
    sys_base_MVA::Float64,
)::DataFrame
    path = joinpath(folder_path, filename)
    isfile(path) || throw(ArgumentError("GFM dynamic data file not found: $path"))
    df_raw = CSV.read(path, DataFrame; delim=';')
    return Parse_GFM_Dynamic_DataFrame(df_raw; sys_base_MVA=sys_base_MVA)
end

# ===================================================================================
# Load table — optional per-load split of a bus's demand, each load with its own ZIP
# ===================================================================================

const _LOAD_TABLE_COLS = (
    "id", "bus", "p_d", "q_d",
    "az_p", "ai_p", "ap_p", "az_q", "ai_q", "ap_q",
)

"""
    Parse_Load_DataFrame(df_raw) -> DataFrame

Parse `load_data.csv`: one row per individual load, several rows allowed on the same
bus. Required headers (case-insensitive): `id`, `bus`, `p_d`, `q_d`, plus the two ZIP
triples `aZ_p`, `aI_p`, `aP_p` (active demand) and `aZ_q`, `aI_q`, `aP_q` (reactive).

`p_d` / `q_d` are MW / MVAr on the same **unscaled** convention as `bus_data.csv`:
`load_system` scales this table and `DBUS` by `load_factor` together, after checking
that the rows of a bus sum to that bus's `Pd` / `Qd`. Coefficient order inside each
triple is (Z, I, P) — the same order as `DynModelConfig.zip_load_p`, and the same trap
(the reference implementation stores (P, I, Z)).
"""
function Parse_Load_DataFrame(df_raw::DataFrame)::DataFrame
    headers = [_normalize_gen_dyn_header(n) for n in names(df_raw)]
    col_index = Dict{String, Int}()
    for (j, h) in enumerate(headers)
        h in _LOAD_TABLE_COLS && (col_index[h] = j)
    end
    for req in _LOAD_TABLE_COLS
        haskey(col_index, req) || throw(ArgumentError(
            "load_data missing required column \"$req\" (headers: $(names(df_raw)))."))
    end

    n = nrow(df_raw)
    out = DataFrame(
        id   = Vector{Int64}(undef, n),
        bus  = Vector{Int64}(undef, n),
        p_d  = Vector{Float64}(undef, n),
        q_d  = Vector{Float64}(undef, n),
        aZ_p = Vector{Float64}(undef, n),
        aI_p = Vector{Float64}(undef, n),
        aP_p = Vector{Float64}(undef, n),
        aZ_q = Vector{Float64}(undef, n),
        aI_q = Vector{Float64}(undef, n),
        aP_q = Vector{Float64}(undef, n),
    )
    for i in 1:n
        out.id[i]   = Int64(df_raw[i, col_index["id"]])
        out.bus[i]  = Int64(df_raw[i, col_index["bus"]])
        out.p_d[i]  = Float64(df_raw[i, col_index["p_d"]])
        out.q_d[i]  = Float64(df_raw[i, col_index["q_d"]])
        out.aZ_p[i] = Float64(df_raw[i, col_index["az_p"]])
        out.aI_p[i] = Float64(df_raw[i, col_index["ai_p"]])
        out.aP_p[i] = Float64(df_raw[i, col_index["ap_p"]])
        out.aZ_q[i] = Float64(df_raw[i, col_index["az_q"]])
        out.aI_q[i] = Float64(df_raw[i, col_index["ai_q"]])
        out.aP_q[i] = Float64(df_raw[i, col_index["ap_q"]])
    end
    return out
end

"""
Read and parse `load_data.csv` from `folder_path`. Returns `nothing` when the file is
absent — the table is optional, and without it demand stays a per-bus aggregate carrying
the system-wide `zip_load_p` / `zip_load_q` splits.
"""
function Read_Load_Data(
    folder_path::String;
    filename::String = "load_data.csv",
)::Union{Nothing, DataFrame}
    path = joinpath(folder_path, filename)
    isfile(path) || return nothing
    df_raw = CSV.read(path, DataFrame; delim=';')
    return Parse_Load_DataFrame(df_raw)
end

# Function used to Read the Input data from the CSV files and store them into Structs
function Read_Input_Data(folder_path::String, TSA::Bool; gen_dynamic_filename::String="gen_dynamic_data.csv")

    # ====================================================================================================
    # If some modification is done in the name of the variables in the CSV files, it must be modified here
    # ==================================================================================================== 

    # Function to read buses data and store in DBUS_Struct
    function read_bus_data()
        df = CSV.read(joinpath(folder_path, "bus_data.csv"), DataFrame; delim=';')  # Read CSV file (absolute path)

        df_bus = deepcopy(df)
        rename!(df_bus, Dict(
            names(df_bus)[1] => :bus,
            names(df_bus)[2] => :type,
            names(df_bus)[3] => :p_d,
            names(df_bus)[4] => :q_d,
            names(df_bus)[5] => :g_sh,
            names(df_bus)[6] => :b_sh,
            names(df_bus)[7] => :area,
            names(df_bus)[8] => :v_spe,
            names(df_bus)[9] => :v_a,
            names(df_bus)[10] => :base_kV,
            names(df_bus)[11] => :zone,
            names(df_bus)[12] => :v_max,
            names(df_bus)[13] => :v_min,
        ))

        df_bus.bus = Int64.(df_bus.bus)
        df_bus.type = Int64.(df_bus.type)
        df_bus.p_d = Float64.(df_bus.p_d)
        df_bus.q_d = Float64.(df_bus.q_d)
        df_bus.g_sh = Float64.(df_bus.g_sh)
        df_bus.b_sh = Float64.(df_bus.b_sh)
        df_bus.area = Int64.(df_bus.area)
        df_bus.v_spe = Float64.(df_bus.v_spe)
        df_bus.v_a = Float64.(df_bus.v_a)
        df_bus.base_kV = Float64.(df_bus.base_kV)
        df_bus.zone = Int64.(df_bus.zone)
        df_bus.v_max = Float64.(df_bus.v_max)
        df_bus.v_min = Float64.(df_bus.v_min)

        return df_bus
    end

    # Function to read generators data and store in DGEN_Struct
    function read_gen_data()
        df = CSV.read(joinpath(folder_path, "generators_data.csv"), DataFrame; delim=';') # Read CSV file (absolute path)
        num_gen = length(df.bus)
        id = collect(1:num_gen)

        df_gen = hcat(DataFrame(id = id), df; makeunique=true)

        rename!(df_gen, Dict(
            names(df_gen)[2] => :bus,
            names(df_gen)[3] => :pg_spe,
            names(df_gen)[4] => :qg_spe,
            names(df_gen)[5] => :qg_max,
            names(df_gen)[6] => :qg_min,
            names(df_gen)[7] => :vg_spe,
            names(df_gen)[8] => :base_MVA,
            names(df_gen)[9] => :g_status,
            names(df_gen)[10] => :pg_max,
            names(df_gen)[11] => :pg_min,
            names(df_gen)[12] => :g_cost_2,
            names(df_gen)[13] => :g_cost_1,
            names(df_gen)[14] => :g_cost_0,
        ))

        df_gen.id = Int64.(df_gen.id)
        df_gen.bus = Int64.(df_gen.bus)
        df_gen.pg_spe = Float64.(df_gen.pg_spe)
        df_gen.qg_spe = Float64.(df_gen.qg_spe)
        df_gen.qg_max = Float64.(df_gen.qg_max)
        df_gen.qg_min = Float64.(df_gen.qg_min)
        df_gen.vg_spe = Float64.(df_gen.vg_spe)
        df_gen.base_MVA = Float64.(df_gen.base_MVA)
        df_gen.g_status = Int64.(df_gen.g_status)
        df_gen.pg_max = Float64.(df_gen.pg_max)
        df_gen.pg_min = Float64.(df_gen.pg_min)
        df_gen.g_cost_2 = Float64.(df_gen.g_cost_2)
        df_gen.g_cost_1 = Float64.(df_gen.g_cost_1)
        df_gen.g_cost_0 = Float64.(df_gen.g_cost_0)
        
        return df_gen
    end

    # Function to read generators dynamic data and store in DGEN_DYNAMIC_Struct
    function read_gen_dynamic_data()
        return Read_Gen_Dynamic_Data(folder_path; filename=gen_dynamic_filename)
    end

    # Function to read circuits data and store in DCIR_Struct
    function read_circuit_data()
        df = CSV.read(joinpath(folder_path, "line_data.csv"), DataFrame; delim=';')  # (absolute path)
        num_circ = length(df.fbus)
        id = collect(1:num_circ)

        df_cir = hcat(DataFrame(circ = id), df; makeunique=true)

        rename!(df_cir, Dict(
            names(df_cir)[1] => :id,
            names(df_cir)[2] => :from_bus,
            names(df_cir)[3] => :to_bus,
            names(df_cir)[4] => :l_res,
            names(df_cir)[5] => :l_reac,
            names(df_cir)[6] => :l_sh_susp,
            names(df_cir)[7] => :l_cap_1,
            names(df_cir)[8] => :l_cap_2,
            names(df_cir)[9] => :l_cap_3,
            names(df_cir)[10] => :t_tap,
            names(df_cir)[11] => :t_shift,
            names(df_cir)[12] => :l_status,
            names(df_cir)[13] => :ang_min,
            names(df_cir)[14] => :ang_max,
        ))

        df_cir.id = Int64.(df_cir.id)
        df_cir.from_bus = Int64.(df_cir.from_bus)
        df_cir.to_bus = Int64.(df_cir.to_bus)
        df_cir.l_res = Float64.(df_cir.l_res)
        df_cir.l_reac = Float64.(df_cir.l_reac)
        df_cir.l_sh_susp = Float64.(df_cir.l_sh_susp)
        df_cir.l_cap_1 = Float64.(df_cir.l_cap_1)
        df_cir.l_cap_2 = Float64.(df_cir.l_cap_2)
        df_cir.l_cap_3 = Float64.(df_cir.l_cap_3)
        df_cir.t_tap = Float64.(df_cir.t_tap)
        df_cir.t_shift = Float64.(df_cir.t_shift)
        df_cir.l_status = Int64.(df_cir.l_status)
        df_cir.ang_min = Float64.(df_cir.ang_min)
        df_cir.ang_max = Float64.(df_cir.ang_max)

        replace!(df_cir.t_tap, 0.0 => 1.0) # t_tap = 0 means ignore transformers in matpower-like files. By setting it to 1 it will has no impact on the model (leaving it to 0 raises zero dividion error in the following otherwise)

        return df_cir
    end

    # Input files are read via absolute paths (joinpath(folder_path, ...)) inside
    # the reader closures above, so no cd into folder_path is needed.

    DBUS     = read_bus_data()          # Generate the Struct with Buses data
    DGEN     = read_gen_data()          # Generate the Struct with Generators data
    if (TSA) DGEN_DYN = read_gen_dynamic_data() end # Generate the Struct with Generators Dynamic data
    DCIR     = read_circuit_data()      # Generate the Struct with Circuits data

    # For the code to work properly, the bus indices must be set in ascending order from 1 to nBUS
    bus_mapping, reverse_bus_mapping = Mapping_Buses_Labels(DBUS) # Map the buses labels from old to new nomeclature
    DBUS.bus = [bus_mapping[b] for b in DBUS.bus]                              # Rename the buses labels from 1 to nBUS

    # Map the buses labels to be in ascending order from 1 to nBUS
    DGEN.bus      = [bus_mapping[b] for b in DGEN.bus]
    if (TSA) DGEN_DYN.bus  = [bus_mapping[b] for b in DGEN_DYN.bus] end
    DCIR.from_bus = [bus_mapping[b] for b in DCIR.from_bus]
    DCIR.to_bus   = [bus_mapping[b] for b in DCIR.to_bus]

    # Invariant the whole pipeline relies on (review R3): after remapping, bus labels
    # are EXACTLY 1:nBUS in order, so a DataFrame row position equals its bus label.
    # Many builders index DBUS columns positionally by bus id (e.g. DBUS.p_d[bus]); if a
    # future change broke the remap, those reads would silently hit the wrong row. Assert
    # it here so the failure is loud and immediate instead of producing wrong results.
    @assert DBUS.bus == collect(1:nrow(DBUS)) "bus-remap invariant violated: DBUS.bus must equal 1:nBUS after Mapping_Buses_Labels (got $(DBUS.bus))."

    if (TSA)
        return DBUS, DGEN, DGEN_DYN, DCIR, bus_mapping, reverse_bus_mapping # Return the data
    else
        return DBUS, DGEN, DCIR, bus_mapping, reverse_bus_mapping # Return the data
    end
end

# Function used to map the from old to new nomeclature
function Mapping_Buses_Labels(DBUS::DataFrame)

    # Given bus numbers
    original_buses = DBUS.bus

    # Create a dictionary that maps original bus labels to new indices
    bus_mapping = OrderedDict(original_buses[i] => i for i in eachindex(original_buses))

    # Reverse mapping (for converting back later)
    reverse_bus_mapping = OrderedDict(i => original_buses[i] for i in eachindex(original_buses))

    return bus_mapping, reverse_bus_mapping
end

# Function that can change the buses labels according to the new nomenclature
function Change_Buses_Labels(DBUS::DataFrame, DGEN::DataFrame, DGEN_DYN::DataFrame, DCIR::DataFrame, bus_mapping::OrderedDict)
    
    # Convert using the reverse mapping
    DBUS.bus      = [bus_mapping[b] for b in DBUS.bus]
    DGEN.bus      = [bus_mapping[b] for b in DGEN.bus]
    DGEN_DYN.bus  = [bus_mapping[b] for b in DGEN_DYN.bus]
    DCIR.from_bus = [bus_mapping[b] for b in DCIR.from_bus]
    DCIR.to_bus   = [bus_mapping[b] for b in DCIR.to_bus]

    return DBUS, DGEN, DGEN_DYN, DCIR
end

# Function that can return the buses labels according to the original nomenclature
function Reverse_Buses_Labels(DBUS::DataFrame, DGEN::DataFrame, DGEN_DYN::DataFrame, DCIR::DataFrame, reverse_bus_mapping::OrderedDict)
        
    # Convert using the reverse mapping
    DBUS.bus      = [reverse_bus_mapping[b] for b in DBUS.bus]
    DGEN.bus      = [reverse_bus_mapping[b] for b in DGEN.bus]
    DGEN_DYN.bus  = [reverse_bus_mapping[b] for b in DGEN_DYN.bus]
    DCIR.from_bus = [reverse_bus_mapping[b] for b in DCIR.from_bus]
    DCIR.to_bus   = [reverse_bus_mapping[b] for b in DCIR.to_bus]

    return DBUS, DGEN, DGEN_DYN, DCIR
end


# Function used to map the generators, circuits connected and adjacent buses for each bus
function Organize_Bus_Gen_Circ(DBUS::DataFrame, DGEN::DataFrame, DCIR::DataFrame)

    # Initialize the dictionary with empty vectors for each bus
    bus_gen_circ_dict = OrderedDict(b => Dict(:gen_ids => Int[], :gen_status => Bool[], :circ => Int[], :adj_buses => Int[], :vg_spe => 1.0, :pg_max => 0.0, :pg_min => 0.0, :qg_max => 0.0, :qg_min => 0.0, :pd_tot => 0.0, :qd_tot => 0.0, :gsh_tot => 0.0, :bsh_tot => 0.0) for b in DBUS.bus)

    bus_gen_circ_dict_ON = OrderedDict(b => Dict(:gen_ids => Int[], :gen_status => Bool[], :circ => Int[], :adj_buses => Int[], :vg_spe => 1.0, :pg_max => 0.0, :pg_min => 0.0, :qg_max => 0.0, :qg_min => 0.0, :pd_tot => 0.0, :qd_tot => 0.0, :gsh_tot => 0.0, :bsh_tot => 0.0) for b in DBUS.bus)

    for bus_id in DBUS.bus
        bus_gen_circ_dict[bus_id][:pd_tot] = DBUS.p_d[bus_id]
        bus_gen_circ_dict[bus_id][:qd_tot] = DBUS.q_d[bus_id]
        bus_gen_circ_dict[bus_id][:gsh_tot] = DBUS.g_sh[bus_id]
        bus_gen_circ_dict[bus_id][:bsh_tot] = DBUS.b_sh[bus_id]

        bus_gen_circ_dict_ON[bus_id][:pd_tot] = DBUS.p_d[bus_id]
        bus_gen_circ_dict_ON[bus_id][:qd_tot] = DBUS.q_d[bus_id]
        bus_gen_circ_dict_ON[bus_id][:gsh_tot] = DBUS.g_sh[bus_id]
        bus_gen_circ_dict_ON[bus_id][:bsh_tot] = DBUS.b_sh[bus_id]
    end

    # Loop through the generators and fill the dictionary
    # Associate the number of the bus with a vector containing
    # the ids of the generators
    for (gen_idx, bus_id) in enumerate(DGEN.bus)
        push!(bus_gen_circ_dict[bus_id][:gen_ids], gen_idx)
        push!(bus_gen_circ_dict[bus_id][:gen_status], DGEN.g_status[gen_idx])
        bus_gen_circ_dict[bus_id][:pg_max] += DGEN.pg_max[gen_idx]
        bus_gen_circ_dict[bus_id][:pg_min] += DGEN.pg_min[gen_idx]
        bus_gen_circ_dict[bus_id][:qg_max] += DGEN.qg_max[gen_idx]
        bus_gen_circ_dict[bus_id][:qg_min] += DGEN.qg_min[gen_idx]

        # All generators connected to the same bus must share the same voltage
        # setpoint. Check the unique set across every gen seen so far at this bus
        # (gen_ids was just appended above), not a single scalar (which `unique`
        # would always reduce to length 1, making the check unreachable).
        vg_at_bus = unique(DGEN.vg_spe[bus_gen_circ_dict[bus_id][:gen_ids]])
        if length(vg_at_bus) > 1
            throw(ArgumentError("The voltages specified for the generators connected to bus $bus_id must be equal."))
        else
            bus_gen_circ_dict[bus_id][:vg_spe] = vg_at_bus[1]
        end

        if DGEN.g_status[gen_idx] == 1
            push!(bus_gen_circ_dict_ON[bus_id][:gen_ids], gen_idx)
            push!(bus_gen_circ_dict_ON[bus_id][:gen_status], DGEN.g_status[gen_idx])
            bus_gen_circ_dict_ON[bus_id][:pg_max] += DGEN.pg_max[gen_idx] 
            bus_gen_circ_dict_ON[bus_id][:pg_min] += DGEN.pg_min[gen_idx]  
            bus_gen_circ_dict_ON[bus_id][:qg_max] += DGEN.qg_max[gen_idx] 
            bus_gen_circ_dict_ON[bus_id][:qg_min] += DGEN.qg_min[gen_idx]  

            # Same equal-setpoint check, restricted to the ON generators at the bus.
            vg_at_bus_ON = unique(DGEN.vg_spe[bus_gen_circ_dict_ON[bus_id][:gen_ids]])
            if length(vg_at_bus_ON) > 1
                throw(ArgumentError("The voltages specified for the generators connected to bus $bus_id must be equal."))
            else
                bus_gen_circ_dict_ON[bus_id][:vg_spe] = vg_at_bus_ON[1]
            end
        end

    end

    # Loop through circuits and populate ids of circuits connectec
    # and adjacent buses
    for (i, from) in enumerate(DCIR.from_bus)
        to = DCIR.to_bus[i]
        
        # Add the number of the circuit
        push!(bus_gen_circ_dict[from][:circ], i)
        push!(bus_gen_circ_dict[to][:circ], i)

        # Add each bus as adjacent to the other
        push!(bus_gen_circ_dict[from][:adj_buses], to)
        push!(bus_gen_circ_dict[to][:adj_buses], from)

        if DCIR.l_status[i] == 1
            # Add the number of the circuit
            push!(bus_gen_circ_dict_ON[from][:circ], i)
            push!(bus_gen_circ_dict_ON[to][:circ], i)

            # Add each bus as adjacent to the other
            push!(bus_gen_circ_dict_ON[from][:adj_buses], to)
            push!(bus_gen_circ_dict_ON[to][:adj_buses], from)
        end

    end
    
    return bus_gen_circ_dict, bus_gen_circ_dict_ON
end