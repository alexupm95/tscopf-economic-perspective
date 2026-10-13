#=
================================================================================
 _manage_inputs/functions_2_parse_matpower.jl
   Standalone MATPOWER (.m) case parser → in-house DataFrames + CSV files
================================================================================
 Purpose
 -------
 Read a typical MATPOWER case file (e.g. INPUT_FILES/PowerModels/case9.m) and
 produce the network data in this codebase's standard form:

   * the in-memory `DBUS` / `DGEN` / `DCIR` DataFrames (exactly the schema that
     `Read_Input_Data` returns: id columns, renamed/typed fields, t_tap 0→1 fix,
     buses remapped to 1..nBUS); and
   * the three positional ';'-delimited CSVs (`bus_data.csv`,
     `generators_data.csv`, `line_data.csv`) written into a timestamped
     `RESULTS/Results - <ts>/Inputs/` folder.

 This parser is deliberately self-contained: it does NOT depend on PowerModels
 (which is a test-only dependency here). It hand-reads the `mpc.*` matrix blocks.

 Column mapping (the in-house CSV order already mirrors MATPOWER):
   mpc.bus                         -> bus_data.csv   (13 cols, verbatim)
   mpc.gen[:, 1:10] + gencost c2,c1,c0 -> generators_data.csv (13 cols)
   mpc.branch[:, 1:13]             -> line_data.csv  (13 cols)

 Scope note: MATPOWER carries no dynamic-machine data, so `DGEN_DYN` is NOT
 produced. Building a transient-ready case still needs a hand-authored
 `gen_dynamic_data.csv` (and `contingencies.csv` for SC faults).

 Entry point: `Import_Matpower_Case(path_m)`.
================================================================================
=#

# ------------------------------------------------------------------ block names
const _MPC_BLOCKS = (:bus, :gen, :branch, :gencost)

# ==============================================================================
#  1. Lexer
# ==============================================================================
"""
    _strip_matpower_comments(text) -> String

Drop MATPOWER `%` comments (to end of line) and normalise CR/LF so the block
extractor sees clean, `\\n`-separated text.
"""
function _strip_matpower_comments(text::AbstractString)::String
    text = replace(text, "\r\n" => "\n", "\r" => "\n")
    cleaned = IOBuffer()
    for line in split(text, '\n')
        c = findfirst('%', line)                 # first comment marker on the line
        write(cleaned, c === nothing ? line : line[1:c-1])
        write(cleaned, '\n')
    end
    return String(take!(cleaned))
end

"""
    _parse_matpower_matrix(text, name) -> Union{Matrix{Float64}, Nothing}

Extract the numeric matrix assigned to `mpc.<name> = [ ... ];`. Returns a
`Matrix{Float64}` (rows × columns) or `nothing` when the block is absent. Rows
are terminated by `;`; cells are whitespace-separated. Rows are allowed to carry
a different number of trailing columns — they are padded/handled by the callers.
"""
function _parse_matpower_matrix(text::AbstractString, name::Symbol)
    # Locate "mpc.<name>" then the opening '[' and matching closing ']'.
    key = "mpc." * String(name)
    kpos = findfirst(key, text)
    kpos === nothing && return nothing
    open_br = findnext('[', text, last(kpos))
    open_br === nothing && return nothing
    close_br = findnext(']', text, open_br)
    close_br === nothing &&
        throw(ArgumentError("MATPOWER block mpc.$name has no closing ']'."))

    body = text[open_br+1:close_br-1]

    rows = Vector{Vector{Float64}}()
    for raw_row in split(body, ';')
        cells = split(strip(raw_row))            # split on any whitespace run
        isempty(cells) && continue               # blank line / trailing terminator
        push!(rows, [parse(Float64, c) for c in cells])
    end
    isempty(rows) && return nothing

    ncol = maximum(length, rows)
    M = Matrix{Float64}(undef, length(rows), ncol)
    fill!(M, NaN)
    for (i, r) in enumerate(rows)
        M[i, 1:length(r)] .= r                   # ragged rows left as NaN past their end
    end
    return M
end

"""
    Parse_Matpower_File(path_m) -> (; baseMVA, bus, gen, branch, gencost)

Pure text parse of a MATPOWER `.m` case. Each `mpc.*` data block is returned as a
`Matrix{Float64}` (or `nothing` if absent — only `gencost` is allowed to be
absent). `baseMVA` is read from `mpc.baseMVA = <value>;`.
"""
function Parse_Matpower_File(path_m::String)
    isfile(path_m) || throw(ArgumentError("MATPOWER case file not found: $path_m"))
    text = _strip_matpower_comments(read(path_m, String))

    # baseMVA scalar
    m = match(r"mpc\.baseMVA\s*=\s*([0-9.eE+\-]+)", text)
    m === nothing && throw(ArgumentError("Could not find `mpc.baseMVA` in $path_m."))
    baseMVA = parse(Float64, m.captures[1])

    blocks = Dict{Symbol, Any}()
    for name in _MPC_BLOCKS
        blocks[name] = _parse_matpower_matrix(text, name)
    end

    for req in (:bus, :gen, :branch)
        blocks[req] === nothing &&
            throw(ArgumentError("MATPOWER file $path_m is missing the mpc.$req block."))
    end

    return (; baseMVA,
              bus     = blocks[:bus]::Matrix{Float64},
              gen     = blocks[:gen]::Matrix{Float64},
              branch  = blocks[:branch]::Matrix{Float64},
              gencost = blocks[:gencost])
end

# ==============================================================================
#  2. Schema mapping
# ==============================================================================
"""
    Matpower_Gencost_To_Quadratic(gencost, ng) -> (c2, c1, c0)

Map the MATPOWER `mpc.gencost` matrix to per-generator quadratic coefficients
`(c2, c1, c0)` (each a length-`ng` vector). Only the polynomial cost model
(model 2) is supported; coefficients are stored highest-degree first:

  ncost = 3 -> [c2, c1, c0]      ncost = 2 -> [0, c1, c0]      ncost = 1 -> [0, 0, c0]

Piecewise-linear (model 1) and cubic-or-higher (ncost > 3) are not representable
by the quadratic OPF objective and raise an error. When `gencost` has `2*ng`
rows (active then reactive cost), only the first `ng` (active) rows are used.
A missing `gencost` yields zero-cost generators (with a warning).
"""
function Matpower_Gencost_To_Quadratic(gencost, ng::Int)
    if gencost === nothing
        @warn "MATPOWER case has no mpc.gencost block; defaulting to zero generation cost."
        return zeros(ng), zeros(ng), zeros(ng)
    end

    nrows = size(gencost, 1)
    if nrows == 2 * ng
        @warn "mpc.gencost has 2*ng rows; using the first ng (active-power) rows and ignoring reactive-power cost."
    elseif nrows != ng
        throw(ArgumentError("mpc.gencost has $nrows rows; expected ng=$ng (or 2*ng=$(2ng))."))
    end

    c2 = zeros(ng); c1 = zeros(ng); c0 = zeros(ng)
    for g in 1:ng
        model = gencost[g, 1]
        model == 2 ||
            throw(ArgumentError("Generator $g uses MATPOWER cost model $(Int(model)); only model 2 " *
                                "(polynomial) is supported by the quadratic OPF objective."))
        ncost = Int(gencost[g, 4])                # number of polynomial coefficients
        coeffs = gencost[g, 5:5+ncost-1]          # highest degree first
        if ncost == 3
            c2[g], c1[g], c0[g] = coeffs[1], coeffs[2], coeffs[3]
        elseif ncost == 2
            c1[g], c0[g] = coeffs[1], coeffs[2]
        elseif ncost == 1
            c0[g] = coeffs[1]
        else
            throw(ArgumentError("Generator $g has ncost=$ncost (cubic or higher); the OPF objective " *
                                "is quadratic. Reduce the cost order to ≤ 3."))
        end
    end
    return c2, c1, c0
end

"""
    matpower_source_dataframes(parsed) -> (bus_df, gen_df, line_df)

Build the on-disk CSV layout (MATPOWER column order, ';'-headers, no `id`
column) consumed by the canonical readers. Extra trailing columns on `mpc.gen`
(e.g. 21-col rows) and `mpc.branch` are dropped; short `mpc.branch` rows are
padded with sensible defaults (angmin/angmax = ±360).
"""
function matpower_source_dataframes(parsed)
    bus     = parsed.bus
    gen     = parsed.gen
    branch  = parsed.branch
    ng      = size(gen, 1)

    # --- bus: 13 columns verbatim ----------------------------------------
    size(bus, 2) >= 13 ||
        throw(ArgumentError("mpc.bus has $(size(bus,2)) columns; expected ≥ 13."))
    bus_df = DataFrame(bus[:, 1:13], [:bus_i, :type, :Pd, :Qd, :Gs, :Bs, :area,
                                      :Vm, :Va, :baseKV, :zone, :Vmax, :Vmin])

    # --- gen: first 10 columns + quadratic cost --------------------------
    size(gen, 2) >= 10 ||
        throw(ArgumentError("mpc.gen has $(size(gen,2)) columns; expected ≥ 10."))
    c2, c1, c0 = Matpower_Gencost_To_Quadratic(parsed.gencost, ng)
    gen_df = DataFrame(gen[:, 1:10], [:bus, :Pg, :Qg, :Qmax, :Qmin, :Vg, :mBase,
                                      :status, :Pmax, :Pmin])
    gen_df.c2 = c2; gen_df.c1 = c1; gen_df.c0 = c0

    # --- branch: 13 columns (pad missing angmin/angmax) ------------------
    nbr_col = size(branch, 2)
    nbr_col >= 11 ||
        throw(ArgumentError("mpc.branch has $nbr_col columns; expected ≥ 11."))
    B = Matrix{Float64}(undef, size(branch, 1), 13)
    B[:, 1:min(13, nbr_col)] = branch[:, 1:min(13, nbr_col)]
    if nbr_col < 12; B[:, 12] .= -360.0; end       # angmin default
    if nbr_col < 13; B[:, 13] .=  360.0; end       # angmax default
    # Any NaN left by ragged rows in the optional tail → neutral defaults.
    for (j, default) in ((12, -360.0), (13, 360.0))
        @inbounds for i in axes(B, 1)
            isnan(B[i, j]) && (B[i, j] = default)
        end
    end
    line_df = DataFrame(B, [:fbus, :tbus, :r, :x, :b, :rateA, :rateB, :rateC,
                            :ratio, :angle, :status, :angmin, :angmax])

    return bus_df, gen_df, line_df
end

# ==============================================================================
#  3. Orchestrator
# ==============================================================================
"""
    Import_Matpower_Case(path_m; path_main, path_results)
        -> (; DBUS, DGEN, DCIR, baseMVA, path_names)

Parse a MATPOWER `.m` case, write `bus_data.csv` / `generators_data.csv` /
`line_data.csv` into a fresh timestamped `RESULTS/Results - <ts>/Inputs/` folder,
and return the in-house DataFrames re-read through the canonical `Read_Input_Data`
(so the schema, typing, `t_tap` fix and 1..nBUS bus remap are identical to every
other case).

`path_main` defaults to `project_root()`; `path_results` defaults to `default_results_dir()`.
"""
function Import_Matpower_Case(path_m::String;
        path_main::String = project_root(),
        path_results::String = default_results_dir())

    parsed = Parse_Matpower_File(path_m)
    bus_df, gen_df, line_df = matpower_source_dataframes(parsed)

    # Timestamped results tree (`Inputs/` only — no dispatch/TSC subfolders).
    path_names = build_import_results_paths(path_main, path_results)
    inputs_dir = path_names[:pf_inputs]

    CSV.write(joinpath(inputs_dir, "bus_data.csv"),        bus_df;  delim = ';')
    CSV.write(joinpath(inputs_dir, "generators_data.csv"), gen_df;  delim = ';')
    CSV.write(joinpath(inputs_dir, "line_data.csv"),       line_df; delim = ';')
    println("MATPOWER case '$(basename(path_m))' written as in-house CSVs to: $inputs_dir")

    # Re-read through the canonical reader → schema-perfect, bus-remapped frames.
    DBUS, DGEN, DCIR, _bus_map, _rev_map = Read_Input_Data(inputs_dir, false)

    println("Parsed MATPOWER case: $(nrow(DBUS)) buses, $(nrow(DGEN)) generators, " *
            "$(nrow(DCIR)) branches, baseMVA=$(parsed.baseMVA).")

    return (; DBUS, DGEN, DCIR, baseMVA = parsed.baseMVA, path_names)
end

# ==============================================================================
#  4. In-memory loader (for load_system() — no file I/O)
# ==============================================================================
"""
    matpower_to_inhouse_dataframes(path_m::String)
        -> (DBUS, DGEN, DCIR, baseMVA, bus_mapping, reverse_bus_mapping)

Pure in-memory MATPOWER → in-house DataFrame conversion. Applies the same
column renames, type casts, `t_tap` 0→1 fix, and 1..nBUS bus remap that
`Read_Input_Data` produces for the CSV path, but reads directly from a `.m`
file without writing any intermediate files.

Intended for use in `load_system()` when `RunConfig.matpower_file` is set.
`Import_Matpower_Case` (the write-to-disk orchestrator) is left unchanged.
"""
function matpower_to_inhouse_dataframes(path_m::String)
    # ── parse raw MATPOWER blocks ──────────────────────────────────────────────
    parsed = Parse_Matpower_File(path_m)
    bus_raw, gen_raw, line_raw = matpower_source_dataframes(parsed)

    nb   = nrow(bus_raw)
    ng   = nrow(gen_raw)
    ncir = nrow(line_raw)

    # ── bus DataFrame: positional rename mirrors Read_Input_Data/read_bus_data ─
    # MATPOWER col → in-house col:
    #   bus_i→:bus, type→:type, Pd→:p_d, Qd→:q_d, Gs→:g_sh, Bs→:b_sh,
    #   area→:area, Vm→:v_spe, Va→:v_a, baseKV→:base_kV, zone→:zone,
    #   Vmax→:v_max, Vmin→:v_min
    DBUS = DataFrame(
        bus     = Int64.(bus_raw.bus_i),
        type    = Int64.(bus_raw.type),
        p_d     = Float64.(bus_raw.Pd),
        q_d     = Float64.(bus_raw.Qd),
        g_sh    = Float64.(bus_raw.Gs),
        b_sh    = Float64.(bus_raw.Bs),
        area    = Int64.(bus_raw.area),
        v_spe   = Float64.(bus_raw.Vm),
        v_a     = Float64.(bus_raw.Va),
        base_kV = Float64.(bus_raw.baseKV),
        zone    = Int64.(bus_raw.zone),
        v_max   = Float64.(bus_raw.Vmax),
        v_min   = Float64.(bus_raw.Vmin),
    )

    # ── generator DataFrame: prepend :id, rename to in-house convention ────────
    # MATPOWER col → in-house col:
    #   bus→:bus, Pg→:pg_spe, Qg→:qg_spe, Qmax→:qg_max, Qmin→:qg_min,
    #   Vg→:vg_spe, mBase→:base_MVA, status→:g_status, Pmax→:pg_max,
    #   Pmin→:pg_min, c2→:g_cost_2, c1→:g_cost_1, c0→:g_cost_0
    DGEN = DataFrame(
        id       = Int64.(collect(1:ng)),
        bus      = Int64.(gen_raw.bus),
        pg_spe   = Float64.(gen_raw.Pg),
        qg_spe   = Float64.(gen_raw.Qg),
        qg_max   = Float64.(gen_raw.Qmax),
        qg_min   = Float64.(gen_raw.Qmin),
        vg_spe   = Float64.(gen_raw.Vg),
        base_MVA = Float64.(gen_raw.mBase),
        g_status = Int64.(gen_raw.status),
        pg_max   = Float64.(gen_raw.Pmax),
        pg_min   = Float64.(gen_raw.Pmin),
        g_cost_2 = Float64.(gen_raw.c2),
        g_cost_1 = Float64.(gen_raw.c1),
        g_cost_0 = Float64.(gen_raw.c0),
    )

    # ── circuit DataFrame: prepend :id, rename to in-house convention ──────────
    # MATPOWER col → in-house col:
    #   fbus→:from_bus, tbus→:to_bus, r→:l_res, x→:l_reac, b→:l_sh_susp,
    #   rateA→:l_cap_1, rateB→:l_cap_2, rateC→:l_cap_3,
    #   ratio→:t_tap, angle→:t_shift, status→:l_status,
    #   angmin→:ang_min, angmax→:ang_max
    t_tap_col = Float64.(line_raw.ratio)
    # MATPOWER convention: tap ratio = 0 means no transformer → treat as 1
    # (avoids zero-division in Ybus builders; mirrors Read_Input_Data behaviour).
    replace!(t_tap_col, 0.0 => 1.0)

    DCIR = DataFrame(
        id        = Int64.(collect(1:ncir)),
        from_bus  = Int64.(line_raw.fbus),
        to_bus    = Int64.(line_raw.tbus),
        l_res     = Float64.(line_raw.r),
        l_reac    = Float64.(line_raw.x),
        l_sh_susp = Float64.(line_raw.b),
        l_cap_1   = Float64.(line_raw.rateA),
        l_cap_2   = Float64.(line_raw.rateB),
        l_cap_3   = Float64.(line_raw.rateC),
        t_tap     = t_tap_col,
        t_shift   = Float64.(line_raw.angle),
        l_status  = Int64.(line_raw.status),
        ang_min   = Float64.(line_raw.angmin),
        ang_max   = Float64.(line_raw.angmax),
    )

    # ── bus remap to 1..nBUS (same logic as Read_Input_Data) ──────────────────
    bus_mapping, reverse_bus_mapping = Mapping_Buses_Labels(DBUS)
    DBUS.bus      = [bus_mapping[b] for b in DBUS.bus]
    DGEN.bus      = [bus_mapping[b] for b in DGEN.bus]
    DCIR.from_bus = [bus_mapping[b] for b in DCIR.from_bus]
    DCIR.to_bus   = [bus_mapping[b] for b in DCIR.to_bus]

    # Invariant: DBUS.bus must be exactly 1:nBUS after remap (many builders rely
    # on positional indexing by bus id). Fail loudly if violated.
    @assert DBUS.bus == collect(1:nrow(DBUS)) (
        "bus-remap invariant violated in matpower_to_inhouse_dataframes: " *
        "DBUS.bus must equal 1:nBUS after remap (got $(DBUS.bus)).")

    println("MATPOWER '$(basename(path_m))' loaded in-memory: " *
            "$(nb) buses, $(ng) generators, $(ncir) branches, baseMVA=$(parsed.baseMVA).")

    return DBUS, DGEN, DCIR, parsed.baseMVA, bus_mapping, reverse_bus_mapping
end
