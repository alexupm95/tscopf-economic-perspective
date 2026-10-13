# Self-contained 3-generator UC toy case (swap with real CSVs via load_system later).

"""Minimal single-bus system for UC demos and tests."""
function build_uc_toy_system(; base_MVA::Float64 = 100.0, demand_MW::Float64 = 100.0)
    DBUS = DataFrame(
        bus = [1],
        type = [3],
        p_d = [demand_MW],
        q_d = [0.0],
        g_sh = [0.0],
        b_sh = [0.0],
        area = [1],
        v_spe = [1.0],
        v_a = [0.0],
        base_kV = [100.0],
        zone = [1],
        v_max = [1.1],
        v_min = [0.9],
    )

    DGEN = DataFrame(
        id = [1, 2, 3],
        bus = [1, 1, 1],
        pg_spe = [0.0, 0.0, 0.0],
        qg_spe = [0.0, 0.0, 0.0],
        qg_max = [0.0, 0.0, 0.0],
        qg_min = [0.0, 0.0, 0.0],
        vg_spe = [1.0, 1.0, 1.0],
        base_MVA = [base_MVA, base_MVA, base_MVA],
        g_status = [1, 1, 1],
        pg_max = [60.0, 80.0, 50.0],
        pg_min = [0.0, 10.0, 0.0],
        g_cost_2 = [0.0, 0.0, 0.0],
        g_cost_1 = [30.0, 50.0, 40.0],
        g_cost_0 = [0.0, 0.0, 0.0],
    )

    DCIR = DataFrame(
        id = Int[],
        from_bus = Int[],
        to_bus = Int[],
        l_res = Float64[],
        l_reac = Float64[],
        l_sh_susp = Float64[],
        l_cap_1 = Float64[],
        l_cap_2 = Float64[],
        l_cap_3 = Float64[],
        t_tap = Float64[],
        t_shift = Float64[],
        l_status = Int[],
        ang_min = Float64[],
        ang_max = Float64[],
    )

    bus_gen_circ_dict_ON = OrderedDict(
        1 => Dict(
            :gen_ids => [1, 2, 3],
            :gen_status => [true, true, true],
            :circ => Int[],
            :adj_buses => Int[],
            :vg_spe => 1.0,
            :pg_max => sum(DGEN.pg_max),
            :pg_min => sum(DGEN.pg_min),
            :qg_max => 0.0,
            :qg_min => 0.0,
            :pd_tot => demand_MW,
            :qd_tot => 0.0,
            :gsh_tot => 0.0,
            :bsh_tot => 0.0,
        ),
    )

    bus_mapping = OrderedDict(1 => 1)
    reverse_bus_mapping = OrderedDict(1 => 1)

    return (;
        DBUS, DGEN, DCIR,
        bus_gen_circ_dict_ON,
        bus_mapping, reverse_bus_mapping,
        base_MVA,
        nBUS = 1,
        nGEN = 3,
        nCIR = 0,
    )
end

"""Wrap the toy case as `SystemData` for `run_case!` (no INPUT_FILES folder)."""
function uc_toy_system_data(; kwargs...)
    toy = build_uc_toy_system(; kwargs...)
    return SystemData(
        toy.DBUS, toy.DGEN, nothing, nothing, nothing, toy.DCIR,
        toy.bus_mapping, toy.reverse_bus_mapping,
        toy.nBUS, toy.nGEN, toy.nCIR,
        toy.bus_gen_circ_dict_ON, toy.bus_gen_circ_dict_ON,
        "",  # no on-disk case folder — run_case! skips input archival
    )
end
