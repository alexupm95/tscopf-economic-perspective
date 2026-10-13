# Loaded by ext/TSCOPFPlotsExt.jl only (not included from TSCOPF.jl).

using Plots
using Measures: mm
using DataStructures: OrderedDict

function _plots_save_tsred!(t_window_total::Vector{Float64},
    Pmt::OrderedDict{Int, Vector{Float64}},
    Pet::OrderedDict{Int, Vector{Float64}},
    Qet::Union{Nothing, OrderedDict{Int, Vector{Float64}}},
    δt::OrderedDict{Int, Vector{Float64}},
    δ_COIt::OrderedDict{Int, Vector{Float64}},
    δCOIt::Vector{Float64},
    Δωt::OrderedDict{Int, Vector{Float64}},
    Δω_COIt::OrderedDict{Int, Vector{Float64}},
    ΔωCOIt::Vector{Float64},
    time_RoCoF::Vector{Float64},
    RoCoF::OrderedDict{Int, Vector{Float64}},
    RoCoFCOI::Vector{Float64},
    Vke::OrderedDict{Int, Vector{Float64}},
    Pacc::OrderedDict{Int, Vector{Float64}},
    Pacc_COI::OrderedDict{Int, Vector{Float64}},
    PaccCOI::OrderedDict{Int, Vector{Float64}},
    Vpe::OrderedDict{Int, Vector{Float64}},
    base_MVA::Float64,
    f_syn::Float64,
    δ_tol::Tuple{Float64, Float64},
    path_names::OrderedDict{Symbol, String};
    t_clear_fault = nothing,
    # Machine-referenced corridor (`bound_style_δ ∈ (:highest_H, :ref_gen)`): the run bounds
    # δ_g − δ_ref, not δ_g − δ_COI, so the ±δ_tol band belongs on the reference machine and the
    # COI-relative figure is not the picture of the constraint. Both stay `nothing` on the
    # COI-referenced styles, where every plot below is unchanged.
    δ_ref_gen::Union{Nothing, Int} = nothing,
    δ_reft::Union{Nothing, OrderedDict{Int, Vector{Float64}}} = nothing,
    )

    pf_figures = path_names[:pf_TS_figures]
    mkpath(pf_figures)

    if !isnothing(t_clear_fault) index_t_ct = findfirst(x -> x >= t_clear_fault,  t_window_total) end # Index of the vector when the time is equal or greater than the clearing time (only for plotting purposes)

    # ====================================
    #              Plots
    # ====================================

    #----------------------------------------
    # Plot the electrical power for each generator in MW
    plot_Pe = plot()  # create an empty plot

    for (gen_id, values) in Pet
        plot!(plot_Pe,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Pe = plot!(xlabel="t (s)", ylabel="P_e (MW)", title="Electrical Power",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_Pe = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end

    plot_Qe = nothing
    if Qet !== nothing
        plot_Qe = plot()
        for (gen_id, values) in Qet
            plot!(plot_Qe,
                t_window_total,
                values,
                lw = 3,
                label = "G$gen_id",
                ls=:solid
            )
        end
        plot_Qe = plot!(xlabel="t (s)", ylabel="Q_e (MVAr)", title="Reactive Electrical Power",
        size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
        fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
        gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
        if !isnothing(t_clear_fault) plot_Qe = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end
    end

    #----------------------------------------
    # Plot rotor angles in degrees
    plot_δ = plot()  # create an empty plot

    for (gen_id, values) in δt
        plot!(plot_δ,
            t_window_total,
            rad2deg.(values),
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    # The corridor band is anchored on whatever the run actually referenced: the COI, or the
    # reference machine's own angle. Drawing the COI band on a machine-referenced run showed a
    # limit that was never imposed. The COI curve itself stays either way — it costs one line and
    # still tells the reader where the fleet's centre went.
    # Keyed off `δ_reft` rather than `δ_ref_gen`: the caller's helper already checked that the
    # reference machine has a trajectory, so this cannot raise a KeyError on `δt`.
    machine_referenced = δ_reft !== nothing
    δ_band_anchor = machine_referenced ? δt[δ_ref_gen] : δCOIt
    band_ref_label = machine_referenced ? "G$(δ_ref_gen)" : "COI"
    plot_δ = plot!(t_window_total, rad2deg.(δCOIt)                     , label="COI",   lw=2, ls=:solid, lc=:black)
    plot_δ = plot!(t_window_total, rad2deg.(δ_band_anchor .+ δ_tol[1]) , label="- tol (ref $band_ref_label)", lw=2, ls=:dash, lc=:red)
    plot_δ = plot!(t_window_total, rad2deg.(δ_band_anchor .+ δ_tol[2]) , label="+ tol (ref $band_ref_label)", lw=2, ls=:dash, lc=:red)
    plot_δ = plot!(xlabel="t (s)", ylabel="δ (deg)", title="Rotor Angles",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_δ = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end


    # Plot rotor angles in degrees in relation to the COI reference frame
    plot_δ_COI = plot()  # create an empty plot

    for (gen_id, values) in δ_COIt
        plot!(plot_δ_COI,
            t_window_total,
            rad2deg.(values),
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_δ_COI = plot!(xlabel="t (s)", ylabel="δ_COI (deg)", title="Rotor Angles (ref. COI)",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_δ_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end


    # Plot rotor angles in degrees in relation to the reference machine (machine-referenced
    # corridor only). The reference is flat at zero by construction — δ_ref − δ_ref — and the
    # ±δ_tol lines are horizontal here, so this figure is the corridor itself: a curve leaving
    # the band is an infeasible machine, read straight off the axes.
    plot_δ_ref = nothing
    if δ_reft !== nothing
        plot_δ_ref = plot()  # create an empty plot

        for (gen_id, values) in δ_reft
            plot!(plot_δ_ref,
                t_window_total,
                rad2deg.(values),
                lw = 3,
                label = "G$gen_id",
                ls=:solid
            )
        end
        plot_δ_ref = plot!([rad2deg(δ_tol[1])], seriestype=:hline, label="- tol", lw=2, ls=:dash, lc=:red)
        plot_δ_ref = plot!([rad2deg(δ_tol[2])], seriestype=:hline, label="+ tol", lw=2, ls=:dash, lc=:red)
        plot_δ_ref = plot!(xlabel="t (s)", ylabel="δ - δ_G$(δ_ref_gen) (deg)", title="Rotor Angles (ref. G$(δ_ref_gen))",
        size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
        fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
        gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
        if !isnothing(t_clear_fault) plot_δ_ref = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end
    end


    #----------------------------------------
    # Plot the speed deviation in p.u.
    plot_Δω = plot()  # create an empty plot

    for (gen_id, values) in Δωt
        plot!(plot_Δω,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Δω = plot!(t_window_total, ΔωCOIt              , label="COI",   lw=2, ls=:solid, lc=:black)
    plot_Δω = plot!(xlabel="t (s)", ylabel="Δω (p.u.)", title="Speed Deviation",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_Δω = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end


    # Plot speed deviation in relation to the COI reference frame
    plot_Δω_COI = plot()  # create an empty plot

    for (gen_id, values) in Δω_COIt
        plot!(plot_Δω_COI,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Δω_COI = plot!(xlabel="t (s)", ylabel="Δω_COI (p.u.)", title="Speed Deviation (ref. COI)",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_Δω_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end

    #--------------------------------------------
    # Plot the rotor speed in p.u.
    plot_Δω_1 = plot()  # create an empty plot

    for (gen_id, values) in Δωt
        plot!(plot_Δω_1,
            t_window_total,
            1 .+ values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Δω_1 = plot!(t_window_total, 1 .+ (ΔωCOIt)              , label="COI",   lw=2, ls=:solid, lc=:black)
    plot_Δω_1 = plot!(xlabel="t (s)", ylabel="Δω (p.u.)", title="Rotor Speed",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_Δω_1 = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end


    # Plot rotor speed in p.u. in relation to the COI reference frame
    plot_Δω_COI_1 = plot()  # create an empty plot

    for (gen_id, values) in Δω_COIt
        plot!(plot_Δω_COI_1,
            t_window_total,
            1 .+ values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Δω_COI_1 = plot!(xlabel="t (s)", ylabel="Δω_COI (p.u.)", title="Rotor Speed (ref. COI)",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_Δω_COI_1 = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end

    #--------------------------------------------
    # Plot the frequency in Hz
    plot_f = plot()  # create an empty plot

    for (gen_id, values) in Δωt
        plot!(plot_f,
            t_window_total,
            f_syn .* (1 .+ values),
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end

    plot_f = plot!(t_window_total, f_syn .* (1 .+ (ΔωCOIt))              , label="COI",   lw=2, ls=:solid, lc=:black)
    plot_f = plot!(xlabel="t (s)", ylabel="f (Hz)", title="Frequency",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_f = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end


    # Plot the frequency in relation to the COI in Hz
    plot_f_COI = plot()  # create an empty plot

    for (gen_id, values) in Δω_COIt
        plot!(plot_f_COI,
            t_window_total,
            f_syn .* (1 .+ values),
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_f_COI = plot!(xlabel="t (s)", ylabel="f (Hz)", title="Frequency (ref. COI)",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_f_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end

    #--------------------------------------------
    # Plot the RoCoF in Hz/s
    plot_RoCoF = plot()  # create an empty plot

    for (gen_id, values) in RoCoF
        plot!(plot_RoCoF,
            time_RoCoF,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end

    plot_RoCoF = plot!(time_RoCoF, RoCoFCOI, label="COI",   lw=2, ls=:solid, lc=:black)
    plot_RoCoF = plot!(xlabel="t (s)", ylabel="df/dt (Hz/sec)", title="RoCoF",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_RoCoF = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end


    #--------------------------------------------
    # Plot δ vs ω
    plot_δω = plot()  # create an empty plot

    for (gen_id, values) in δ_COIt
        plot!(plot_δω,
            rad2deg.(values),
            1.0 .+ Δω_COIt[gen_id],
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_δω = plot!(xlabel="δ (deg)", ylabel="Δω (p.u.)", title="Rotor speed vs Angle",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)

    if !isnothing(t_clear_fault)
        for (gen_id, values) in δ_COIt
            scatter!(plot_δω,
                [rad2deg.(values[index_t_ct])],
                [1.0 .+ Δω_COIt[gen_id][index_t_ct]],
                shape=:circle, 
                msize= 7, 
                color=:black, 
                label=""
            )
        end
    end

    #--------------------------------------------
    # Plot δ vs Pe and Pm
    plot_δ_COIPePm = plot()  # create an empty plot

    for (gen_id, values) in δ_COIt
        # Pe
        plot!(plot_δ_COIPePm,
            rad2deg.(values),
            Pet[gen_id],
            lw = 3,
            label = "Pe G$gen_id",
            ls=:solid
        )
        # Pm (constant or governor trajectory)
        plot!(plot_δ_COIPePm,
            rad2deg.(values),
            Pmt[gen_id],
            lw = 3,
            label = "Pm G$gen_id",
            ls=:dash
        )
    end
    
    plot_δ_COIPePm = plot!(xlabel="δ (deg)", ylabel="P (MW)", title="P vs Angle",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)

    #--------------------------------------------
    # Plot t vs Pacc 
    plot_Pacc = plot()  # create an empty plot
    # Pacc
    for (gen_id, values) in Pacc
        plot!(plot_Pacc,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Pacc = plot!(t_window_total, PaccCOI[1], label="COI", lw=2, ls=:solid, lc=:black)
    plot_Pacc = plot!(xlabel="t (s)", ylabel="P_acc (MW)", title="Accelerating Power",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_Pacc = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end

    #--------------------------------------------
    # Plot t vs Pacc (in relation to COI)
    plot_Pacc_COI = plot()  # create an empty plot
    # Pacc
    for (gen_id, values) in Pacc_COI
        plot!(plot_Pacc_COI,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Pacc_COI = plot!(xlabel="t (s)", ylabel="P_acc_COI (MW)", title="Accelerating Power (ref. COI)",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_Pacc_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end

    #--------------------------------------------
    # Plot t vs Vke
    plot_Vke = plot()  # create an empty plot
    # Vke
    for (gen_id, values) in Vke
        plot!(plot_Vke,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Vke = plot!(xlabel="t (s)", ylabel="V_ke (p.u. ⋅ rad)", title="Vke vs Time",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_Vke = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end


    #--------------------------------------------
    # Plot t vs Vpe
    plot_Vpe = plot()  # create an empty plot
    # Vpe
    for (gen_id, values) in Vpe
        plot!(plot_Vpe,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Vpe = plot!(xlabel="t (s)", ylabel="V_pe (p.u. ⋅ rad)", title="Vpe vs Time",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_Vpe = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end

    #--------------------------------------------
    # Plot t vs Ve (SG only — Vke omits GFM / H=0)
    plot_Ve = plot()  # create an empty plot
    # Ve
    for (gen_id, values) in Vpe
        haskey(Vke, gen_id) || continue
        plot!(plot_Ve,
            t_window_total,
            Vke[gen_id] .+ values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Ve = plot!(xlabel="t (s)", ylabel="V_e (p.u. ⋅ rad)", title="Total Energy vs Time",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault) plot_Ve = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct") end

    # ====================================
    #            Save Plots
    # ====================================
    savefig(plot_Pe,        joinpath(pf_figures, "Electrical Power vs time.svg"))
    if plot_Qe !== nothing
        savefig(plot_Qe, joinpath(pf_figures, "Reactive Electrical Power vs time.svg"))
    end
    savefig(plot_δ,         joinpath(pf_figures, "Delta vs time.svg"))
    savefig(plot_δ_COI,     joinpath(pf_figures, "Delta (ref COI) vs time.svg"))
    if plot_δ_ref !== nothing
        savefig(plot_δ_ref, joinpath(pf_figures, "Delta (ref G$(δ_ref_gen)) vs time.svg"))
    end
    savefig(plot_Δω,        joinpath(pf_figures, "Speed Deviation vs time.svg"))
    savefig(plot_Δω_COI,    joinpath(pf_figures, "Speed Deviation (ref COI) vs time.svg"))
    savefig(plot_Δω_1,      joinpath(pf_figures, "Rotor speed vs time.svg"))
    savefig(plot_Δω_COI_1,  joinpath(pf_figures, "Rotor speed (ref COI) vs time.svg"))
    savefig(plot_f,         joinpath(pf_figures, "Frequency vs time.svg"))
    savefig(plot_f_COI,     joinpath(pf_figures, "Frequency (ref COI) vs time.svg"))
    savefig(plot_RoCoF,     joinpath(pf_figures, "RoCoF vs time.svg"))
    savefig(plot_δω,        joinpath(pf_figures, "Rotor speed vs delta.svg"))
    savefig(plot_δ_COIPePm, joinpath(pf_figures, "Power vs delta.svg"))
    savefig(plot_Pacc,      joinpath(pf_figures, "Accelerating power vs time.svg"))
    savefig(plot_Pacc_COI,  joinpath(pf_figures, "Accelerating power (ref COI) vs time.svg"))
    savefig(plot_Vke,       joinpath(pf_figures, "Kinetic Energy vs time.svg"))
    savefig(plot_Vpe,       joinpath(pf_figures, "Potential Energy vs time.svg"))
    savefig(plot_Ve,        joinpath(pf_figures, "Total Energy vs time.svg"))

    println("Figures of the dynamic model successfully saved as SVG files in: ", pf_figures)
end

"""
    _plots_save_dual_ts_svgs!(t_window_total, families, pf_ts_figures_duals)

One SVG per generator per bound side, for every stability-corridor dual family the run built.

`families` is a vector of NamedTuples `(; title, side, tag, ylabel, data)` assembled by the
caller, which is the only place that knows which corridor style produced the multipliers. A
family whose `data` is `nothing` was not built and is skipped — that is how the δ-COI, δ-ref,
Δω-COI and Δω-abs corridors all reach this function without a style branch here.

`tag` disambiguates the filenames across families (`""` for the angle corridor, which keeps the
historical name; `"Domega "` for the speed one). Only one angle family and one speed family can
exist in a given run, so the untagged name never collides.
"""
function _plots_save_dual_ts_svgs!(
    t_window_total::Vector{Float64},
    families,
    pf_ts_figures_duals::String,
)
    # The pairwise corridor is keyed by a `(i, j)` tuple rather than a generator id; a raw
    # tuple would land in both the legend and the file name as "G(1, 2)".
    key_label(k::Tuple{Int, Int}) = "G$(k[1])-G$(k[2])"
    key_label(k) = "G$(k)"

    for fam in families
        fam.data === nothing && continue
        for (gen_id, values) in fam.data
            plt = plot(t_window_total, values,
                lw = 3, label = key_label(gen_id), ls=:solid,
                size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
                fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
                gridlinewidth=2, gridalpha=0.05, gridstyle=:dash
            )
            # Axis labels and a title naming the frame: without them these figures were
            # unreadable out of context, and COI-framed vs machine-framed duals are different
            # quantities that used to land in identically named files.
            plt = plot!(plt, xlabel="t (s)", ylabel=fam.ylabel, title=fam.title)
            savefig(plt, joinpath(pf_ts_figures_duals,
                "Duals Trans. Stab. Const $(fam.tag)$(key_label(gen_id)) $(fam.side).svg"))
        end
    end
    return nothing
end

"""
Per-series time plot.

`label_prefix` names the series ("Bus " for nodal quantities, "GFM " for converters).
`overlay` draws a second dict dashed on the same axes — used for raw-vs-clipped pairs,
where the gap between the two curves is the clip acting.
"""
function _plots_save_fullbus_svg!(
    t_window::Vector{Float64},
    data::OrderedDict{Int, Vector{Float64}},
    title::String,
    ylabel::String,
    pf_figures::String,
    filename::String;
    t_clear_fault=nothing,
    label_prefix::String="Bus ",
    overlay::Union{Nothing, OrderedDict{Int, Vector{Float64}}}=nothing,
    overlay_label_suffix::String=" raw",
)
    # Same contract as `_write_fullbus_csv`: a family spanning no units produces no file
    # rather than an empty axes object. Reachable on a converter-only fleet, where the
    # SG-only dq/AVR/governor families are present but hold no rows.
    isempty(data) && return nothing
    plt = plot()
    for (bus_id, values) in data
        plot!(plt, t_window, values, lw=3, label="$(label_prefix)$(bus_id)", ls=:solid)
    end
    if overlay !== nothing
        for (bus_id, values) in overlay
            plot!(plt, t_window, values, lw=2,
                label="$(label_prefix)$(bus_id)$(overlay_label_suffix)", ls=:dash)
        end
    end
    plt = plot!(xlabel="t (s)", ylabel=ylabel, title=title,
        size=(1200, 800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35),
        guidefont=font(35), legendfont=font(15), fontfamily="Times New Roman",
        left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
        gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    if !isnothing(t_clear_fault)
        plt = plot!(plt, [t_clear_fault], seriestype=:vline, line=:dash,
            color=:magenta, label="t_ct")
    end
    savefig(plt, joinpath(pf_figures, filename))
    return nothing
end
