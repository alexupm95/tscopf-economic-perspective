# ==============================================================================
#  Optional Plots.jl callbacks (registered by package extension TSCOPFPlotsExt)
# ==============================================================================

const _SAVE_TSRED_PLOTS_FN = Ref{Any}(nothing)
const _SAVE_DUAL_TS_SVG_FN = Ref{Any}(nothing)
const _SAVE_FULLBUS_SVG_FN = Ref{Any}(nothing)

"""`true` when `TSCOPFPlotsExt` has registered plotting callbacks."""
plots_extension_loaded() = _SAVE_TSRED_PLOTS_FN[] !== nothing

function register_plots_callbacks!(;
    save_tsred_plots!,
    save_dual_ts_constraint_svgs!,
    save_fullbus_svg!,
)
    _SAVE_TSRED_PLOTS_FN[] = save_tsred_plots!
    _SAVE_DUAL_TS_SVG_FN[] = save_dual_ts_constraint_svgs!
    _SAVE_FULLBUS_SVG_FN[] = save_fullbus_svg!
    return nothing
end

function _plots_extension_source_path()::String
    root = pkgdir(@__MODULE__)
    root === nothing &&
        throw(ArgumentError("load_plots_extension!(): TSCOPF pkgdir is unavailable."))
    return joinpath(root, "ext", "TSCOPFPlotsExt.jl")
end

"""
    _ensure_headless_gr_backend!()

Set `GKSwstype=100` (file-only GR output) when there is no display to draw on.

GR opens a window by default; on a machine with no display server that is an error
rather than a fallback, and it takes the whole run down at the first `savefig`. CI sets
this in the workflow `env:`, but a headless local run — an SSH session, a container, a
batch sweep — had nothing setting it. Never overrides an existing `GKSwstype`, and on
Windows (where GR always has a device) it is a no-op.
"""
function _ensure_headless_gr_backend!()
    haskey(ENV, "GKSwstype") && return nothing
    Sys.iswindows() && return nothing
    (haskey(ENV, "DISPLAY") || haskey(ENV, "WAYLAND_DISPLAY")) && return nothing
    ENV["GKSwstype"] = "100"
    @info "No display detected; setting GKSwstype=100 so GR writes figures to file only."
    return nothing
end

"""
    load_plots_extension!() -> Bool

Load Plots.jl and Measures.jl in `Main` so `TSCOPFPlotsExt` registers trajectory
callbacks. Safe to call multiple times. Returns `true` when
[`plots_extension_loaded`](@ref) is `true`.

`run_case!` calls this automatically when `save_ts_plots=true`.
"""
function load_plots_extension!()
    plots_extension_loaded() && return true
    _ensure_headless_gr_backend!()
    if !isdefined(Main, :Plots) || !isdefined(Main, :Measures)
        try
            Core.eval(Main, :(using Plots, Measures))
        catch err
            @warn "Plots.jl could not be loaded; trajectory figures will be skipped." exception=err
            return false
        end
    end
    if !plots_extension_loaded()
        try
            Base.retry_load_extensions()
        catch err
            @warn "TSCOPFPlotsExt auto-load failed; trying manual include." exception=err
        end
    end
    if !plots_extension_loaded()
        ext_path = _plots_extension_source_path()
        if isfile(ext_path) && !isdefined(Main, :TSCOPFPlotsExt)
            try
                Main.include(ext_path)
            catch err
                @warn "Manual TSCOPFPlotsExt include failed." exception=err
                return false
            end
        end
    end
    return plots_extension_loaded()
end

function _warn_plots_missing(feature::String)
    @warn "TSC trajectory plots require Plots.jl. Call `load_plots_extension!()` " *
          "(or `using Plots` after `using TSCOPF`) to enable $feature."
end

function _invoke_plots_callback!(fn, args...; kwargs...)
    return Base.invokelatest(fn, args...; kwargs...)
end

function Save_Dynamic_Results_Plots_tsred(args...; kwargs...)
    fn = _SAVE_TSRED_PLOTS_FN[]
    if fn === nothing
        _warn_plots_missing("trajectory figures")
        return nothing
    end
    return _invoke_plots_callback!(fn, args...; kwargs...)
end

function invoke_save_dual_ts_constraint_svgs!(args...; kwargs...)
    fn = _SAVE_DUAL_TS_SVG_FN[]
    if fn === nothing
        _warn_plots_missing("transient-stability dual SVG figures")
        return nothing
    end
    return _invoke_plots_callback!(fn, args...; kwargs...)
end

function _save_fullbus_svg(args...; kwargs...)
    fn = _SAVE_FULLBUS_SVG_FN[]
    if fn === nothing
        _warn_plots_missing("full-bus trajectory SVG figures")
        return nothing
    end
    return _invoke_plots_callback!(fn, args...; kwargs...)
end
