module TSCOPFPlotsExt

using TSCOPF
using Plots
using Measures: mm
using DataStructures: OrderedDict

const _PLOTS_IMPL = joinpath(dirname(@__DIR__), "src", "_manage_outputs", "_plots_impl_tsred.jl")
include(_PLOTS_IMPL)

function __init__()
    TSCOPF.register_plots_callbacks!(
        save_tsred_plots! = _plots_save_tsred!,
        save_dual_ts_constraint_svgs! = _plots_save_dual_ts_svgs!,
        save_fullbus_svg! = _plots_save_fullbus_svg!,
    )
    return nothing
end

end
