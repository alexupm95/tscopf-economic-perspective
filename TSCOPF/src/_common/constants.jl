# ==============================================================================
#  TSC-OPF — global numerical / physical constants
# ==============================================================================
# Single place to tune sentinel values and other shared literals used across
# dispatch, transient-stability, and I/O layers.  Included early from
# `src/TSCOPF.jl` so every entry point (main, tests, engine) sees the same
# definitions in `Main`.  Phase 6 will wrap this file in the package module and
# export the public constants for `using TSCOPF` workflows.
#
# Registry tables (dual export specs, variable order, …) stay in their own
# files — they are schema metadata, not user-tunable physical constants.
# ==============================================================================

# Magnitude, pu on base_MVA, of a solid three-phase short (|V_k| near 1e-5 on
# FULL_BUS). The default SC stamp is the inductive admittance
# `0 - im*FAULT_BUS_SHUNT`, added to the faulted-bus diagonal (`:add`): no active
# draw, residual angle near the pre-fault phasor. The resistive solid short
# remains available as the real value `FAULT_BUS_SHUNT + 0im` with `:overwrite`;
# that stamp is an artificial active-power sink during the fault and puts a small
# shadow price on the fault-window P balance. The old real 1e10 overwrite pushed
# |V_k| to ~1e-9 and scaled that bus's KCL rows so badly that the classical
# FULL_BUS 9-bus solve took 100+ Ipopt iterations with MUMPS reallocations, and
# ran past its 600 s limit on the Linux CI runners.
const FAULT_BUS_SHUNT = 1e6
