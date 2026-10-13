# ===================================================================================
# Function to calculate the bus admittance matrix (Ybus) using the sparse COO method.
# Only the structural nonzeros are built (≈4 per branch + the per-bus shunts) and then
# assembled with SparseArrays.sparse, so the cost is O(nnz) in both time and memory —
# scalable to large systems, unlike a dense zeros(nBUS, nBUS) buffer. `sparse` sums
# duplicate (i, j) triplets, which correctly accumulates parallel branches and shunts.
# This is now the single Ybus builder (the dense Calculate_Ybus was removed).
# ===================================================================================
function Calculate_Ybus_sparse(DBUS::DataFrame, DCIR::DataFrame, nBUS::Int64, nCIR::Int64, base_MVA::Float64)
    # In MATPOWER and POWERMODELS, the TAP and SHIFT of the transformers are treated as "from" to "to"
    # This is the reason why they divide by TAP and not multiply when building the Ybus matrix
    # This approach is different from the one use during my undergraduate studies
    row_idx = Int[]
    col_idx = Int[]
    values = ComplexF64[]

    # Add line admittances
    for i in 1:nCIR
        k = DCIR.from_bus[i]
        m = DCIR.to_bus[i]

        if DCIR.l_status[i] == true
            ykm = 1 / (DCIR.l_res[i] + 1im * DCIR.l_reac[i])
            bkm_sh = DCIR.l_sh_susp[i] / 2
            tap = DCIR.t_tap[i]
            shift = deg2rad(DCIR.t_shift[i])

            # Ykk
            push!(row_idx, k)
            push!(col_idx, k)
            push!(values, (1/tap)^2 * ykm + 1im * bkm_sh)

            # Ykm
            push!(row_idx, k)
            push!(col_idx, m)
            push!(values, -(1/tap) * ykm * exp(1im * shift))

            # Ymk
            push!(row_idx, m)
            push!(col_idx, k)
            push!(values, -(1/tap) * ykm * exp(-1im * shift))

            # Ymm
            push!(row_idx, m)
            push!(col_idx, m)
            push!(values, ykm + 1im * bkm_sh)
        end
    end

    # Add shunt admittances from DBUS
    for i in 1:nBUS
        if DBUS.g_sh[i] != 0 || DBUS.b_sh[i] != 0
            push!(row_idx, i)
            push!(col_idx, i)
            push!(values, (DBUS.g_sh[i] + 1im * DBUS.b_sh[i]) / base_MVA)
        end
    end

    # Create the sparse matrix
    Ybus = sparse(row_idx, col_idx, values, nBUS, nBUS)
    return Ybus
end

# =========================================================================
# FULL_BUS dynamic-period admittance — NETWORK ONLY.
# Takes an already-built base network Ybus (Calculate_Ybus_sparse output = branches +
# bus shunts) and, for an SC bus fault, stamps a solid short on the faulted-bus
# diagonal (V ≈ 0). Loads are deliberately NOT folded in: in the FULL_BUS path the nodal
# power balance models them explicitly via the ZIP model (eq_const_fullbus_Pbalance!),
# so stamping load admittances here would double-count them. Post-fault and GL pass
# bus_fault=nothing (post-fault rebuilds the base on the tripped topology; GL reuses the
# pre-fault base, since Calculate_Ybus_sparse ignores p_d/q_d/g_status).
#
# `y_fault` / `stamp` apply only when `bus_fault` is set. The defaults are the
# inductive solid short `0 - im*FAULT_BUS_SHUNT` added to the network diagonal.
# The resistive option overwrites that diagonal with the real `FAULT_BUS_SHUNT`.
# =========================================================================
function Calculate_Ybus_fullbus_dynamics(Ybus_base::SparseMatrixCSC;
                                         bus_fault::Union{Nothing, Int64}=nothing,
                                         y_fault=-im * FAULT_BUS_SHUNT,
                                         stamp::Symbol=:add)
    Y = copy(Ybus_base)
    if bus_fault !== nothing
        _stamp_fault_admittance!(Y, bus_fault, y_fault, stamp)
    end
    return Y
end

# Stamp `y_fault` on one diagonal. `:add` accumulates, which is the default and
# what a realistic fault impedance must do so the network's own `Y_kk` is not
# dropped. `:overwrite` replaces the entry (the historical resistive stamp, when
# `y_fault` is the real `FAULT_BUS_SHUNT`).
function _stamp_fault_admittance!(Y::SparseMatrixCSC, bus_fault::Int, y_fault, stamp::Symbol)
    stamp in (:overwrite, :add) ||
        throw(ArgumentError("fault stamp must be :overwrite or :add (got $(stamp))."))
    if stamp == :overwrite
        Y[bus_fault, bus_fault] = y_fault
    else
        Y[bus_fault, bus_fault] += y_fault
    end
    return Y
end

# =========================================================================
# Fault-window start point for the shorted bus (FULL_BUS SC).
# With the fault shunt on the diagonal and no injection of its own, current balance at
# the shorted bus k gives  V̄_k ≈ −Σ_{j≠k} Y_kj V̄_j / Y_kk,  evaluated at the pre-fault
# phasors of the other buses. Two properties make this the seed to use:
#   * |V̄_k| scales as 1/|y_f| (≈ 1e-5 pu at |y_f| = 1e6 on 9-bus), so it follows the constant;
#   * the angle follows y_f. The default inductive shunt sits near the neighbours;
#     a real shunt sits near quadrature (≈ −90°). Seeding only a small |V_k| at the
#     pre-fault angle, under that real shunt, starts Ipopt ~90° away, next to the
#     polar singularity V ≈ 0 where the Q row has no slope in θ_k; that start stalled
#     the DQ_4TH + AVR + governor 9-bus case.
# Returns new Dicts; the caller's pre-fault val_V / val_θ are left untouched.
# =========================================================================
function fault_bus_voltage_seed(Ybus_fault::SparseMatrixCSC, val_V::Dict, val_θ::Dict,
                                bus_fault::Int)
    s = sum(Ybus_fault[bus_fault, j] * val_V[j] * cis(val_θ[j])
            for j in keys(val_V) if j != bus_fault)
    V_k = -s / Ybus_fault[bus_fault, bus_fault]
    return merge(val_V, Dict(bus_fault => abs(V_k))),
           merge(val_θ, Dict(bus_fault => angle(V_k)))
end

# =========================================================================
# Augment a bus-only network admittance with the classical generator internal
# nodes. Each active generator adds one internal node behind its transient
# reactance X′d (admittance y = 1/(jX′d)), coupled to its terminal bus, i.e.
# stamped through the incidence A_gen as A_gen·diag(y)·A_genᵀ. Returns the
# (nBUS+nGEN_ON)×(nBUS+nGEN_ON) sparse augmented matrix.
#
# Shared by Calculate_Ybus_fault_SC_Kron / _fault_LG / _postf_ClearFault (was copied
# verbatim in all three). Stays fully sparse: spdiagm for the gen admittances,
# blockdiag to embed the network block (no dense buffer, no CSC index-assignment).
# =========================================================================
function augment_with_gen_internals_Kron(Y_network::SparseMatrixCSC, DGEN_DYN::DataFrame,
                                    active_gen::Vector{Int64}, nBUS::Int64)
    nGEN_ON = length(active_gen)
    Y_gen = [1 / (1im * DGEN_DYN.Xd_tr[gen]) for gen in active_gen]   # y = -j / X′d

    # Generator incidence: +1 at the terminal bus, -1 at the internal node.
    A_gen = SparseArrays.sparse(
        vcat(DGEN_DYN.bus[active_gen], collect(1:nGEN_ON) .+ nBUS),   # row indices (terminal + internal)
        vcat(1:nGEN_ON, 1:nGEN_ON),                                   # column indices
        vcat(ones(nGEN_ON), -ones(nGEN_ON)),                          # +1 terminal, -1 internal
        nBUS + nGEN_ON, nGEN_ON)
    Y_gen_matrix = A_gen * SparseArrays.spdiagm(Y_gen) * A_gen'

    # Embed the network block top-left and add the generator contribution (all sparse).
    return SparseArrays.blockdiag(Y_network, SparseArrays.spzeros(ComplexF64, nGEN_ON, nGEN_ON)) +
           Y_gen_matrix
end

# =========================================================================
# Function to calculate the admittance matrix for the "during fault" period
# =========================================================================
function Calculate_Ybus_fault_SC_Kron(Ybus::SparseMatrixCSC, DBUS::DataFrame, DGEN::DataFrame, DGEN_DYN::DataFrame, nBUS::Int64, active_gen::Vector{Int64}, base_MVA::Float64, bus_fault::Int64;
                                      y_fault=-im * FAULT_BUS_SHUNT, stamp::Symbol=:add)

    # Work on a copy so the caller's base Ybus is never mutated (matches the
    # FULL_BUS path's Calculate_Ybus_fullbus_dynamics; lets callers build one base Ybus
    # and reuse it across fault windows instead of rebuilding it each time).
    Y = copy(Ybus)

    # Considering the loads as constant admittances
    for bus in eachindex(DBUS.bus)
        Y[bus, bus] = Y[bus, bus] + ((DBUS.p_d[bus] - 1im*DBUS.q_d[bus]) / base_MVA)
    end

    # Solid 3-phase short at the faulted bus. The default adds the inductive
    # shunt on top of the diagonal (network and the load just added). `:overwrite`
    # replaces that diagonal — the resistive option, real `FAULT_BUS_SHUNT`.
    _stamp_fault_admittance!(Y, bus_fault, y_fault, stamp)

    return augment_with_gen_internals_Kron(Y, DGEN_DYN, active_gen, nBUS)
end

# =======================================================================
# Function to calculate the admittance matrix for the "post-fault" period
# =======================================================================
# function Calculate_Ybus_postf(DBUS::DataFrame, DCIR::DataFrame, DGEN::DataFrame, DGEN_DYN::DataFrame, nBUS::Int64, nCIR::Int64, nGEN::Int64, base_MVA::Float64, l_status::Vector{Int64})

function Calculate_Ybus_postf_ClearFault_Kron(DBUS::DataFrame, DCIR::DataFrame, DGEN::DataFrame, DGEN_DYN::DataFrame, nBUS::Int64, nCIR::Int64, active_gen::Vector{Int64}, base_MVA::Float64)
    # Cleared-fault topology. Rebuild the network Ybus from the (post-fault) DCIR via
    # the single sparse builder — no dense zeros(nBUS, nBUS) buffer and no second copy
    # of the branch-stamping logic — then add the loads as constant admittances and the
    # generator internal nodes. Calculate_Ybus_sparse already includes the bus shunts.
    Y = Calculate_Ybus_sparse(DBUS, DCIR, nBUS, nCIR, base_MVA)

    # Loads as constant admittances on the bus diagonal
    for bus in 1:nBUS
        Y[bus, bus] = Y[bus, bus] + ((DBUS.p_d[bus] - 1im*DBUS.q_d[bus]) / base_MVA)
    end

    return augment_with_gen_internals_Kron(Y, DGEN_DYN, active_gen, nBUS)
end

# ===================================================================================================
# Function to calculate the admittance matrix considering the loss of load and/or generation in a bus
# ===================================================================================================
function Calculate_Ybus_fault_LG_Kron(Ybus::SparseMatrixCSC, DBUS::DataFrame, DGEN::DataFrame, DGEN_DYN::DataFrame, nBUS::Int64, active_gen::Vector{Int64}, base_MVA::Float64)

    # Work on a copy so the caller's base Ybus is never mutated (see the note in
    # Calculate_Ybus_fault_SC_Kron).
    Y = copy(Ybus)

    # Considering the loads as constant admittances
    for bus in DBUS.bus
        Y[bus, bus] = Y[bus, bus] + ((DBUS.p_d[bus] - 1im*DBUS.q_d[bus]) / base_MVA)
    end

    return augment_with_gen_internals_Kron(Y, DGEN_DYN, active_gen, nBUS)
end

# ============================================================================
# Per-branch DC susceptance "b" for the flow relation  P = b·Δθ.
#   SIMPLE      : b = 1/x                  (textbook DC; ignores series resistance)
#   POWERMODELS : b = x/(r²+x²) = -imag(1/(r+jx))   (PowerModels/MATPOWER variant)
# Both reduce to 1/x when r = 0. Centralising the choice here keeps every DC flow
# constraint and the post-solve flow report consistent with the selected model.
# ============================================================================
function dc_branch_susceptance(r::Real, x::Real, m::SusceptanceModel)
    return m == POWERMODELS ? x / (r^2 + x^2) : 1 / x
end

# ==========================================================================================
# Function to calculate the bus susceptance matrix  B = A · diag(b_0) · Aᵀ
#
# The per-branch kernel is the one selected by `susceptance_model`, evaluated through
# `dc_branch_susceptance` above so the matrix, the DC flow constraints and the post-solve
# flow report can never drift apart:
#   SIMPLE      : b_0 = -status / x                (textbook DC)
#   POWERMODELS : b_0 = -status · x / (r² + x²)    (imag(1/(r+jx)); equals -1/x when r = 0)
# ==========================================================================================
function Calculate_Matrix_B(DBUS::DataFrame, DCIR::DataFrame, nBUS::Int64, nCIR::Int64;
                            susceptance_model::SusceptanceModel = SIMPLE)
    # DBUS is the array related to the bus data
    # DCIR is the array related to the circuit data
    # nBUS is the number of buses
    # nCIR is the number of circuits

    # From the terminal nodes of each line (from_bus and to_bus), we create the incidence matrix,
    # where we assign 1 to from_bus nodes and -1 to to_bus nodes.
    # For the sparse function in SparseArrays, the arguments are:
    # sparse([Row Indices], [Column Indices], [Value], [Total Number of Rows], [Total Number of Columns])
    A = SparseArrays.sparse(DCIR.from_bus, 1:nCIR, 1, nBUS, nCIR) + SparseArrays.sparse(DCIR.to_bus, 1:nCIR, -1, nBUS, nCIR)

    # Create a vector with the susceptance values of each line (out-of-service branches drop
    # out through l_status). `Ref` keeps the enum scalar under broadcasting.
    B_0 = .- DCIR.l_status .* dc_branch_susceptance.(DCIR.l_res, DCIR.l_reac, Ref(susceptance_model))

    # Once we have the Incidence Matrix "A" and the Susceptance vector "B_0",
    # we can construct the Susceptance Matrix "B":
    B = A * SparseArrays.spdiagm(B_0) * A'
    # Here, spdiagm creates a sparse matrix and assigns the elements of vector B to the main diagonal

    # Return the susceptance matrix
    return B
end

# Thin alias kept because the PowerModels variant is referenced by name in the tests and docs.
Calculate_Matrix_B_PowerModels(DBUS::DataFrame, DCIR::DataFrame, nBUS::Int64, nCIR::Int64) =
    Calculate_Matrix_B(DBUS, DCIR, nBUS, nCIR; susceptance_model = POWERMODELS)

# =====================================================================================
# Function to calculate the Inverse of the Susceptance Matrix using Sparsity techniques
# =====================================================================================
function Calculate_Inverse_Matrix_B(B::SparseMatrixCSC, nBUS::Int64, ref_bus::Int64)
    S  = deepcopy(B)

    if !(ref_bus > 0 && ref_bus <= nBUS)
        throw(ArgumentError("invalid ref_bus in Calculate_Inverse_Matrix_B"))
    end

    S[ref_bus, :] .= 0.0
    S[:, ref_bus] .= 0.0
    S[ref_bus, ref_bus] = 1.0
    
    F = LinearAlgebra.ldlt(Symmetric(S); check=false)

    if !LinearAlgebra.issuccess(F)
        throw(ArgumentError("Failed factorization in Calculate_Inverse_Matrix_B"))
    end

    B_inv = F \ Matrix(1.0I, nBUS, nBUS)
    B_inv[ref_bus, :] .= 0.0  # zero-out the row of the slack bus
    
    # return Admittance Matrix Inverse
    return B_inv
end

# ================================================================
# Function to calculate the reduced matrix based on a given matrix
# ================================================================
function Reduce_Matrix(_matrix::SparseMatrixCSC, nBUS::Int64)
    if (size(_matrix)[1] ≤ nBUS)  && (size(_matrix)[2] ≤ nBUS )
        throw(ArgumentError("The size of the matrix must be greater than the number of buses."))
    end
    # Kron reduction = Schur complement eliminating the first nBUS (network) nodes:
    #   Y_red = Y22 − Y21 · Y11⁻¹ · Y12.
    # Keep the large network block Y11 SPARSE and factorise it with a sparse LU
    # (UMFPACK; works on ComplexF64). Only the thin RHS Y12 (nBUS × nGEN) is made dense
    # for the solve, so memory is O(nnz + nBUS·nGEN) instead of the O(nBUS²) that
    # densifying Y11 (Matrix(Y11) or inv(Y11)) would cost. Result is the small dense
    # nGEN × nGEN reduced matrix.
    Y11 = _matrix[1:nBUS, 1:nBUS]
    Y12 = _matrix[1:nBUS, nBUS+1:end]
    Y21 = _matrix[nBUS+1:end, 1:nBUS]
    Y22 = _matrix[nBUS+1:end, nBUS+1:end]
    F = LinearAlgebra.lu(Y11)                # sparse LU factorisation of the network block
    X = F \ Matrix(Y12)                      # X = Y11 \ Y12   (RHS is thin: nBUS × nGEN)
    red_matrix = Y22 - Y21 * X

    return red_matrix
end

# ======================================
# Function to save matrices in csv files
# ======================================
function Save_Matrix_CSV(path_main::String, path_results::String, file_name::String, _matrix)
    df_Ybus = DataFrame(Matrix(_matrix), :auto)        # Convert the admittance matrix into a DataFrame to save it
    CSV.write(joinpath(path_results, file_name*".csv"), df_Ybus; delim=';')    # Save the admittance matrix in a CSV file (absolute path)
end

# ======================================
# Function to save the the matrix in XLSX
# ======================================
function Save_Admittance_Matrix_XLSX(path_main::String, path_results::String, file_name::String, _matrix)
    mkpath(path_results)
    # 1. Define the file path (absolute; no cd needed)
    save_path = joinpath(path_results, file_name*".xlsx")
    
    # 2. Prepare the data
    # Convert to dense for Excel (Excel doesn't support sparse format)
    Y_dense = Matrix(_matrix)
    G = real.(Y_dense)
    B = imag.(Y_dense)
    
    # Create a string representation for the "Complex" sheet (e.g., "0.5 + 2.0im")
    Y_string = string.(Y_dense)

    # 3. Create DataFrames for each sheet
    # We use 'header=false' or generic names since Ybus is usually indexed by bus number
    df_complex = DataFrame(Y_string, :auto)
    df_real    = DataFrame(G, :auto)
    df_imag    = DataFrame(B, :auto)

    # 4. Write to XLSX with multiple sheets
    XLSX.writetable(save_path, 
        "Complex_Ybus" => df_complex, 
        "Real_G"       => df_real, 
        "Imaginary_B"  => df_imag;
        overwrite = true
    )

    println("Ybus successfully saved to: ", save_path)
end

# ======================================
# Function to save the the matrix in XLSX
# ======================================
function Save_Susceptance_Matrix_XLSX(path_main::String, path_results::String, file_name::String, _matrix)
    mkpath(path_results)
    # 1. Define the file path (absolute; no cd needed)
    save_path = joinpath(path_results, file_name*".xlsx")
    
    # 2. Prepare the data
    # Convert to dense for Excel (Excel doesn't support sparse format)
    B_dense = Matrix(_matrix)
    
    # Create a string representation
    B_string = string.(B_dense)

    # 3. Create DataFrames for each sheet
    # We use 'header=false' or generic names since Bbus is usually indexed by bus number
    df = DataFrame(B_string, :auto)

    # 4. Write to XLSX with multiple sheets
    XLSX.writetable(save_path, 
        "Susceptance Matrix" => df;
        overwrite = true
    )

    println("Susceptance Matrix successfully saved to: ", save_path)
end

# ==========================================================================
#  store_and_save_matrix!  —  Phase 2 (MF-4) infrastructure
# ==========================================================================
# Decouples the scenario-dependent matrix COMPUTATION (always needed: each
# fault / post-fault topology yields a different Ybus/Yred) from disk I/O
# (opt-in). Stores `matrix` in `matrices_dict[name]` unconditionally, and writes
# it to XLSX only when `save_matrices` is true, so a large δ_tol or contingency
# sweep can skip the redundant matrix dumps.
#
# `saver` selects the writer: Save_Admittance_Matrix_XLSX for complex Y matrices
# (default), Save_Susceptance_Matrix_XLSX for the real Bbus.
#
# NOTE: reserved for the Phase 3 fault redesign; the builders are not wired to it
# yet (their inline Save_*_XLSX calls are unchanged for now).
function store_and_save_matrix!(matrices_dict::AbstractDict,
                                path_names::OrderedDict,
                                name::String,
                                matrix;
                                save_matrices::Bool = true,
                                saver = Save_Admittance_Matrix_XLSX)
    matrices_dict[name] = matrix
    if save_matrices
        saver(path_names[:pf_main], path_names[:pf_bus_matrices], name, matrix)
    end
    return matrix
end

