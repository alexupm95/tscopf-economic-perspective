
# Function to calculate the DC power flow
# Auxiliar function to calculate the DC power flow at each branch
function Calculate_DC_Power_Flow(DCIR::DataFrame, nCIR::Int64, θ::Vector{Float64}, base_MVA::Float64, susc_model::SusceptanceModel = SIMPLE)

    Pik = zeros(Float64, nCIR)           # Initializing vector Pik
    Qik = zeros(Float64, nCIR)           # Initializing vector Qik
    Sik = zeros(Float64, nCIR)           # Initializing vector Sik
    Pki = zeros(Float64, nCIR)           # Initializing vector Pki
    Qki = zeros(Float64, nCIR)           # Initializing vector Qki
    Ski = zeros(Float64, nCIR)           # Initializing vector Ski
    Plosses = zeros(Float64, nCIR)       # Initializing vector Plosses
    Qlosses = zeros(Float64, nCIR)       # Initializing vector Qlosses
    circ_loading = zeros(Float64, nCIR)  # Initializing vector of percentage loading

    # Loop to calculate the power flow in the lines and transformers

    for lin = 1:nCIR
        if DCIR.l_status[lin] == 1 # Check if the branch is ON

            i       = DCIR.from_bus[lin]       # Bus i (from)
            k       = DCIR.to_bus[lin]         # Bus k (to)

            b     = dc_branch_susceptance(DCIR.l_res[lin], DCIR.l_reac[lin], susc_model)  # Series susceptance (SIMPLE 1/x | POWERMODELS x/(r²+x²))

            angik = θ[i] - θ[k]                # Angular difference between bus i and k
            angki = θ[k] - θ[i]                # Angular difference between bus k and i

            # Active power flow from i to k
            Pik[lin] = b * (angik)

            # Active power flow from k to i
            Pki[lin] = b * (angki)

            # Apparent power flow from i to k
            Sik[lin] = Pik[lin]

            # Apparent power flow from k to i
            Ski[lin] = Pki[lin] 

            circ_loading[lin] = DCIR.l_status[lin] * Sik[lin] / (DCIR.l_cap_1[lin] / base_MVA)

        end
    end

    return Pik, Qik, Sik, Pki, Qki, Ski, Plosses, Qlosses, circ_loading
end

# Function to calculate the DC power flow
# Auxiliar function to calculate the DC power flow at each branch
function Calculate_DC_Power_Flow(Bbus::SparseMatrixCSC, DCIR::DataFrame, nCIR::Int64, θ::Vector{Float64}, base_MVA::Float64)

    Pik = zeros(Float64, nCIR)           # Initializing vector Pik
    Qik = zeros(Float64, nCIR)           # Initializing vector Qik
    Sik = zeros(Float64, nCIR)           # Initializing vector Sik
    Pki = zeros(Float64, nCIR)           # Initializing vector Pki
    Qki = zeros(Float64, nCIR)           # Initializing vector Qki
    Ski = zeros(Float64, nCIR)           # Initializing vector Ski
    Plosses = zeros(Float64, nCIR)       # Initializing vector Plosses
    Qlosses = zeros(Float64, nCIR)       # Initializing vector Qlosses
    circ_loading = zeros(Float64, nCIR)  # Initializing vector of percentage loading

    # Loop to calculate the power flow in the lines and transformers

    for lin = 1:nCIR
        if DCIR.l_status[lin] == 1 # Check if the branch is ON

            i     = DCIR.from_bus[lin]       # Bus i (from)
            k     = DCIR.to_bus[lin]         # Bus k (to)

            B_ik  = Bbus[i,k]                # Series admittance

            angik = θ[i] - θ[k]                # Angular difference between bus i and k
            angki = θ[k] - θ[i]                # Angular difference between bus k and i

            # Active power flow from i to k
            Pik[lin] = B_ik * (angik)

            # Active power flow from k to i
            Pki[lin] = B_ik * (angki)

            # Apparent power flow from i to k
            Sik[lin] = Pik[lin]

            # Apparent power flow from k to i
            Ski[lin] = Pki[lin] 

            circ_loading[lin] = DCIR.l_status[lin] * Sik[lin] / (DCIR.l_cap_1[lin] / base_MVA)

        end
    end

    return Pik, Qik, Sik, Pki, Qki, Ski, Plosses, Qlosses, circ_loading
end

