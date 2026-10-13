# Function to calculate the AC power flow
# Auxiliar function to calculate the AC power flow at each branch
function Calculate_AC_Power_Flow(DCIR::DataFrame, nCIR::Int64, V::Vector{Float64}, θ::Vector{Float64}, base_MVA::Float64)
    # In MATPOWER and POWERMODELS, the TAP and SHIFT of the transformers are treated as "from" to "to"
    # This is the reason why they divide by TAP and the signal of the shift is switched

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

            i       = DCIR.from_bus[lin]                           # Bus i (from)
            k       = DCIR.to_bus[lin]                             # Bus k (to)

            yik     = 1 / (DCIR.l_res[lin] + 1im*DCIR.l_reac[lin]) # Series admittance
            bik_sh  = DCIR.l_sh_susp[lin] / 2                      # Shunt suscpetance
            g       = real(yik)                                    # Branch series conductance
            b       = imag(yik)                                    # Branch series susceptance

            t_tap   = DCIR.t_tap[lin]                              # Transformer tap ratio (tap:1)
            t_shift = deg2rad(DCIR.t_shift[lin])                   # Transformer shift angle

            angik = θ[i] - θ[k]                                    # Angular difference between bus i and k
            angki = θ[k] - θ[i]                                    # Angular difference between bus k and i

            # Active power flow from i to k
            Pik[lin] = g*((1/t_tap) * V[i])^2 - ((1/t_tap) * V[i] * V[k]) * (g*cos(angik - t_shift) + b*sin(angik - t_shift))

            # Reactive power flow from i to k
            Qik[lin] = -(b + bik_sh)*((1/t_tap) * V[i])^2 + ((1/t_tap) * V[i] * V[k]) * (b*cos(angik - t_shift) - g*sin(angik - t_shift))

            # Active power flow from k to i
            Pki[lin] = g*(V[k]^2) - ((1/t_tap) * V[i] * V[k]) * (g*cos(angki + t_shift) + b*sin(angki + t_shift))

            # Reactive power flow from k to i
            Qki[lin] = -(b + bik_sh)*(V[k]^2) + ((1/t_tap) * V[i] * V[k]) * (b*cos(angki + t_shift) - g*sin(angki + t_shift))

            # Apparent power flow from i to k
            Sik[lin] = abs(Pik[lin] + 1im*Qik[lin])

            # Apparent power flow from k to i
            Ski[lin] = abs(Pki[lin] + 1im*Qki[lin])

            # Active power losses in the branch
            Plosses[lin] = Pik[lin] + Pki[lin]

            # Reactive power losses in the branch
            Qlosses[lin] = Qik[lin] + Qki[lin]

            circ_loading[lin] = DCIR.l_status[lin] * abs(max(Sik[lin], Ski[lin])) / (DCIR.l_cap_1[lin] / base_MVA)

        end
    end

    return Pik, Qik, Sik, Pki, Qki, Ski, Plosses, Qlosses, circ_loading
end

# Function to calculate the AC power flow
# Auxiliar function to calculate the AC power flow at each branch
function Calculate_AC_Power_Flow(Ybus::SparseMatrixCSC, DCIR::DataFrame, nCIR::Int64, V::Vector{Float64}, θ::Vector{Float64}, base_MVA::Float64)

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

            i       = DCIR.from_bus[lin]                           # Bus i (from)
            k       = DCIR.to_bus[lin]                             # Bus k (to)

            G_ik    = real(Ybus[i,k])                              # Conductance
            B_ik    = imag(Ybus[i,k])                              # Susceptance

            angik = θ[i] - θ[k]                                    # Angular difference between bus i and k
            angki = θ[k] - θ[i]                                    # Angular difference between bus k and i

            # Active power flow from i to k
            Pik[lin] = V[i] * V[k] * (G_ik*cos(angik) + B_ik*sin(angik))

            # Reactive power flow from i to k
            Qik[lin] = V[i] * V[k] * (G_ik*sin(angik) - B_ik*cos(angik))

            # Active power flow from k to i
            Pki[lin] = V[i] * V[k] * (G_ik*cos(angki) + B_ik*sin(angki))

            # Reactive power flow from k to i
            Qki[lin] = V[i] * V[k] * (G_ik*sin(angki) - B_ik*cos(angki))

            # Apparent power flow from i to k
            Sik[lin] = abs(Pik[lin] + 1im*Qik[lin])

            # Apparent power flow from k to i
            Ski[lin] = abs(Pki[lin] + 1im*Qki[lin])

            # Active power losses in the branch
            Plosses[lin] = Pik[lin] + Pki[lin]

            # Reactive power losses in the branch
            Qlosses[lin] = Qik[lin] + Qki[lin]

            circ_loading[lin] = DCIR.l_status[lin] * abs(max(Sik[lin], Ski[lin])) / (DCIR.l_cap_1[lin] / base_MVA)

        end
    end

    return Pik, Qik, Sik, Pki, Qki, Ski, Plosses, Qlosses, circ_loading
end

