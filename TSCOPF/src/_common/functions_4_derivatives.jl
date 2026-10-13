""" ===============================================================================================================
# Main References : 
[1] Steven C. Chapra and Raymond P. Canale, "Numerical Methods for Engineering", 5th edition, 2011, McGrawHill
# ==============================================================================================================="""

# =========================================================================================
# Function that returns the values in x and y axis of derivatives of two vectors
# It's based on numerical derivatives with data in unevenly spaced intervals
# Technique: Lagrange Polynomial Interpolator (see subsection 23.3, page 549, ref. [1])
# =========================================================================================
function Calculate_Derivatives_Lagrange_Interpolator(x_input::Vector{Float64}, y_input::Vector{Float64})
    """ Input parameters
    x_input = Number of points for the curve
    y_input = Mean value of the true distribution
    """
    # Check if the vectors x and y have the same length and if the length is at least 3
	if length(x_input) != length(y_input) ||  length(x_input) < 3
        throw(ArgumentError("Input vectors x and y must have the same length of at least 3."))
    end

    n = Int64(length(x_input))  # Number of elements of the vector x_input
    x_output = copy(x_input)    # Creates a copy of the vector input x
    x_output = x_output[2:n-1]  # Eliminate the first and last elements to perform the method

    y_output = zeros(Float64, n-2) # Initialize the vector that will receive the derivatives of the function related to point x
    aux_count = Int64(2) # Auxiliar counter
    for i in 1:n-2
        local aux1, aux2, aux3
        aux1 = y_input[aux_count-1] * ((2.0*x_input[aux_count] - x_input[aux_count] - x_input[aux_count+1])/((x_input[aux_count-1] - x_input[aux_count])*(x_input[aux_count-1] - x_input[aux_count+1])))
        aux2 = y_input[aux_count] * ((2.0*x_input[aux_count] - x_input[aux_count-1] - x_input[aux_count+1])/((x_input[aux_count] - x_input[aux_count-1])*(x_input[aux_count] - x_input[aux_count+1])))
        aux3 = y_input[aux_count+1] * ((2.0*x_input[aux_count] - x_input[aux_count-1] - x_input[aux_count])/((x_input[aux_count+1] - x_input[aux_count-1])*(x_input[aux_count+1] - x_input[aux_count])))
        y_output[i] =  aux1 + aux2 + aux3
        aux_count += 1
    end 

    return copy(x_output), copy(y_output) 
    """ Output parameters
    x_output = The value in x-axis of the derivatives (it's the same as x_input)
    y_output = The value in y-axis that represents the derivatives
    """
end

# =========================================================================================
# Function that returns the values in x and y axis of derivatives of two vectors
# It's based on numerical derivatives with data in evenly spaced intervals
# Technique: Centered finite divided difference (see subsection 23.1, page 547, ref. [1])
# =========================================================================================
function Calculate_Derivatives_CFDD(x_input::Vector{Float64}, y_input::Vector{Float64})
    """ Input parameters
    x_input = Number of points for the CDF curve
    y_input = Mean value of the true distribution
    """
    # Check if the vectors x and y have the same length and if the length is at least 5
	if length(x_input) != length(y_input) ||  length(x_input) < 5
        throw(ArgumentError("Input vectors x and y must have the same length of at least 5."))
    end
    # Check if the vector is evenly spaced (it's mandatory for this stencil).
    # The centred 5-point formula below uses a single step h = x[2]-x[1], so a
    # non-uniform grid would silently produce wrong derivatives.
    steps = diff(x_input)
    if iszero(steps[1]) || !all(isapprox.(steps, steps[1]; rtol=1e-9, atol=0.0))
        throw(ArgumentError("Input vector x is not evenly spaced."))
    end

    n = Int64(length(x_input))  # Number of elements of the vector x_input
    x_output = copy(x_input)    # Creates a copy of the vector input x
    x_output = x_output[3:n-2]  # Eliminate the first and last two elements to perform the method
    step_der = x_input[2] - x_input[1] # Calculate the step of derivation

    y_output = zeros(Float64, n-4) # Initialize the vector that will receive the derivatives of the function related to point x
    aux_count = Int64(3) # Auxiliar counter

    for i in 1:n-4
        y_output[i] = (-y_input[aux_count+2] + 8.0*y_input[aux_count+1] - 8.0*y_input[aux_count-1] + y_input[aux_count-2]) / (12.0*step_der)
        aux_count += 1
    end 

    return copy(x_output), copy(y_output) 
    """ Output parameters
    x_output = The value in x-axis of the derivatives (it's the same as x_input regarding two positions at the begining and at the end that were eliminated)
    y_output = The value in y-axis that represents the derivatives
    """
end

# =========================================================================================
# Function that returns the values in x and y axis of derivatives of two vectors
# Based on central finite differences
function Calculate_Derivative_CFD(x::Vector{<:Real}, y::Vector{<:Real})
    # Ensure the input vectors are of the same length
    @assert length(x) == length(y) "x and y must have the same length"
    @assert length(x) > 2 "vector length must be greater than 2"

    n = length(x)
    dydx = zeros(Float64, n)

    # Use central differences for interior points
    for i in 2:n-1
        dydx[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
    end

    # Use forward and backward difference at the ends
    dydx[1] = (y[2] - y[1]) / (x[2] - x[1])
    dydx[n] = (y[n] - y[n-1]) / (x[n] - x[n-1])

    return x, dydx
end
