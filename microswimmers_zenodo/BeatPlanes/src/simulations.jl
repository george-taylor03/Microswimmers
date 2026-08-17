abstract type Simulation end

struct TwoParameterSweep{S <: Simulation}
    param_1::Symbol
    param_2::Symbol
    values_1::Vector{Float64}
    values_2::Vector{Float64}
    results::Matrix{S}
end