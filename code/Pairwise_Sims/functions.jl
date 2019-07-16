##Helper functions/types for bacterial interaction simulations


#Types
struct Boltz
    B0::Float64
    E::Float64
end


struct Com
    U::Vector{Float64}
    γ::Array{Float64,2}
    R::Vector{Float64}
    N::Int64
end


#Functions

#Parameters
function boltz(B::Boltz,ΔT::Float64)
    return(B.B0 * exp(-B.E * ΔT))
end

function make_params_temp(U::Vector{Boltz},γ::Array{Boltz,2},R::Vector{Boltz},
                          N::Int64,ΔT::Float64 = 0.0)

    @assert all(N == size(U,1) == size(R,1) .== size(γ))

    U_T = boltz.(U,Ref(ΔT))
    γ_T = boltz.(γ,Ref(ΔT))
    R_T = boltz.(R,Ref(ΔT))

    return(Com(U_T,γ_T,R_T,N))
end

#Model
function dC(dC,C,p,t)
    dC[end] = 0.0 #set carbon resource inflow

    for i = 1:p.N
        dC[i] = C[i] * ( (p.U[i] * C[end]) - p.R[i] ) #uptake and respiration
        dC[end] -= C[end] * C[i] * p.U[i] #uptake from carbon source
        for j = 1:p.N
            dC[i] += C[end] * C[i] * C[j] * p.γ[i,j] #interactions
            dC[end] -= C[end] * C[i] * C[j] * p.γ[i,j]
        end
    end

    for i = 1:p.N+1
        if C[i] + dC[i] < eps()
            dC[i] = -C[i]
        end
    end

end

#Analysis
function get_biomass(sol,t_vec = nothing)
    if t_vec == nothing
        t_vec = sol.t
    end

    return sum(hcat(sol.(t_vec)...)[1:5,:],dims = 1)'
end

function get_respiration(sol,t_vec = nothing)
    if t_vec == nothing
        t_vec = sol.t
    end

    return  hcat(sol.(t_vec)...)[1:5,:]' * sol.prob.p.R
end
