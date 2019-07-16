#Types
struct Boltz
    B0::Float64
    E::Float64
end

struct Param
    r::Float64
end

##Functions

#Parameters
function boltz(B::Boltz,ΔT::Float64)
    return(B.B0 * exp(-B.E * ΔT))
end

function make_params_temp(r::Boltz,ΔT::Float64 = 0.0)

    r_T = boltz(r,ΔT)

    return(Param(r_T))
end

function dC(dC,C,p,t)
    dC[1] = p.r * C[1]
end

function get_t_ϵ(p,C0,ϵ)
    return(log(ϵ/C0) / p.r)
end

function get_biomass(p,C0,t)
    C0 * exp(p.r*t)
end
