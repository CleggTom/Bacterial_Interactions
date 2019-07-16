#Types
struct Boltz
    B0::Float64
    E::Float64
end

struct Param
    r::Float64
    K::Float64
end

##Functions

#Parameters
function boltz(B::Boltz,ΔT::Float64)
    return(B.B0 * exp(-B.E * ΔT))
end

function make_params_temp(r::Boltz,K::Boltz,ΔT::Float64 = 0.0)

    r_T = boltz(r,ΔT)
    K_T = boltz(K,ΔT)

    return(Param(r_T,K_T))
end

function dC(dC,C,p,t)
    dC[1] = p.r * C[1] * ( 1 - (C[1]/p.K) )
end

function get_t_exp(p,C0)
    A = (p.K - C0) / C0
    t_ϵ = - log(1/A)  / p.r
    return(t_ϵ)
end


function get_t_ϵ(p,C0,ϵ)
    A = (p.K - C0) / C0
    t_ϵ = - log((1/A)*(ϵ/(ϵ+A)) )  / p.r
    return(t_ϵ)
end

function get_t_τ(p,C0,τ,ϵ)
    A = (p.K - C0) / C0
    t_τ  = -τ * log(1/A * (ϵ/(ϵ+A)) ) / p.r
    return(t_τ)
end


function get_biomass(p,C0,t)
    A = (p.K - C0)/C0
    return (p.K) / (1 + (A * exp(-p.r * t)))
end

function get_biomass_τ(p,C0,τ,ϵ)
    A = (p.K - C0)/C0
    return (p.K) / (1 + exp(τ) * (ϵ / (ϵ + A)) )
end

```
Approximates the temperature sensitivity of a system in exponetial phase.
```
function exponential_approx(t,T)
    -Er*r0*t * exp(-Er*T)
end

```
    fit_E(t)

fits E values
```
function fit_E(t)
    biomass = vcat([res[i](t) for i = 1:Ntemp]...)
    @. model(x, p) = p[1] + (p[2] * x)
    p0 = [-1, 1.0]
    fit = curve_fit(model, collect(T_vec), log.(biomass), p0)
    p = coef(fit)
    return(fit)
end
