function one_pass_cov(x::Vector{Float64},y::Vector{Float64})
    n = length(x)
    
    sum_xy = sum(x.*y)
    sum_x = sum(x)
    sum_y = sum(y)
    
    return (sum_xy - (sum_x*sum_y)/n)/n

end

function cov_test(x::Vector{Float64},y::Vector{Float64})
    mean(x .* y) - mean(x)*mean(y)
end

#define parameter structure
struct params
    r::Vector{Float64}
    R::Vector{Float64}
    a::Array{Float64,2}
    Nsp::Int8
end

#define GLV model
function dC!(dC,C,p,t)
    for i = 1:p.Nsp
        dC[i] = C[i] * p.r[i]
        for j = 1:p.Nsp
            dC[i] += C[i] * C[j] * p.a[i,j]
        end
    end
end



const k =  8.617e-5

kT(Temp) = 1 / (k * (Temp+273.15))

function boltz(B0,E,temp)
    B0 * exp((-E / 8.617e-5)  * ((1/(temp + 273.15))- (1/(20 + 273.15))) )
end


#approximate biomass solution
function approx(p,t,C0)
    ψ = mean(sum(a,dims = 2) ./ (Nsp - 1))
    
    top = mean(p.r) * mean(C0) * exp(mean(p.r) * t)
    bot = mean(p.r) - ((mean(C0) * ψ) * (1 + exp(mean(p.r) * t)))
    
    return( top / bot)
end

#approximate sensitvity
function app_E(C0,r0,E,Tdiff,t,ψ)
    r = r0*exp(E*Tdiff)
    top = (ψ * C0 * t * r) + (ψ*C0) + (t*r^2) - (ψ*C0*exp(t*r))
    bot = (ψ * C0 * exp(t*r)) - (ψ*C0) - (r)
    
    return( -E * (top/bot))
end


###
# analysis
###

#Respiration
function resp(sol,t)
    sol.prob.p.R .* sol(t)
end



###
#approximations
###

#define parameter structure
struct params_mean
    r::Vector{Float64}
    R::Vector{Float64}
    a::Array{Float64,2}
    ψ::Vector{Float64}
    Nsp::Int64
end

# Testing Approximations
function dC_meanfield!(dC,C,p::params_mean,t)
    for i = 1:p.Nsp
        dC[i] = (C[i] * p.r[i]) + (C[i] * C[i] * p.a[i,i]) + (C[i] * mean(C) * p.ψ[i])
    end
end

#average biomass
#for non-meanfield
function uC_equi(p::params)
    a = [p.a[i,i] for i = 1:p.Nsp]
    ψ = sum(p.a,dims=2) .- a
    -(mean(p.r) / mean(a)) * ( 1 / (1 + ( mean(ψ)/mean(a) ) ) )
end

#for Rv
function uC_equi(ur,ua,up)
    -(ur/ua) * (1 / (1 + (up/ua)))
end

#for meanfield
function uC_equi(p::params_mean)
    a = [p.a[i,i] for i = 1:p.Nsp]
    -(mean(p.r) / mean(a)) * ( 1 / (1 + ( mean(p.ψ)/mean(a) ) ) )
end

#equilibrial biomass
#for non-meanfield
function C_equi(p)
    a = [p.a[i,i] for i = 1:p.Nsp]
    ψ = sum(p.a,dims=2) .- a
    return (-p.r .- (ψ .* uC_equi(p))) ./ a
end



#equilibrial cov r-C
function cov_rc(p)
    a = [p.a[i,i] for i = 1:p.Nsp]
    
     (mean(p.r)^2 - mean(p.r .^ 2)) / mean(a)
end

###
#Simulation functions
###

#simulate a meanfield system and approximate average biomass
function meanfield_sim(n, ratio)
    Random.seed!(n)
    
    #Number of species
    Nsp = 100
    
    #Interactions
    a = -rand(Nsp,Nsp) ./ (Nsp) 
    [a[i,i] = 0.0 for i = 1:Nsp]
    
    #Intra
    intra = mean(sum(a,dims=2)) * ratio
    [a[i,i] = intra for i = 1:Nsp]
    
    #growth rate
    r = rand(Nsp)
    dR = Normal(rand()*10,100)
    R = rand(dR,Nsp)
    #params object
    p = params(r,R,a,Nsp)
    
    C0 = rand(Nsp) ./ 1000
    tspan = (0.0,50.0)
    tsave = tspan[1]:0.01:tspan[2]

    prob = ODEProblem(dC!, C0, tspan, p)
    sol = solve(prob ,saveat = tsave)

    mean_obs = mean(sol[end])
    mean_ana = uC_equi(p)
    
    return((mean_obs,mean_ana))
end

#Simulate at temperature ΔT
struct Temp_params
    r0::Array{Float64,1}
    Er::Array{Float64,1}
    R0::Array{Float64,1}
    ER::Array{Float64,1}
    
    ua::Float64
    sa::Float64
    ratio::Float64
    
    Nsp::Int64
    seed::Int64
    T::Float64
end

function sim_T(pT::Temp_params,C0, tspan, tsave)
    @assert Nsp == length(pT.r0) == length(pT.Er)
    
    #set seed
    Random.seed!(pT.seed)

    #Interactions
    da = Normal(pT.ua,pT.sa)
    a_T = rand(da,pT.Nsp,pT.Nsp)   
    a_T .+ boltz.(pT.r0,pT.Er,pT.T)./10
    [a_T[i,i] = pT.ratio for i = 1:pT.Nsp]
    
    
    #get r
    r_T = boltz.(pT.r0,pT.Er,pT.T)
    #resp
    R_T = boltz.(pT.R0,pT.ER,pT.T)
    
    p = params(r_T,R_T,a_T,pT.Nsp)
    prob = ODEProblem(dC!, C0, tspan, p)
    sol = DifferentialEquations.solve(prob , AutoTsit5(Rosenbrock23()), saveat = tsave)
    
    return(sol)
end

function pred_ΔT(pT::Temp_params)
    #get end biomass
    r_T = boltz.(pT.r0,pT.Er,pT.T)
    R_T = boltz.(pT.R0,pT.ER,pT.T)
    
    ur = mean(r_T)
    ua = -1
    up = pT.ua * pT.Nsp
    
    uC = uC_equi(ur,ua,up)
    
    return(sum(uC .* R_T))
end

function pred_ΔT(pT::Temp_params,sol)
    #get end biomass
    r_T = boltz.(pT.r0,pT.Er,pT.T)
    R_T = boltz.(pT.R0,pT.ER,pT.T)
    
    ur = mean(r_T)
    ua = -1
    up = pT.ua * pT.Nsp
    
    uC = uC_equi(sol.prob.p)
    
    return(sum(uC .* R_T))
end


#Pair simulations
struct pair_params
    r::Vector{Float64}
    a::Array{Float64,2}
end

function dC_pair!(dC,C,p,t)
    dC[1] = C[1] * (p.r[1] + p.a[1,1] * C[1] + p.a[1,2]*C[2])
    dC[2] = C[2] * (p.r[2] + p.a[2,2] * C[2] + p.a[2,1]*C[1])
end