using Pkg
Pkg.activate(".")
using Plots, DifferentialEquations, LinearAlgebra

include("functions.jl")

N = 5
U0 = fill(0.1,N) .+ rand(N)/10 ; R0 = rand(N) ./ 10
EU = fill(0.65,N)  ; ER = fill(0.65,N)

γ0 = rand(N,N)  .-  0.75
γ0[diagind(γ0)] .= -1.0
Eγ = fill(1.0,N,N)

U = Boltz.(U0,EU)
R = Boltz.(R0,ER)
γ = Boltz.(γ0,Eγ)

tspan = (0.0,1e6)
C0 = fill(0.01,N)
append!(C0,1.0)
# p = make_params_temp(U,γ,R,N,0.0)
# prob = ODEProblem(dC,C0,tspan,p)
# sol = solve(prob)

cb = CallbackSet(TerminateSteadyState())
res = Array{Any,1}(undef,10)
for i = 1:10
    p = make_params_temp(U,γ,R,N,range(-1.0,1.0,length=10)[i])
    prob = ODEProblem(dC,C0,tspan,p)
    res[i] = solve(prob,callback = cb, alg_hints=[:auto])
end

t_max = [res[i].t[end] for i = 1:10]
# t_vec = range.(eps(),t_max,length = 250)
t_vec = [range(eps(),5e2, length = 250) for i = 1:10]
col = range(colorant"red",stop=colorant"blue",length = 10)'

biomass = hcat(get_respiration.(res,t_vec)...)
plot(t_vec,biomass,linecolor = col,legend = false, size = (2000,1000))

#normalising to equilibrium time
t_scaled = t_vec .* (1 ./ t_max)
plot(t_scaled,biomass,linecolor = col,legend = false, size = (2000,1000))
