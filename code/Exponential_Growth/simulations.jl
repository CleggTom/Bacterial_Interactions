using Pkg
Pkg.activate(".")
using Plots, DifferentialEquations, LinearAlgebra, LsqFit

include("functions.jl")

r0 = .1; Er = 0.5 ; r = Boltz(r0,Er)
C0 = 0.01

tspan = (0.0,10)
cb = CallbackSet(PositiveDomain())
Ntemp = 25; Ntime = 250
T_vec = range(-1.0,1.0,length = Ntemp)
res = Array{Any,1}(undef,Ntemp)

for i = 1:Ntemp
 p = make_params_temp(r,T_vec[i])
 prob = ODEProblem(dC,[C0],tspan,p)
 res[i] = solve(prob,callback = cb, alg_hints=[:auto])
end

#plotting growth
t_vec = [range(eps(),tspan[2],length = Ntime) for i = 1:Ntemp]
col = range(colorant"red",stop=colorant"blue",length = Ntemp)'
# col = range(colorant"green",stop=colorant"blue",length = Ntime)'
biomass = hcat([hcat(res[i].(t_vec[i])...)' for i = 1:Ntemp]...)
plot(t_vec,log.(biomass),linecolor = col,legend = false, size = (2000,1000), linewidth = 4)

###
#Working out the maximal E value
###

function fit_expo(t)
    biomass = vcat([res[i](t) for i = 1:Ntemp]...)
    @. model(x, p) = p[1] + (p[2] * x)
    p0 = [-1, 1.0]
    fit = curve_fit(model, collect(T_vec), log.(biomass), p0)
    p = coef(fit)
    return(fit)
end

expo_fits = fit_expo.(range(eps(),tspan[2],length = 100))
plot(hcat(coef.(expo_fits)...)[2,:])

#Approximation of exponential E values
function f(t,T)
    -Er*r0*t * exp(-Er*T)
end

x = [f.(range(eps(),tspan[2],length = 100),T_vec[i]) for i = 1:Ntemp]
plot!(hcat(x...),legend = false, c = :black)

plot!(f.(range(eps(),tspan[2],length = 100),0.0),legend = false, c = :red)
