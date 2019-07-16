using Pkg
Pkg.activate(".")
using Plots, DifferentialEquations, LinearAlgebra, LsqFit

include("functions.jl")

r0 = .1; Er = 0.5 ; r = Boltz(r0,Er)
K0 = 10; EK = -0.5; K = Boltz(K0,EK)
C0 = 0.01

tspan = (0.0,1e6)
cb = CallbackSet(PositiveDomain())
Ntemp = 25; Ntime = 250
T_vec = range(-1.0,1.0,length = Ntemp)
res = Array{Any,1}(undef,Ntemp)

for i = 1:Ntemp
 p = make_params_temp(r,K,T_vec[i])
 prob = ODEProblem(dC,[C0],tspan,p)
 res[i] = solve(prob,callback = cb, alg_hints=[:auto])
end
#calculate times
t_equ = [get_t_ϵ(res[i].prob.p , C0 , eps()) for i = 1:Ntemp] #equilibrium
t_equ_ϵ = [get_t_ϵ(res[i].prob.p , C0 , 1e-10) for i = 1:Ntemp] #equilibrium - ϵ
t_exp = [get_t_exp(res[i].prob.p , C0) for i = 1:Ntemp] #exponetial

#plotting growth of all
t_vec = [range(eps(),t_equ[i],length = Ntime) for i = 1:Ntemp]
col_T = range(colorant"red",stop=colorant"blue",length = Ntemp)'
col_t = range(colorant"green",stop=colorant"blue",length = Ntime)'
# col = range(colorant"green",stop=colorant"blue",length = Ntime)'
biomass = hcat([hcat(res[i].(t_vec[i])...)' for i = 1:Ntemp]...)
p_all = plot(t_vec,biomass,linecolor = col_T,legend = false, size = (2000,1000), linewidth = 4)
vline!(t_exp) #adding expo point

#plot only pre exponential
t_vec = range(eps(),minimum(t_exp),length = Ntime)
biomass = hcat([hcat(res[i].(t_vec)...)' for i = 1:Ntemp]...)
p_pre = plot(T_vec,log.(biomass'),linecolor = col_t,legend = false, size = (2000,1000), linewidth = 4)

#intermediate period
t_vec = range(minimum(t_exp),maximum(t_equ_ϵ),length = Ntime)
biomass = hcat([hcat(res[i].(t_vec)...)' for i = 1:Ntemp]...)
p_int = plot(T_vec,log.(biomass'),linecolor = col_t,legend = false, size = (2000,1000), linewidth = 4)

#equilibrium period
t_vec = range(maximum(t_equ_ϵ),maximum(t_equ),length = Ntime)
biomass = hcat([hcat(res[i].(t_vec)...)' for i = 1:Ntemp]...)
p_post = plot(T_vec,log.(biomass'),linecolor = col_t,legend = false, size = (2000,1000), linewidth = 4)

l = @layout [a b c]
p_Temp = plot(p_pre,p_int,p_post, layout = l, ylim = )

l = @layout [a ; b]
plot(p_all,p_Temp,layout = l)

###
#Working out the  E in exponential
###
function fit_E(t)
    biomass = vcat([res[i](t) for i = 1:Ntemp]...)
    @. model(x, p) = p[1] + (p[2] * x)
    p0 = [-1, 1.0]
    fit = curve_fit(model, collect(T_vec), log.(biomass), p0)
    p = coef(fit)
    return(fit)
end

#get time scale to first texpo
t_vec_exp = range(eps(),minimum(t_exp),length = 100)
#caluculate biomass across this time
biomass = hcat([hcat(res[i].(t_vec_exp)...)' for i = 1:Ntemp]...)
#plot the evoluciton over time
p1 = plot(t_vec_exp,(biomass),legend = false, c = col_T)

#fit the tpcs at each timepoint
expo_fits = fit_E.(t_vec_exp)
#get the E values and the apporximation
E_expo = hcat(coef.(expo_fits)...)[2,:]
E_expo_app = exponential_approx.(t_vec_exp,0.0)
p2 = plot(t_vec_exp,E_expo)
plot!(p2,t_vec_exp,E_expo_app)
plot(p1,p2)

####
#Between E values
####
#
# function fit_E_sum(t)
#     biomass = vcat([res[i](t) for i = 1:Ntemp]...)
#     @. model(x, p) = (p[1]*exp(-p[2]*x)) - (p[3]*exp(-p[4]*x))
#     p0 = [0.0, 1.0, 0.0, 1.0]
#     lb = [0.0, 0.0, 0.0, 0.0]
#     ub = [Inf, 2.0, Inf, 2.0]
#
#     fit = curve_fit(model, collect(T_vec), (biomass), p0, lower = lb, upper = ub)
#     p = coef(fit)
#     return(fit)
# end
#
# function plot_fit(fits,i,biomass,T_vec)
#     p = coef(fits[i])
#     pred = p[1] .* exp.(-p[2] .* T_vec) .- p[3] .* exp.(-p[4] .* T_vec)
#     actual = biomass'[:,i]
#     p1 = plot(T_vec,pred)
#     scatter!(p1,T_vec,actual)
#     return(p1)
# end
#
#
# #get time scale to first texpo
# t_vec_int = range(minimum(t_exp),maximum(t_equ_ϵ),length = 100)
# #caluculate biomass across this time
# biomass = hcat([hcat(res[i].(t_vec_int)...)' for i = 1:Ntemp]...)
# #plot the evoluciton over time
# p1 = plot(t_vec_int,(biomass),legend = false, c = col_T)
# #fit the tpcs at each timepoint
# expo_fits = fit_E_sum.(t_vec_int)
#
# plot_fit(expo_fits,10,biomass,T_vec)
# plot(rss.(expo_fits))
# #get the E values and the apporximation
# E_expo = hcat(coef.(expo_fits)...)[2,:]
# E_expo_app = exponential_approx.(t_vec_int,0.0)
# p2 = plot(t_vec_int,E_expo)
# plot!(p2,t_vec_int,E_expo_app)
# plot(p1,p2)


####
# E in carrying capacity
####
t_vec_K = range(maximum(t_equ_ϵ),maximum(t_equ),length = 100)
#caluculate biomass across this time
biomass = hcat([hcat(res[i].(t_vec_K)...)' for i = 1:Ntemp]...)
#plot the evoluciton over time
p1 = plot(t_vec_K,log.(biomass),legend = false, c = col)

#fit the tpcs at each timepoint
expo_fits = fit_E.(t_vec_K)
#get the E values and the apporximation
E_expo = hcat(coef.(expo_fits)...)[2,:]
p2 = plot(t_vec_K,E_expo, ylim = (0.4,0.6))
plot!(p2,t_vec_K,fill(-EK,100))
plot(p1,p2)


#combine
#plotting growth
t_vec = [range(eps(),t_equ[i],length = Ntime) for i = 1:Ntemp]
col = range(colorant"red",stop=colorant"blue",length = Ntemp)'
# col = range(colorant"green",stop=colorant"blue",length = Ntime)'
biomass = hcat([hcat(res[i].(t_vec[i])...)' for i = 1:Ntemp]...)
p1 = plot(t_vec,biomass,linecolor = col,legend = false, size = (2000,1000), linewidth = 2)

#Evalues
t_vec_all = range(eps(),maximum(t_equ) + 100,length = Ntime)
fits_all = fit_E.(t_vec_all)
E_all    = hcat(coef.(fits_all)...)[2,:]
p2 = plot(t_vec_all,E_all)
plot!(p2,t_vec_exp,E_expo_app, lw = 4)
plot!(p2,t_vec_K,fill(-EK,100), lw = 4)

l = @layout [a ; b]

plot(p1,p2, layout = l)


#get time scale to first texpo
t_vec_exp = range(eps(),minimum(t_exp),length = 100)
#caluculate biomass across this time
biomass = hcat([hcat(res[i].(t_vec_exp)...)' ./ res[i].prob.p.K for i = 1:Ntemp]...)
#plot the evoluciton over time
p1 = plot(t_vec_exp,(biomass),legend = false, c = col_T)

#fit the tpcs at each timepoint
expo_fits = fit_E.(t_vec_exp)
#get the E values and the apporximation
E_expo = hcat(coef.(expo_fits)...)[2,:]
E_expo_app = exponential_approx.(t_vec_exp,0.0)
p2 = plot(t_vec_exp,E_expo)
plot!(p2,t_vec_exp,E_expo_app)
plot(p1,p2)



####
#Normalisation of the process
###

#Normalised to K
t_vec = [range(eps(),minimum(t_equ),length = Ntime) for i = 1:Ntemp]

biomass = hcat([hcat(res[i].(t_vec[i])...)' ./ boltz.(Ref(K),T_vec)[i] for i = 1:Ntemp]...)
plot(biomass ,legend = false, lc = col)


#Normalise to r time
t_vec = [range(eps(),minimum(t_equ),length = Ntime) for i = 1:Ntemp]
t_vec_n = [range(eps(),(t_equ[i]),length = Ntime) for i = 1:Ntemp]
t_τ = [get_t_τ.(Ref(res[i].prob.p),C0,range(eps(),1-eps(),length=Ntime),eps()) for i = 1:Ntemp]


biomass = hcat([hcat(res[i].(t_vec[i])...)' for i = 1:Ntemp]...)
p1 = plot(biomass ,legend = false, lc = col,linewidth = 2)
biomass = hcat([hcat(res[i].(t_vec_n[i])...)' for i = 1:Ntemp]...)
p2 = plot(biomass ,legend = false, lc = col,linewidth = 2)
biomass = hcat([hcat(res[i].(t_τ[i])...)' for i = 1:Ntemp]...)
p3 = plot(biomass ,legend = false, lc = col,linewidth = 2)

plot(p1,p2,p3)

n = 2
b_ana_real = get_biomass.(Ref(res[n].prob.p),C0,t_vec[n])
plot!(p1,b_ana_real , c = :black)

b_ana_num = get_biomass.(Ref(res[n].prob.p),C0,t_vec_n[n])
plot!(p2,b_ana_num , c = :black)

b_ana = get_biomass.(Ref(res[n].prob.p),C0,t_τ[n])
plot!(p3,b_ana , c = :black)

plot(p1,p2,p3)

get_biomass.(Ref(res[n].prob.p),C0,t_τ[n])
get_biomass_τ.(Ref(res[n].prob.p),C0,range(eps(),1-eps(),length=Ntime),eps())


t_test = get_t_ρ.(Ref(res[n].prob.p),C0,range(0.1,0.9,length = 250))
plot(collect(t_vec[n]),t_test)
plot!(t_vec[n],t_vec_n[n],legend = false , c = col)


#checking relative time calculation

get_t_ρ(res[1].prob.p,C0,0.999)

plot(biomass[:,1] ./ biomass[end,1])

plot(t_vec,t_vec_n,legend = false)
plot(T_vec,log.(t_equ),legend = false)

@. model(x, p) = p[1] + (p[2] * x)
p0 = [-1, 1.0]
fit = curve_fit(model, collect(T_vec), log.(t_equ), p0)
coef(fit)
