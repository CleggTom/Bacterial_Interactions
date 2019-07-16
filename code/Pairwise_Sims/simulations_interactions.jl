using Pkg
Pkg.activate(".")
using Plots, OrdinaryDiffEq, LinearAlgebra, LsqFit, Statistics

include("functions.jl")

N = 5
U0 = rand(N) ; R0 = rand(N) ./ 100
EU = fill(0.32,N)  ; ER = fill(0.65,N)
Eγ = fill(0.0,N,N)

U = Boltz.(U0,EU)
R = Boltz.(R0,ER)

tspan = (0.0,25)
C0 = fill(0.01,N)
append!(C0,1.0)

# cb = CallbackSet(TerminateSteadyState())

Nrep = 25
Nint = 25
Ntemp = 25
T_vec = range(-1.0,1.0,length=Ntemp)

res = Array{Any,2}(undef,Ntemp,Nrep*Nint)

for j = 1:Nrep*Nint
    γ0 = -rand(N,N)  .+  range(0,1,length = Nint)[Int(ceil(j / Nrep))]
    γ = Boltz.(γ0,Eγ)
    for i = 1:Ntemp
        print("\r : ",i," : ",j)
        p = make_params_temp(U,γ,R,N,T_vec[i])
        prob = ODEProblem(dC,C0,tspan,p)
        res[i,j] = solve(prob, alg_hints=[:auto], Tsit5())
    end
end

#No interactions
no_int_res = Vector{Any}(undef,Ntemp)
γ0 = fill(0.0,N,N)
γ = Boltz.(γ0,Eγ)
for i = 1:Ntemp
    p = make_params_temp(U,γ,R,N,T_vec[i])
    prob = ODEProblem(dC,C0,tspan,p)
    no_int_res[i] = solve(prob, alg_hints=[:auto],Tsit5())
end

#for a single temperature across interactions
T_ind = 5
res_single_T = res[T_ind,:]
t_vec = range.(eps(),tspan[2],length = 250)
col = repeat(range(colorant"black",stop=colorant"red",length = Nint),inner = Nrep)'
biomass = hcat(get_biomass.(res_single_T,Ref(t_vec))...)
p1 = plot(t_vec,biomass, linecolor = col,legend = false, size = (1000,500),xlab="Time",ylab="Biomass")
png(p1,"notebooks/figures/interaction_bio.png")

#respiration without interaction
resp = hcat(get_respiration.(res_single_T,Ref(t_vec))...)
resp_no_int = get_respiration(no_int_res[T_ind],Ref(t_vec))

p2 = plot(t_vec,resp,linecolor = col,legend = false, size = (1000,500),xlab="Time",ylab="Respiration")
plot!(p2,t_vec,resp_no_int,color = :blue, linewidth = 4)
png(p2,"notebooks/figures/interaction_resp.png")

#Looking across temperature
a_ind = 1; t_int = 2.0
res_single_int = res[:,a_ind]
t_vec = range.(eps(),tspan[2],length = 250)
col = range(colorant"red",stop=colorant"blue",length = Ntemp)'
resp = hcat(get_respiration.(res_single_int,Ref(t_vec))...)
p3 = plot(t_vec,resp,linecolor = col,legend = false, size = (1000,500),xlab="Time",ylab="Respiration")
vline!(p3,[t_int], color = :black, linewidth = 5)

a_ind = 1
res_single_int = res[:,a_ind]
t_vec = range.(eps(),t_int,length = 10)
col = range(colorant"black",stop=colorant"green",length = 10)'
resp = hcat(get_respiration.(res_single_int,Ref(t_vec))...)
p4 = plot(T_vec,resp',linecolor = col,legend = false, size = (1000,500),xlab="Temperature",ylab="Respiration")

png(plot(p3,p4),"notebooks/figures/temp_single.png")

##Across interactions
R_int = get_respiration.(res,Ref([t_int]))
R_int = map(x -> x[1] , R_int)

R_no_int = get_respiration.(no_int_res,Ref([t_int]))
R_no_int = map(x -> x[1] , R_no_int)

col = repeat(range(colorant"black",stop=colorant"red",length = Nint),inner = Nrep)'

p5 =plot(T_vec,(R_int),color = col, legend = false, size = (1000,500),xlab="Temperature",ylab="Respiration")
plot!(p5,T_vec,R_no_int,color = :blue, linewidth = 4)
png(p5,"notebooks/figures/temp_int.png")

#plotting E values
@. model(x, p) = p[1] + (p[2] * x)
p0 = [0.0, 1.0]
E_res = Array{Float64,2}(undef,Nrep*Nint,3)
for i = 1:Nrep*Nint
    #Evalue
    x = curve_fit(model, collect(T_vec), log.(R_int[:,i]), p0)
    E_res[i,1] = coef(x)[2]
    E_res[i,2] = rss(x)
    #average interaction
    E_res[i,3] = mean(res[13,i].prob.p.γ)
end

p6 =scatter(E_res[:,3],-E_res[:,1], legend = false, xlab = "Average Interaction Strength", ylab = "Respiration E value")
png(p6,"notebooks/figures/E_int")
