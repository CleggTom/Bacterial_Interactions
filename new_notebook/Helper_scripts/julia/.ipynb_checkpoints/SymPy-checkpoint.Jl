#this script calculates the aproximations for biomass in a GLV system

#exponential phase

###
#Step 1) Define Variables
###

#define symbols
aii, ua   = symbols("a_ii,μ_a")
t,T,N = symbols("t,T,N")
#Time dependent variables
C,uC,sC,x = symbols("C,μ_C,σ_C,x")
C = SymPy.Function(C)(t)
uC = SymPy.Function(uC)(t)
sC = SymPy.Function(sC)(t)

#non-dependent variables for lambdify
r_ , ur_, sr_ = symbols("r, μ_r,sigma_r")
C0,uC0,sC0 = symbols("C_0,μ_{C_0},σ_{C_0}")
ua_ = symbols("u_a")

#growth rates
r, ur, sr = symbols("r, μ_r, sigma_r")
r = SymPy.Function(r)(T)
ur = SymPy.Function(ur)(T)
sr = SymPy.Function(sr)(T) 
# ua = SymPy.Function(ua)(T)

###
#Step 2) get approximation
###

#get series approx
C_app = C.series(t,0,3).removeO()
logC_app = log(C).series(t,0,3).removeO()

###
#Step 3) get expressions for dC d2C and duC
###

#define dC system
dC = C*(r + (aii * C) + N*ua*uC)
dC2_1 = dC.diff(t)

#get avg approximaiton E[dC]
fμ = subs(dC,zip([r , C],[ur, uC])...) 

#second term requires dμC / dC
C_oth = symbols("C_{oth}")
dC_expanded = dC.subs(uC,(C+C_oth)/N)
duC = (fμ + (sC^2 / 2) * dC_expanded.diff(C,C))

#subs back into dC2
dC2 = dC2_1.subs(C.diff(t),dC).subs(uC.diff(t),duC)

###
#Step 4) substitute back into the approximations
###

###
#4.1) for non-log average biomass
###
f2μ = subs(dC2,zip([r , C],[ur, uC])...) 
f2μ

mean_C = uC.subs(t,0) + t * duC.subs(t,0) + t^2 * f2μ.subs(t,0)

r_ , ur_  = symbols("r, μ_r")
C0,uC0,sC0 = symbols("C_0,μ_{C_0},σ_{C_0}")

final_uC_sol = mean_C.subs(C.subs(t,0),C0).subs(uC.subs(t,0),uC0).subs(sC.subs(t,0),sC0).subs(r,r_).subs(ur,ur_)

# # #lambdify parameter order -> (N aii t ua ur uC sC)
uC_sol_func = lambdify(final_uC_sol)

###
# 4.2) Log-biomass 
###
dC_sub = logC_app.args[1].args[2].args[2].args[3].args[1]
dC2_sub = logC_app.args[1].args[2].args[1].args[3]

# #sub into the approximation
log_C_sol = logC_app.subs(dC_sub,dC.subs(t,0)).subs(dC2_sub,dC2.subs(t,0))

#sub dependent vars
into = [r,ur,ua,C.subs(t,0),uC.subs(t,0),sC.subs(t,0)]
outo = [r_,ur_,ua_,C0,uC0,sC0]

final_log_C_sol = subs(log_C_sol,zip(into,outo)...) 

# # #lambdify parameter order -> (C0 N aii r t ua ur uC0 sC0)
C_log_func = lambdify(final_log_C_sol)

#Equilibirum
