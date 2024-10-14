using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using LinearAlgebra
using OrdinaryDiffEq
using StaticArrays
using ForwardDiff
using LaTeXStrings
using Plots
using BenchmarkTools
using Quaternionic

include(joinpath(@__DIR__, "../kuypers.jl"))
include(joinpath(@__DIR__, "../schoemer.jl"))
include(joinpath(@__DIR__, "../naive_lagrangian.jl"))
include(joinpath(@__DIR__, "../naive_hamilton.jl"))


tol = 1e-9
# inital conditions
y0_all_angles = [0.1, 0.1, 0.1, 0.2, 0.2, -1.0]


################################
#plot short range
################################


#kuypers 
_, parameter0 = set_up_original_problem_for_torque()
tspan_short = (0.0, 5.0) 


prob_k1 = ODEProblem(dydt_torque!, y0_all_angles, tspan_short, parameter0)
sol_k1 = solve(prob_k1, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);


p1 = plot_kuypers_omega(sol_k1)

# lagrange
y0_L, _, parameterL = set_up_original_problem_for_Lagrange(y0_all_angles)
prob_L1 = ODEProblem(dydt_Lagrange, y0_L, tspan_short, parameterL)
sol_L1 = solve(prob_L1, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);

p2 = plot_lagrange_omega(sol_L1, ylims(p1))


P_1 = plot(p1, p2, size = (1000, 400), dpi = 300, xlabel = "", ylabel = "")
saved_at1 = joinpath(@__DIR__, "../plots", "vergleich_kuypers_lagrange_short.pdf")
@info "Plot comparing Kuypers and Lagrange (short) saved at", savefig(P_1, saved_at1)



################################
#plot long range
################################

tspan_long = (0.0, 50.0)

# kuypers
prob_k2 = ODEProblem(dydt_torque!, y0_all_angles, tspan_long, parameter0)
sol_k2 = solve(prob_k2, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);

p3 = plot_kuypers_omega(sol_k2)



# lagrange
y0_L, _, parameterL = set_up_original_problem_for_Lagrange(y0_all_angles)
prob_L2 = ODEProblem(dydt_Lagrange, y0_L, tspan_long, parameterL)
sol_L2 = solve(prob_L2, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);

p4 = plot_lagrange_omega(sol_L2, ylims(p3))

P_2 = plot(p3, p4, size = (1000, 400), dpi = 300, xlabel = "", ylabel = "")

saved_at2 = joinpath(@__DIR__, "../plots", "vergleich_kuypers_lagrange_long.pdf")
@info "Plot comparing Kuypers and Lagrange (long) saved at", savefig(P_2, saved_at2)







################################
#plot long range Hamilton
################################


tspan_long = (0.0, 50.0)

# Hamilton
y0_H, _, parameterH = set_up_original_problem_for_Hamilton(y0_all_angles)
prob_H2 = ODEProblem(dydt_Hamilton_AD, y0_H, tspan_long, parameterH)
sol_H2 = solve(prob_H2, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01)

P_3 = plot_hamilton_omega(sol_H2, ylims(p3), parameterH)

save_at3 = joinpath(@__DIR__, "../plots", "vergleich_hamilton.pdf")
@info "Plot looking at Hamilton, to compair it to Lagrange (long) saved at", savefig(P_3, save_at3)