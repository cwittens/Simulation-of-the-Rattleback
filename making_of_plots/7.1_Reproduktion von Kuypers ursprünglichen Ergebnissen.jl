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

# reproduce original kuypers simulation to check the right implementation of the equations
# initial conditions
y0_kuypers, tspan_kuypers, parameter_kuypers = set_up_problem_kuypers_wackelstein_1()
prob_kuypers = ODEProblem(dydt_torque, y0_kuypers, tspan_kuypers, parameter_kuypers)
tol = 1e-9
sol_kuypers = solve(prob_kuypers, VCABM(), reltol=tol, abstol=tol);

P = plot_kuypers_wackelstein_1_Abb4(sol_kuypers, parameter_kuypers)

save_at = joinpath(@__DIR__, "../plots", "kuypers0_reproduced.pdf")
@info "Plot comparing my Kuypers implementation with the original Kuypers simulation saved at", savefig(P, save_at)


