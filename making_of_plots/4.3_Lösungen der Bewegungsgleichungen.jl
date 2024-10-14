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


include(joinpath(@__DIR__, "../schoemer.jl"))
include(joinpath(@__DIR__, "../kuypers.jl"))

tol = 1e-12
tspan = (0.0, 2*34.35) # = (0.0, 2 * sol_k.t[argmin(abs.(sol_k[6,:]))]) to be symmetric 
#set up kuypers 
y0, parameter0 = set_up_original_problem_for_torque()
prob_k = ODEProblem(dydt_torque!, y0, tspan, parameter0)
sol_k = solve(prob_k, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);

p1 = plot_kuypers_omega(sol_k)

stored_at_1 =  joinpath(@__DIR__, "../plots", "vergleich_kuypers_omega.pdf")
@info "Kuypers Plot for comparing Kuypers and Schoemer saved at", savefig(p1,stored_at_1)




#set up schoemer
ys, parameter_s = setup_rattleback_schoemer()
prob_s = ODEProblem(rattleback!, ys, tspan, parameter_s)
sol_s = solve(prob_s, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);


p2 = plot_schoemer_omega(sol_s)
stored_at_2 =  joinpath(@__DIR__, "../plots", "vergleich_schoemer_omega.pdf")
@info "Schoemer Plot for comparing Kuypers and Schoemer saved at", savefig(p2, "$stored_at_2")
