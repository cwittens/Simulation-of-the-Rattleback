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

# pushed other direction
ω0 = [0.01, -0.02, 2.0]


#set up kuypers 
y0, parameter0 = set_up_original_problem_for_torque(ω0)

prob_k = ODEProblem(dydt_torque!, y0, tspan, parameter0)
sol_k = solve(prob_k, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);


p1 = plot_kuypers_omega(sol_k)

stored_at_1 =  joinpath(@__DIR__, "../plots", "kuypers_omega_pushed_other_direction.pdf")
@info "Pushing the rattleback in the other direction - using Kuypers. Saved at" , savefig(p1, stored_at_1)




#set up schoemer
ys, parameter_s = setup_rattleback_schoemer(ω0)
prob_s = ODEProblem(rattleback!, ys, tspan, parameter_s)
sol_s = solve(prob_s, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);


p2 = plot_schoemer_omega(sol_s)

stored_at_2 = joinpath(@__DIR__, "../plots", "schoemer_omega_pushed_other_direction.pdf")
@info "Pushing the rattleback in the other direction - using Schoemer. Saved at", savefig(p2, stored_at_2)





# pushed down
ω0 = [0.01, -1.0, -0.02]


#set up kuypers 
y0, parameter0 = set_up_original_problem_for_torque(ω0)

prob_k = ODEProblem(dydt_torque!, y0, tspan, parameter0)
sol_k = solve(prob_k, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);


p3 = plot_kuypers_omega(sol_k, :bottomright)



stored_at_3 =  joinpath(@__DIR__, "../plots", "kuypers_omega_pushed_down.pdf")
@info "Pushing the rattleback down - using Kuypers. Saved at", savefig(p3, stored_at_3)




#set up schoemer
ys, parameter_s = setup_rattleback_schoemer(ω0)
prob_s = ODEProblem(rattleback!, ys, tspan, parameter_s)
sol_s = solve(prob_s, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);


p4 = plot_schoemer_omega(sol_s, :bottomright)

stored_at_4 =  joinpath(@__DIR__, "../plots", "schoemer_omega_pushed_down.pdf")
@info "Pushing the rattleback down - using Schoemer. Saved at", savefig(p4, stored_at_4)






# pushed to hard

ω0 = [0.01, -10.0, -0.02]


#set up kuypers 
y0, parameter0 = set_up_original_problem_for_torque(ω0)

prob_k = ODEProblem(dydt_torque!, y0, tspan, parameter0)
sol_k = solve(prob_k, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);


p5 = plot_kuypers_omega(sol_k)

stored_at5 =  joinpath(@__DIR__, "../plots", "kuypers_omega_pushed_to_hard.pdf")
@info "Pushing the rattleback to hard - using Kuypers. Saved at", savefig(p5, stored_at5)



#set up schoemer
ys, parameter_s = setup_rattleback_schoemer(ω0)
prob_s = ODEProblem(rattleback!, ys, tspan, parameter_s)
sol_s = solve(prob_s, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01);


p6 = plot_schoemer_omega(sol_s)
ylims!(p6, ylims(p5))

stored_at6 =  joinpath(@__DIR__, "../plots", "schoemer_omega_pushed_to_hard.pdf")
@info "Pushing the rattleback to hard - using Schoemer. Saved at" , savefig(p6, stored_at6)



