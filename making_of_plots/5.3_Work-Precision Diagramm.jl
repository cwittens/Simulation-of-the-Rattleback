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
using Printf
using DiffEqDevTools


include(joinpath(@__DIR__, "../schoemer.jl"))
include(joinpath(@__DIR__, "../kuypers.jl"))


# this code is based on this example:
# https://docs.sciml.ai/SciMLBenchmarksOutput/html/NonStiffODE/FitzhughNagumo_wpd.html


# define parameters for plotting
lw = 3
fontsize = 14
fontsize2 = 18

# define parameters for the problem
ω0 = [0.01, -0.02, -2.0]
tspan0 = (0.0, 2*34.35) # = (0.0, 2 * sol_k.t[argmin(abs.(sol_k[6,:]))]) to be symmetric

# define parameters for the work-precision runs
N_runs = 100

abstols = 1.0 ./ 10.0 .^ (3:13)
reltols = 1.0 ./ 10.0 .^ (3:13)

setups = [Dict(:alg=>Tsit5()),
          Dict(:alg=>DP5()),
          Dict(:alg=>Rodas5()),
          Dict(:alg=>RadauIIA5())
]


##############################
# kuypers
##############################

# set up kuypers
y0, parameter0 = set_up_original_problem_for_torque(ω0)
prob_k = ODEProblem(dydt_torque!, y0, tspan0, parameter0)



sol1 = solve(prob_k,Vern9(),abstol=1e-14,reltol=1e-14);
test_sol = TestSolution(sol1)



wp_l2 = WorkPrecisionSet(prob_k,abstols,reltols,setups;appxsol=test_sol,save_everystep=false,numruns=N_runs,maxiters=1000, error_estimate = :l2);

p1 = plot(wp_l2, xlabel = L"L^2"*" - Fehler", ylabel = "Rechenzeit [s]",   
    title = "Modell von Kuypers",
    titlefont = font(fontsize2,"Computer Modern"),
    guidefont = font("Computer Modern"),
    dpi = 300, size = (700, 500), 
    lw = lw,
    grid = true,
    legendfontsize=fontsize,
    legend = :topright,
    xtickfontsize=fontsize,
    ytickfontsize=fontsize,
    ylabelfontsize=fontsize2,
    xlabelfontsize=fontsize2,
    yformatter = x -> @sprintf("%.0e", x),
)


save_at1 =  joinpath(@__DIR__, "../plots", "work_precision_kuypers_l2.pdf")
@info "Work-Precision-Diagramm for Kuypers with l2 norm saved at", savefig(p1, save_at1)




wp_max = WorkPrecisionSet(prob_k,abstols,reltols,setups;appxsol=test_sol,save_everystep=false,numruns=N_runs,maxiters=1000, error_estimate = :l∞)

p2 = plot(wp_max, xlabel = L"L^\infty"*" - Fehler", ylabel = "Rechenzeit [s]",   
    title = "Modell von Kuypers",
    titlefont = font(fontsize2,"Computer Modern"),
    guidefont = font("Computer Modern"),
    #titlesize = fontsize2,
    dpi = 300, size = (700, 500), 
    lw = lw,
    grid = true,
    legendfontsize=fontsize,
    legend = :topright,
    xtickfontsize=fontsize,
    ytickfontsize=fontsize,
    ylabelfontsize=fontsize2,
    xlabelfontsize=fontsize2,
    yformatter = x -> @sprintf("%.0e", x),
)

save_at2 = joinpath(@__DIR__, "../plots", "work_precision_kuypers_max.pdf")
@info "Work-Precision-Diagramm for Kuypers with max norm saved at", savefig(p2, save_at2)



##############################
# reduced schoemer
##############################
#set up reduced schoemer
y_rs, parameter_rs = setup_rattleback_schoemer_reduced(ω0)
prob_rs = ODEProblem(rattleback_reduced!, y_rs, tspan0, parameter_rs)

sol_rs = solve(prob_rs,Vern9(),abstol=1e-14,reltol=1e-14);
test_sol_rs = TestSolution(sol_rs)




wp_l2_rs = WorkPrecisionSet(prob_rs,abstols,reltols,setups;appxsol=test_sol_rs,save_everystep=false,numruns=N_runs,maxiters=1000, error_estimate = :l2);

p3 = plot(wp_l2_rs, xlabel = L"L^2"*" - Fehler", ylabel = "Rechenzeit [s]",   
    title = "Modell von Schömer (red.)",
    titlefont = font(fontsize2,"Computer Modern"),
    guidefont = font("Computer Modern"),
    #titlesize = fontsize2,
    dpi = 300, size = (700, 500), 
    lw = lw,
    grid = true,
    legendfontsize=fontsize,
    legend = :topleft,
    xtickfontsize=fontsize,
    ytickfontsize=fontsize,
    ylabelfontsize=fontsize2,
    xlabelfontsize=fontsize2,
    yformatter = x -> @sprintf("%.0e", x),
)

save_at3 = joinpath(@__DIR__, "../plots", "work_precision_schoemer_red_l2.pdf")
@info "Work-Precision-Diagramm for Schoemer (red.) with l2 norm saved at", savefig(p3, save_at3)



wp_max_rs = WorkPrecisionSet(prob_rs,abstols,reltols,setups;appxsol=test_sol_rs,save_everystep=false,numruns=N_runs, maxiters=1000, error_estimate = :l∞)

p4 = plot(wp_max_rs, xlabel = L"L^\infty"*" - Fehler", ylabel = "Rechenzeit [s]",   
    title = "Modell von Schömer (red.)",
    titlefont = font(fontsize2,"Computer Modern"),
    guidefont = font("Computer Modern"),
    #titlesize = fontsize2,
    dpi = 300, size = (700, 500), 
    lw = lw,
    grid = true,
    legendfontsize=fontsize,
    legend = :topleft,
    xtickfontsize=fontsize,
    ytickfontsize=fontsize,
    ylabelfontsize=fontsize2,
    xlabelfontsize=fontsize2,
    yformatter = x -> @sprintf("%.0e", x),
)
save_at4 = joinpath(@__DIR__, "../plots", "work_precision_schoemer_red_max.pdf")
@info "Work-Precision-Diagramm for Schoemer (red.) with max norm saved at", savefig(p4, save_at4)
