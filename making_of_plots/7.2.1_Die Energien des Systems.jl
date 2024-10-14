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

#set up kuypers 
tol = 1e-9
y0, parameter0 = set_up_original_problem_for_torque()
tspan0 = (0.0, 2*34.35) # = (0.0, 2 * sol_k.t[argmin(abs.(sol_k[6,:]))]) to be symmetric
prob_k = ODEProblem(dydt_torque!, y0, tspan0, parameter0)
sol_k = solve(prob_k, Tsit5(), reltol=tol, abstol=tol, saveat = 0.01); 


begin
    Eges, Ekin, Erot, Epot = calculate_energy_kuypers(sol_k,parameter0)

    epot_rescale = first(Epot)

    lw = 2
    fontsize = 14
    fontsize2 = 18


    p = plot(sol_k.t, Eges .- epot_rescale, label = L"E_{ges}", xlabel = "t", ylabel = "Energie [a.u.]",   
        dpi = 300, size = (700, 500), 
        lw = lw,
        grid = true,
        legend = :right,
        titlefont = font(fontsize2,"Computer Modern"),
        guidefont = font("Computer Modern"),
        legendfont = font(fontsize,"Computer Modern"),
        xtickfontsize=fontsize,
        ytickfontsize=fontsize,
        ylabelfontsize=fontsize2,
        xlabelfontsize=fontsize2,
        )

    plot!(sol_k.t, Epot .- epot_rescale, label = L"E_{pot}", alpha = 0.8, lw=lw)
    plot!(sol_k.t, Ekin, label = L"E_{trans}",  alpha = 0.8, lw=lw)
    plot!(sol_k.t, Erot, label = L"E_{rot}",  alpha = 0.8, lw=lw)
end

save_at = joinpath(@__DIR__, "../plots", "all_energies.pdf")
@info "Plot looking at all energies saved at" , savefig(p, save_at)


