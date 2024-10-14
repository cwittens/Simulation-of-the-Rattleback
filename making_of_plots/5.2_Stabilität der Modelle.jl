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

include(joinpath(@__DIR__, "../schoemer.jl"))
include(joinpath(@__DIR__, "../kuypers.jl"))



function plot_stability(with_schoemer_original, tol, method=Tsit5, ω0 = [0.01, -0.02, -2.0], tspan0 = (0.0, 60.0))


    #set up schoemer
    y_s, parameter_s = setup_rattleback_schoemer(ω0)
    prob_s = ODEProblem(rattleback!, y_s, tspan0, parameter_s)
    sol_s = solve(prob_s, method(), reltol=tol, abstol=tol);
    # calculate energy
    Eges_s, _, _, _ = calculate_energy_schoemer(sol_s, parameter_s);
    Eges_s = Eges_s .- first(Eges_s) # rescale


    #set up reduced schoemer
    y_rs, parameter_rs = setup_rattleback_schoemer_reduced(ω0)
    prob_rs = ODEProblem(rattleback_reduced!, y_rs, tspan0, parameter_rs)
    sol_rs = solve(prob_rs, method(), reltol=tol, abstol=tol);
    # calculate energy
    Eges_rs, _, _, _ = calculate_energy_schoemer_reduced(sol_rs, parameter_rs);
    Eges_rs = Eges_rs .- first(Eges_rs) # rescale



    # set up kuypers
    y0, parameter0 = set_up_original_problem_for_torque(ω0)
    prob_k = ODEProblem(dydt_torque!, y0, tspan0, parameter0)
    sol_k = solve(prob_k, method(), reltol=tol, abstol=tol); 
    # calculate energy
    Eges_k, _, _, _ = calculate_energy_kuypers(sol_k,parameter0)
    Eges_k = Eges_k .- first(Eges_k) # rescale


    lw = 2
    fontsize = 14
    fontsize2 = 18


    p = plot(sol_k.t, Eges_k, label = L"\Delta E  \: \: Kuypers", 
        xlabel = "t", ylabel = L"\Delta E"*" [a.u.]",   
        dpi = 300, size = (700, 500), 
        lw = lw,
        grid = true,
        titlefont = font(fontsize2,"Computer Modern"),
        guidefont = font("Computer Modern"),
        legendfont = font(fontsize,"Computer Modern"),
        xtickfontsize=fontsize,
        ytickfontsize=fontsize,
        ylabelfontsize=fontsize2,
        xlabelfontsize=fontsize2,
        xlims=(first(sol_k.t), last(sol_k.t)),
        yformatter = x -> @sprintf("%.0e", x),
    )
    if with_schoemer_original
        plot!(sol_s.t, Eges_s, label = L"\Delta E \: \: Schömer", lw=lw)
    end
    
    plot!(sol_rs.t, Eges_rs, label = L"\Delta E \: \: Schömer \; red.", lw=lw)
    return p
end

function plot_stability_different_methods(tol, method1, method2,  ω0 = [0.01, -0.02, -2.0], tspan0 = (0.0, 60.0))


    #set up reduced schoemer
    y_rs, parameter_rs = setup_rattleback_schoemer_reduced(ω0)
    prob_rs = ODEProblem(rattleback_reduced!, y_rs, tspan0, parameter_rs)
    sol_rs = solve(prob_rs, method1(), reltol=tol, abstol=tol);
    # calculate energy
    Eges_rs, _, _, _ = calculate_energy_schoemer_reduced(sol_rs, parameter_rs);
    Eges_rs = Eges_rs .- first(Eges_rs) # rescale


    sol_rs2 = solve(prob_rs, method2(), reltol=tol, abstol=tol);
    # calculate energy
    Eges_rs2, _, _, _ = calculate_energy_schoemer_reduced(sol_rs2, parameter_rs);
    Eges_rs2 = Eges_rs2 .- first(Eges_rs2) # rescale


    # set up kuypers
    y0, parameter0 = set_up_original_problem_for_torque(ω0)
    prob_k = ODEProblem(dydt_torque!, y0, tspan0, parameter0)
    sol_k = solve(prob_k, method1(), reltol=tol, abstol=tol); 
    # calculate energy
    Eges_k, _, _, _ = calculate_energy_kuypers(sol_k,parameter0)
    Eges_k = Eges_k .- first(Eges_k) # rescale


    sol_k2 = solve(prob_k, method2(), reltol=tol, abstol=tol); 
    # calculate energy
    Eges_k2, _, _, _ = calculate_energy_kuypers(sol_k2,parameter0)
    Eges_k2 = Eges_k2 .- first(Eges_k2) # rescale


    lw = 2
    fontsize = 14
    fontsize2 = 18


    p = plot(sol_k.t, Eges_k, label = L"Kuyper" * " mit $method1", xlabel = "t", ylabel = L"\Delta E"*" [a.u.]",   
        dpi = 300, size = (700, 500), 
        lw = lw,
        grid = true,
        titlefont = font(fontsize2,"Computer Modern"),
        guidefont = font("Computer Modern"),
        legendfont = font(fontsize,"Computer Modern"),
        xtickfontsize=fontsize,
        ytickfontsize=fontsize,
        ylabelfontsize=fontsize2,
        xlabelfontsize=fontsize2,
        xlims=(first(sol_rs.t), last(sol_rs.t)),
        legend=true,
        yformatter = x -> @sprintf("%.e", x),
    )
    plot!(sol_k2.t, Eges_k2, label = L"Kuyper" * " mit $method2", lw=lw)
    plot!(sol_rs.t, Eges_rs, label = L"Schömer \; red."*" mit $method1", lw=lw)
    plot!(sol_rs2.t, Eges_rs2, label = L"Schömer \; red."*" mit $method2", lw=lw)

    return p
end

# vgl schoemer, schoemer_reduced und kuypers
tol = 1e-5
with_schoemer_original = true
p1 = plot_stability(with_schoemer_original, tol, Tsit5)

save_at1 =  joinpath(@__DIR__, "../plots", "vergleich_energie_schoemer_schoemer_red_kuypers_tol_$tol.pdf")
@info "Plot comparing energies differences of schoemer, schoemer_red and kuypers at tol = $tol saved at", savefig(p1, save_at1)



############################################

# only look at schoemer_red and kuypers now
tol = 1e-5
with_schoemer_original = false
p_for_ylim = plot_stability(with_schoemer_original, tol, Tsit5) # berechne ohne schoemer_original um ylims zu bekommen und nutze diese für den plot mit schoemer_original
p2 = ylims!(p1, ylims(p_for_ylim))

save_at2 = joinpath(@__DIR__, "../plots", "vergleich_energie_schoemer_red_kuypers_tol_$tol.pdf")
@info "Plot comparing energies differences of schoemer_red and kuypers at tol = $tol saved at", savefig(p2, save_at2)


############################################

# look at differences between implicit and explicit methods
tol = 1e-10

# Tsit and Rodas
method1 = Tsit5
method2 = Rodas5


p3 = plot_stability_different_methods(tol, method1, method2)

save_at3 = joinpath(@__DIR__, "../plots", "vergleich_energie_implicit_methods1_tol_$tol.pdf")
@info "Plot comparing energies differences of schoemer_red and kuypers using $method1 and $method2 at tol = $tol saved at", savefig(p3, save_at3)


############################################


# DP5 and RadauIIA5
method1 = DP5
method2 = RadauIIA5


p4 = plot_stability_different_methods(tol, method1, method2)

save_at4 = joinpath(@__DIR__, "../plots", "vergleich_energie_implicit_methods2_tol_$tol.pdf")
@info "Plot comparing energies differences of schoemer_red and kuypers using $method1 and $method2 at tol = $tol saved at", savefig(p4, save_at4)

