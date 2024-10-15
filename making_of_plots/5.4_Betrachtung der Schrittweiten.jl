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


function plot_stability_and_ratio_with_konstant_kuypers(tol, method=Tsit5, ω0 = [0.01, -0.02, -2.0], tspan0 = (0.0, 60.0))

    # set up kuypers
    y0, parameter0 = set_up_original_problem_for_torque(ω0)
    prob_k = ODEProblem(dydt_torque!, y0, tspan0, parameter0)
    sol_k = solve(prob_k, method(), reltol=tol, abstol=tol); 
    # calculate energy
    Eges_k, _, _, _ = calculate_energy_kuypers(sol_k,parameter0)
    Eges_k = Eges_k .- first(Eges_k) # rescale




    # set up kuypers with konstante Schrittweite
    dt = minimum(diff(sol_k.t[10:end-10])) # take the minimum stepsize (ignoring start and end)
    sol_k2 = solve(prob_k, method(), adaptive = false, dt = dt);
    # calculate energy
    Eges_k2, _, _, _ = calculate_energy_kuypers(sol_k2,parameter0)
    Eges_k2 = Eges_k2 .- first(Eges_k2) # rescale


    lw = 2
    fontsize = 14
    fontsize2 = 18


    p = plot(sol_k.t, Eges_k, label = L"\Delta E"*" Kuypers", xlabel = "t", ylabel = L"\Delta E"*" [a.u.]",   
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

    plot!(sol_k2.t, Eges_k2, label = L"\Delta E"*" Kuypers konst.", lw=lw)
    n = sol_k.stats.naccept + sol_k.stats.nreject
    n2 =length(sol_k2) - 1
    verhaeltnis = n2/n
    # println("constant method used $verhaeltnis times more steps")
    return p, verhaeltnis
end

function plot_stability_and_ratio_with_konstant_schoemer_red(tol, method=Tsit5, ω0 = [0.01, -0.02, -2.0], tspan0 = (0.0, 60.0))
    
    #set up reduced schoemer
    y_rs, parameter_rs = setup_rattleback_schoemer_reduced(ω0)
    prob_rs = ODEProblem(rattleback_reduced!, y_rs, tspan0, parameter_rs)
    sol_rs = solve(prob_rs, method(), reltol=tol, abstol=tol);
    # calculate energy
    Eges_rs, _, _, _ = calculate_energy_schoemer_reduced(sol_rs, parameter_rs);
    Eges_rs = Eges_rs .- first(Eges_rs) # rescale


    # set up reduced schoemer with konstante Schrittweite
    dt = minimum(diff(sol_rs.t[10:end-10])) # take the minimum stepsize (ignoring start and end)
    # dt_mean = mean(diff(sol_rs.t))
    # println("mean zu kleinster" , dt_mean / dt)
    sol_rs2 = solve(prob_rs, method(), adaptive = false, dt = dt);
    # calculate energy
    Eges_rs2, _, _, _ = calculate_energy_schoemer_reduced(sol_rs2, parameter_rs);
    Eges_rs2 = Eges_rs2 .- first(Eges_rs2) # rescale


    lw = 2
    fontsize = 14
    fontsize2 = 18


    p = plot(sol_rs.t, Eges_rs, label = L"\Delta E"*" Schömer (red.)", xlabel = "t", ylabel = L"\Delta E"*" [a.u.]",   
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

    plot!(sol_rs2.t, Eges_rs2, label = L"\Delta E"*" Schömer konst.", lw=lw)
    n = sol_rs.stats.naccept + sol_rs.stats.nreject
    n2 =length(sol_rs2) - 1
    verhaeltnis = n2/n
    # println("constant method used $verhaeltnis times more steps")
    return p, verhaeltnis
end


tol = 1e-5
method = Tsit5


begin
    tspan = (0.0, 2*34.35) # = (0.0, 2 * sol_k.t[argmin(abs.(sol_k[6,:]))]) to be symmetric 
    ω0 = [0.01, -0.02, -2.0]

    #set up kuypers 
    y0, parameter0 = set_up_original_problem_for_torque(ω0)
    prob_k = ODEProblem(dydt_torque!, y0, tspan, parameter0)
    sol_k = solve(prob_k, method(), reltol=tol, abstol=tol);



    # set up schoemer_red
    y_rs, parameter_rs = setup_rattleback_schoemer_reduced(ω0)
    prob_rs = ODEProblem(rattleback_reduced!, y_rs, tspan, parameter_rs)
    sol_rs = solve(prob_rs, method(), reltol=tol, abstol=tol);

    lw = 2
    fontsize = 14
    fontsize2 = 18



    p1 = plot(sol_k.t[1:end-1], diff(sol_k.t), 
        xlabel = "t", ylabel = "Schrittweiten "*L"\Delta t",   
        label = L"Kuypers",
        # title = "Schrittweiten",
        titlefont = font(fontsize2,"Computer Modern"),
        guidefont = font("Computer Modern"),
        #titlesize = fontsize2,
        dpi = 300, size = (700, 500), 
        lw = lw,
        grid = true,
        legendfontsize=fontsize,
        legend = :bottom,
        xtickfontsize=fontsize,
        ytickfontsize=fontsize,
        ylabelfontsize=fontsize2,
        xlabelfontsize=fontsize2,
        
    
    )
    plot!(sol_rs.t[1:end-1], diff(sol_rs.t), label = L"Schömer \; red.", lw = 2)
    hline!([mean(diff(sol_k.t))], label = "durchschn. Schrittweite "* L"Kuypers", lw = 2, linestyle = :dash)
    hline!([mean(diff(sol_rs.t))], label = "durchschn. Schrittweite "*L"Schömer \; red.", lw = 2, linestyle = :dash)
end



save_at1 = joinpath(@__DIR__, "../plots", "step_sizes.pdf")
@info "Plot looking at step sizes saved at" , savefig(p1, save_at1)







begin

    tols = [1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5]
    tols = reverse(tols)
    x = [plot_stability_and_ratio_with_konstant_schoemer_red(tol)[2] for tol in tols]
    y = [plot_stability_and_ratio_with_konstant_kuypers(tol)[2] for tol in tols]

    p2 = plot(x, xticks = (2:2:length(tols), tols[2:2:end]), 
        xrotation =0,
        label = L"Schömer \; red.",
        xlabel = "Toleranzen "*L"\tau", 
        # ylabel = "# Schritte adaptive / # Schritte konstant", 
        ylabel = "Verhältnis durchschnittlicher \n zu kleinster Schrittweite",
        # title = "Faktor mehr Schritte bei konstanter Schrittweite", 
        titlefont = font(fontsize2,"Computer Modern"),
        guidefont = font("Computer Modern"),
        dpi = 300, size = (700, 500), 
        lw = lw,
        grid = true,
        legendfontsize=fontsize,
        legend = :bottomright,
        xtickfontsize=fontsize,
        ytickfontsize=fontsize,
        ylabelfontsize=fontsize2,
        xlabelfontsize=fontsize2,
    )

    plot!(y, label = L"Kuypers", lw=lw)

end

save_at2 = joinpath(@__DIR__, "../plots", "step_sizes_factor_constant.pdf")
@info "Plot looking at ratio of how much more steps are needed for constant stepsize" , savefig(p2, save_at2)
