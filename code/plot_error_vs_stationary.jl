# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

include("fokker_planck.jl")
pyplot()

tspan = (0.0, 10.0) # time span
# Create Fokker-Planck problem
prob = ConservativePDSProblem(P!, u0, tspan, p; p_prototype)
ode_prob = ODEProblem(rhs!, u0, tspan, p)
saveat = range(tspan..., length = 100)

dts = [dw^2/(2*σ2), dw, 10 * dw]

for dt in dts
    @time sol_euler = solve(ode_prob, Euler(), dt = dt, adaptive = false, save_everystep = false, saveat = saveat)
    @time sol_impeuler = solve(ode_prob, ImplicitEuler(autodiff = false), dt = dt, adaptive = false, save_everystep = false, saveat = saveat)
    @time sol_mpe = solve(prob, MPE(), dt = dt, adaptive = false, save_everystep = false, saveat = saveat)
    @time sol_mprk = solve(prob, MPRK22(1.0), dt = dt, adaptive = false, save_everystep = false, saveat = saveat)
    @time sol_heun = solve(ode_prob, Heun(), dt = dt, adaptive = false, save_everystep = false, saveat = saveat)

    l1_error_euler = []
    for i in 1:length(sol_euler.t)
        push!(l1_error_euler, sum(abs.(sol_euler.u[i] .- u_stat)))
    end
    l1_error_impeuler = []
    for i in 1:length(sol_impeuler.t)
        push!(l1_error_impeuler, sum(abs.(sol_impeuler.u[i] .- u_stat)))
    end
    l1_error_mpe = []
    for i in 1:length(sol_mpe.t)
        push!(l1_error_mpe, sum(abs.(sol_mpe.u[i] .- u_stat)))
    end
    l1_error_mprk = []
    for i in 1:length(sol_mprk.t)
        push!(l1_error_mprk, sum(abs.(sol_mprk.u[i] .- u_stat)))
    end
    l1_error_heun = []
    for i in 1:length(sol_heun.t)
        push!(l1_error_heun, sum(abs.(sol_heun.u[i] .- u_stat)))
    end

    ylim = (1e-4, 10.0)
    plot(palette = :matter,xlabel = "t", ylabel = "L1-error", xtickfontsize = 12, ytickfontsize = 12, xguidefontsize = 18, yguidefontsize = 18, legendfontsize = 10)
    plot!(sol_mpe.t, dw .* l1_error_mpe, yaxis = :log, label = "MPE", ylim = ylim, color = 1, linestyle = :solid, linewidth = 2)
    plot!(sol_mprk.t, dw .* l1_error_mprk, yaxis = :log, label = "MPRK", ylim = ylim, color = 65, linestyle = :dash, linewidth = 2)
    plot!(sol_impeuler.t, dw .* l1_error_impeuler, yaxis = :log, label = "Implicit Euler", color = 129, linestyle = :dashdot, linewidth = 2)
    plot!(sol_euler.t, dw .* l1_error_euler, yaxis = :log, label = "Euler", ylim = ylim, color = 180, linestyle = :dot, linewidth = 2)
    plot!(sol_heun.t, dw .* l1_error_heun, yaxis = :log, label = "Heun", ylim = ylim, color = 256, linestyle = :dot, linewidth = 2)

    dt_rounded = round(dt, digits = 4)
    savefig(joinpath(OUT, "error_vs_stationary_$dt_rounded.pdf"))
end

