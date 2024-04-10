# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using BenchmarkTools

include("fokker_planck.jl")
gr()

tspan = (0.0, 10.0) # time span
# Create Fokker-Planck problem
prob = ConservativePDSProblem(P!, u0, tspan, p; p_prototype)
ode_prob = ODEProblem(rhs!, u0, tspan, p)

error_times_euler = []
error_times_impeuler = []
error_times_mpe = []
error_times_mprk = []
error_times_heun = []
dts = 0.7.^(0:18)

BenchmarkTools.DEFAULT_PARAMETERS.samples = 5
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60
for dt in dts
    bench_euler = @benchmark (global sol_euler = solve($ode_prob, Euler(), dt = $dt, adaptive = false, save_everystep = false)) evals = 1
    # error is the same for all runs, so just take the result from the last run
    l1_error_euler = dw * sum(abs.(sol_euler.u[end] .- u_stat))
    push!(error_times_euler, (l1_error_euler, bench_euler.times))

    bench_impeuler = @benchmark (global sol_impeuler = solve($ode_prob, ImplicitEuler(autodiff=false), dt = $dt, adaptive = false, save_everystep = false)) evals = 1
    l1_error_impeuler = dw * sum(abs.(sol_impeuler.u[end] .- u_stat))
    push!(error_times_impeuler, (l1_error_impeuler, bench_impeuler.times))

    bench_mpe = @benchmark (global sol_mpe = solve($prob, MPE(), dt = $dt, adaptive = false, save_everystep = false)) evals = 1
    l1_error_mpe = dw * sum(abs.(sol_mpe.u[end] .- u_stat))
    push!(error_times_mpe, (l1_error_mpe, bench_mpe.times))

    bench_mprk = @benchmark (global sol_mprk = solve($prob, MPRK22(1.0), dt = $dt, adaptive = false, save_everystep = false)) evals = 1
    l1_error_mprk = dw * sum(abs.(sol_mprk.u[end] .- u_stat))
    push!(error_times_mprk, (l1_error_mprk, bench_mprk.times))

    bench_heun = @benchmark (global sol_heun = solve($ode_prob, Heun(), dt = $dt, adaptive = false, save_everystep = false)) evals = 1
    l1_error_heun = dw * sum(abs.(sol_heun.u[end] .- u_stat))
    push!(error_times_heun, (l1_error_heun, bench_heun.times))
end

median_from_times(error_times) = median.(last.(error_times)) ./ 1e9 # in s

median_mpe = median_from_times(error_times_mpe)
median_mprk = median_from_times(error_times_mprk)
median_impeuler = median_from_times(error_times_impeuler)
median_euler = median_from_times(error_times_euler)
median_heun = median_from_times(error_times_heun)

p = plot(axis = :log, palette = :matter, xrange = (4e-5, 5.0), yrange = (2.55e-4, 2.8e-4), xlabel = "computation time [s]", ylabel = "L1-error at last timestep",
         xtickfontsize = 12, ytickfontsize = 12, xguidefontsize = 18, yguidefontsize = 18, legendfontsize = 10, legend = :topright,
         yticks = ([2.6e-4, 2.7e-4], ["2.6⋅10⁻⁴", "2.7⋅10⁻⁴"]), bottom_margin = 0.3 * Plots.cm)
scatter!(p, median_mpe, first.(error_times_mpe), color = 1, label = "MPE", marker = :heptagon, markerstrokewidth = 0, markersize = 8)
scatter!(p, median_mprk, first.(error_times_mprk), color = 65, label = "MPRK", marker = :dtriangle, markerstrokewidth = 0, markersize = 8)
scatter!(p, median_impeuler, first.(error_times_impeuler), color = 129, label = "Implicit Euler", marker = :star7, markerstrokewidth = 0, markersize = 8)
scatter!(p, median_euler, first.(error_times_euler), color = 180, label = "Euler", marker = :star4, markerstrokewidth = 0, markersize = 8)
scatter!(p, median_heun, first.(error_times_heun), color = 256, label = "Heun", marker = :xcross, markerstrokewidth = 2, markersize = 6)


savefig(p, joinpath(OUT, "error_last_vs_computation_time.pdf"))
