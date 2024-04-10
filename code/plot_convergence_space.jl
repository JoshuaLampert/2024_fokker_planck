# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Interpolations
using CSV
using DataFrames: DataFrame
using TrixiBase: trixi_include
using Plots: gr

gr()

tspan = (0.0, 10.0) # time span
dt = 0.001

errors_euler = []
errors_impeuler = []
errors_mpe = []
errors_mprk = []
errors_heun = []

Ns = 10 * 2 .^ (0:3)
dws = []

refsol = CSV.read("out/refsol.csv", DataFrame)
I = (-1.0, 1.0)
N_ref = size(refsol)[2] - 1
w_interfaces_ref = LinRange(first(I), last(I), N_ref + 1)
dw_ref = w_interfaces_ref[2] - w_interfaces_ref[1]
w_cells_ref = w_interfaces_ref[1:end - 1] .+ dw_ref/2

dt_ref = refsol[2, 1] - refsol[1, 1]
t_ref = refsol[1, 1]:dt_ref:refsol[end, 1]
refsol_interp = cubic_spline_interpolation((t_ref, w_cells_ref), Matrix(refsol[1:(end - 2), 2:end]), extrapolation_bc = Line())
for N in Ns
    trixi_include("fokker_planck.jl", N = N)
    push!(dws, dw)
    # Create Fokker-Planck problem
    prob = ConservativePDSProblem(P!, u0, tspan, p; p_prototype)
    ode_prob = ODEProblem(rhs!, u0, tspan, p)

    sol_euler = solve(ode_prob, Euler(), dt = dt, adaptive = false, save_everystep = true)
    l1_error_euler = dw * sum(abs.(stack(sol_euler.u) .- transpose(refsol_interp(sol_euler.t, w_cells)))) / length(sol_euler.t)
    push!(errors_euler, (l1_error_euler, sol_euler.retcode))

    sol_impeuler = solve(ode_prob, ImplicitEuler(autodiff=false), dt = dt, adaptive = false, save_everystep = true)
    l1_error_impeuler = dw * sum(abs.(stack(sol_impeuler.u) .- transpose(refsol_interp(sol_impeuler.t, w_cells)))) / length(sol_impeuler.t)
    push!(errors_impeuler, (l1_error_impeuler, sol_impeuler.retcode))

    sol_mpe = solve(prob, MPE(), dt = dt, adaptive = false, save_everystep = true)
    l1_error_mpe = dw * sum(abs.(stack(sol_mpe.u) .- transpose(refsol_interp(sol_mpe.t, w_cells)))) / length(sol_mpe.t)
    push!(errors_mpe, (l1_error_mpe, sol_mpe.retcode))

    sol_mprk = solve(prob, MPRK22(1.0), dt = dt, adaptive = false, save_everystep = true)
    l1_error_mprk = dw * sum(abs.(stack(sol_mprk.u) .- transpose(refsol_interp(sol_mprk.t, w_cells)))) / length(sol_mprk.t)
    push!(errors_mprk, (l1_error_mprk, sol_mprk.retcode))

    sol_heun = solve(ode_prob, Heun(), dt = dt, adaptive = false, save_everystep = true)
    l1_error_heun = dw * sum(abs.(stack(sol_heun.u) .- transpose(refsol_interp(sol_heun.t, w_cells)))) / length(sol_heun.t)
    push!(errors_heun, (l1_error_heun, sol_heun.retcode))
end

p = plot(axis = :log, palette = :matter, yrange = (1e-5, 0.1), xlabel = "Δw", ylabel = "average L1-error",
         xtickfontsize = 12, ytickfontsize = 12, xguidefontsize = 18, yguidefontsize = 18, legendfontsize = 10, legend = :bottomright)
plot!(p, dws, first(errors_mprk[2]) * dws.^2 / dws[2]^2, color = :grey, label = "p = 2", linestyle = :dot, linewidth = 3)
scatter!(p, dws, first.(errors_mpe), color = 1, label = "MPE", marker = :heptagon, markerstrokewidth = 0, markersize = 8)
scatter!(p, dws, first.(errors_mprk), color = 65, label = "MPRK", marker = :dtriangle, markerstrokewidth = 0, markersize = 8)
scatter!(p, dws, first.(errors_impeuler), color = 129, label = "Implicit Euler", marker = :star7, markerstrokewidth = 0, markersize = 8)
scatter!(p, dws, first.(errors_euler), color = 180, label = "Euler", marker = :star4, markerstrokewidth = 0, markersize = 8)
scatter!(p, dws, first.(errors_heun), color = 256, label = "Heun", marker = :xcross, markerstrokewidth = 2, markersize = 6)
savefig(p, joinpath(OUT, "order_space.pdf"))

using PrettyTables

orders_str = Vector{String}[]

for errors in [errors_mpe, errors_mprk, errors_impeuler, errors_euler, errors_heun]
    orders_str_ = String[]
    for i in 2:length(errors)
        if last(errors[i]) == ReturnCode.Unstable
            push!(orders_str_, "--")
       else
            push!(orders_str_, string(round(-log2(first(errors[i]) / first(errors[i - 1])), digits = 4)))
       end
    end

    push!(orders_str, orders_str_)
end

pretty_table(stack(orders_str),
             header = ["MPE", "MPRK", "Implicit Euler", "Euler", "Heun"], backend = Val(:latex), tf = tf_latex_booktabs, alignment = :c)
