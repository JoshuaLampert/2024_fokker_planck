# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Interpolations
using CSV
using DataFrames: DataFrame
using TrixiBase: trixi_include

trixi_include("fokker_planck.jl", N = 160)
gr()

tspan = (0.0, 10.0) # time span
# Create Fokker-Planck problem
prob = ConservativePDSProblem(P!, u0, tspan, p; p_prototype)
ode_prob = ODEProblem(rhs!, u0, tspan, p)
dt_ref = dw^2/(2*σ2)
sol_ref = solve(ode_prob, Heun(), dt = dt_ref, adaptive = false, save_everystep = true)

errors_impeuler = []
errors_mpe = []
errors_mprk = []
dts = 0.5 .^ (3:10)

for dt in dts
    sol_impeuler = solve(ode_prob, ImplicitEuler(autodiff=false), dt = dt, adaptive = false, save_everystep = true)
    l1_error_impeuler = dw * sum(abs.(stack(sol_impeuler.u) .- sol_ref(sol_impeuler.t))) / length(sol_impeuler.t)
    push!(errors_impeuler, l1_error_impeuler)

    sol_mpe = solve(prob, MPE(), dt = dt, adaptive = false, save_everystep = true)
    l1_error_mpe = dw * sum(abs.(stack(sol_mpe.u) .- sol_ref(sol_mpe.t))) / length(sol_mpe.t)
    push!(errors_mpe, l1_error_mpe)

    sol_mprk = solve(prob, MPRK22(1.0), dt = dt, adaptive = false, save_everystep = true)
    l1_error_mprk = dw * sum(abs.(stack(sol_mprk.u) .- sol_ref(sol_mprk.t))) / length(sol_mprk.t)
    push!(errors_mprk, l1_error_mprk)
end

p = plot(axis = :log, palette = :matter, xlabel = "Δt", ylabel = "average L1-error",
         xtickfontsize = 12, ytickfontsize = 12, xguidefontsize = 18, yguidefontsize = 18, legendfontsize = 10, legend = :bottomright)
plot!(p, dts, errors_mpe[4] * dts / dts[4], color = :grey, label = "p = 1", linestyle = :dot, linewidth = 2)
plot!(p, dts, errors_mprk[4] * dts.^2 / dts[4]^2, color = :grey, label = "p = 2", linestyle = :dot, linewidth = 2)
scatter!(p, dts, errors_mpe, color = 1, label = "MPE", marker = :heptagon, markerstrokewidth = 0, markersize = 8)
scatter!(p, dts, errors_mprk, color = 65, label = "MPRK", marker = :dtriangle, markerstrokewidth = 0, markersize = 8)
scatter!(p, dts, errors_impeuler, color = 129, label = "Implicit Euler", marker = :star7, markerstrokewidth = 0, markersize = 8)
savefig(p, joinpath(OUT, "order_time.pdf"))

using PrettyTables

orders_str = Vector{String}[]

for errors in [errors_mpe, errors_mprk, errors_impeuler]
    orders_str_ = String[]
    for i in 2:length(errors)
        push!(orders_str_, string(round(-log2(errors[i] / errors[i - 1]), digits = 4)))
    end
    push!(orders_str, orders_str_)
end

pretty_table(stack(orders_str),
             header = ["MPE", "MPRK", "Implicit Euler"], backend = Val(:latex), tf = tf_latex_booktabs, alignment = :c)
