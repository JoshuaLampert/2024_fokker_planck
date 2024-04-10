# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DataFrames: DataFrame
using CSV
using TrixiBase: trixi_include

trixi_include("fokker_planck.jl", N = 640)

tspan = (0.0, 10.0) # time span
dt = dw^2/(2*σ2)
# Create Fokker-Planck problem
ode_prob = ODEProblem(rhs!, u0, tspan, p)

@time sol = solve(ode_prob, Euler(), dt = dt, adaptive = false, save_everystep = true)
CSV.write(joinpath(OUT, "refsol.csv"), DataFrame(sol))
