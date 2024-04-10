# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using LaTeXStrings
using PrettyTables
using BenchmarkTools

include("fokker_planck.jl")

tspan = (0.0, 10.0) # time span
# Create Fokker-Planck problem
prob = ConservativePDSProblem(P!, u0, tspan, p; p_prototype)
ode_prob = ODEProblem(rhs!, u0, tspan, p)

times_euler = []
times_impeuler = []
times_mpe = []
times_mprk = []
times_heun = []
dts = [dw^2.5/(2*σ2), dw^2/(2*σ2), dw, dw/(2*σ2), 10*dw]

BenchmarkTools.DEFAULT_PARAMETERS.samples = 5
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 300
for dt in dts
    bench_euler = @benchmark (global sol = solve($ode_prob, Euler(), dt = $dt, adaptive = false, save_everystep = false)) evals = 1
    push!(times_euler, (bench_euler.times, sol.retcode))

    bench_impeuler = @benchmark (global sol = solve($ode_prob, ImplicitEuler(autodiff=false), dt = $dt, adaptive = false, save_everystep = false)) evals = 1
    push!(times_impeuler, (bench_impeuler.times, sol.retcode))

    bench_mpe = @benchmark (global sol = solve($prob, MPE(), dt = $dt, adaptive = false, save_everystep = false)) evals = 1
    push!(times_mpe, (bench_mpe.times, sol.retcode))

    bench_mprk = @benchmark (global sol = solve($prob, MPRK22(1.0), dt = $dt, adaptive = false, save_everystep = false)) evals = 1
    push!(times_mprk, (bench_mprk.times, sol.retcode))

    bench_heun = @benchmark (global sol = solve($ode_prob, Heun(), dt = $dt, adaptive = false, save_everystep = false)) evals = 1
    push!(times_heun, (bench_heun.times, sol.retcode))
end

# To compute the first non-zero digit of a float
ilog10(x) = Int(floor(log10(x)))
alog10(x) = abs(ilog10(x))
times_str = Vector{AbstractString}[]

for times in [times_mpe, times_mprk, times_impeuler, times_euler, times_heun]
    times_str_ = AbstractString[]
    for (time, retcode) in times
       if retcode == ReturnCode.Unstable
            push!(times_str_, "--")
       else
            mean_ms = mean(time) / 10^6
            std_ms = std(time) / 10^6
            n_digits = max(2, alog10(std_ms))
            std_ms = round(std_ms, digits = n_digits)
            if mean_ms < 1000.0
                push!(times_str_, latexstring("(", round(mean_ms, digits = 2), " ", L"\pm", " ", std_ms, ")~ms"))
            else
                push!(times_str_, latexstring(round(mean_ms / 1000, digits = 2), "~s ", L"\pm", " ", std_ms, "~ms"))
            end
       end
    end
    push!(times_str, times_str_)
end

pretty_table([[L"\frac{\Delta_w^{2.5}}{2\sigma^2}", L"\frac{\Delta_w^2}{2\sigma^2}", L"\Delta_w", L"\frac{\Delta_w}{2\sigma^2}", L"10\Delta_w"] times_str...],
             header = [L"\Delta_t", "MPE", "MPRK", "Implicit Euler", "Euler", "Heun"], backend = Val(:latex), tf = tf_latex_booktabs, alignment = :c)
