# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

include("fokker_planck.jl")
gr()

tspan = (0.0, 10.0) # time span
dt = dw^2/(2*σ2)
# Create Fokker-Planck problem
prob = ConservativePDSProblem(P!, u0, tspan, p; p_prototype)

saveat = range(tspan..., length = 70)
#method = MPE()
method = MPRK22(1.0)
#method = ImplicitEuler(autodiff = false)
#method = Euler()
#method = Heun()
@time sol = solve(prob, method, dt = dt, adaptive = false, save_everystep = false, saveat = saveat)

method_name = nameof(typeof(method))
if method_name == :MPRK22
   method_name = :MPRK
end

p = plot3d(xlabel = "w", ylabel = "t", zlabel = "f(w, t)", legend = :top, camera = (45, 20), palette = :matter)
plot3d!(p, w_cells, sol.t[2] * ones(N), sol[2], color = :grey, label = "numerical solution ($method_name)")
for i in 3:length(sol.t)
   plot3d!(p, w_cells, sol.t[i] * ones(N), sol[i], color = :grey, label = :none)
end
plot3d!(p, w_cells, sol.t[1] * ones(N), sol[1], color = 100, label = "initial condition", linestyle = :dash)
scatter3d!(p, w_cells, sol.t[end] * ones(N), u_stat, color = 200, label = "stationary solution", markersize = 1.5, markerstrokewidth = 0)

savefig(joinpath(OUT, "solution_$method_name.pdf"))
