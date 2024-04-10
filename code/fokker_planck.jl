# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using PositiveIntegrators
using OrdinaryDiffEq
using QuadGK: quadgk
using SimpleUnPack: @unpack
using LinearAlgebra: Tridiagonal
using Plots

const OUT = "out/"
ispath(OUT) || mkpath(OUT)

# parameters
I = (-1.0, 1.0)
N = 80
w_interfaces = LinRange(first(I), last(I), N + 1) |> collect
dw = w_interfaces[2] - w_interfaces[1]
w_cells = w_interfaces[1:end - 1] .+ dw/2
σ2 = 0.2 # = σ^2
# functions in the model
D = @. σ2/2 * (1 - w_interfaces^2)^2 # diffusion term
Dprime = @. -2 * σ2 * w_interfaces * (1 - w_interfaces^2) # D', i.e. first derivative of diffusion term
# aggregation dynamics operator (use simple quadrature rule to compute integral)
# To avoid allocations, also pass `w_cells` and `dw` to avoid global variables
function B(f_cells, w_cells, w, dw)
    res = 0.0
    for i in eachindex(w_cells, f_cells)
        res += (w - w_cells[i]) * f_cells[i]
    end
    res *= dw
    return res
end

λ = zeros(N) # on interfaces, first entry will not be used, don't need to allocate for last entry
C = zeros(N) # on interfaces, first entry will not be used, don't need to allocate for last entry
δ = zeros(N) # on interfaces, first entry will not be used, don't need to allocate for last entry
BB = zeros(N) # on interfaces, first entry will not be used, don't need to allocate for last entry
u_tilde = zeros(N)
fluxes = zeros(N + 1)
p = (w_interfaces = w_interfaces, w_cells = w_cells, dw = dw, D = D, Dprime = Dprime, B = B, λ = λ, C = C, δ = δ, BB = BB, u_tilde = u_tilde, fluxes = fluxes)

# initial condition
c = 30.0
f0_1(w) = exp(-c*(w + 0.5)^2) + exp(-c*(w - 0.5)^2)
oneoverβ, _ = quadgk(f0_1, first(I), last(I), rtol = 1e-12)
β = 1/oneoverβ
f0(w) = β * f0_1(w)
u0 = f0.(w_cells)

# stationary solution for reference
f0w(w) = w * f0(w)
# call this uu instead of u because u is used as unknown in the `ODEProblem`
uu, _ = quadgk(f0w, first(I), last(I), rtol = 1e-12)
f_stat_1(w) = 1/(1 - w^2)^2 * ((1 + w)/(1 - w))^(uu/(2 * σ2))*exp(-(1 - uu*w)/(σ2 * (1 - w^2)))
oneoverK, _ = quadgk(f_stat_1, first(I), last(I), rtol = 1e-12)
K = 1/oneoverK
f_stat(w) = K * f_stat_1(w)
u_stat = f_stat.(w_cells)


function calculate_values!(BB, λ, C, δ, u_tilde, u, w_cells, w_interfaces, dw, D, Dprime, i)
    BB[i] = B(u, w_cells, w_interfaces[i], dw)
    λ[i] = dw * (BB[i] + Dprime[i])/D[i]
    C[i] = λ[i] * D[i]/dw
    δ[i] = isapprox(λ[i], 0.0) ? 0.5 : 1/λ[i] - 1/expm1(λ[i])
    u_tilde[i] = (1 - δ[i]) * u[i] + δ[i] * u[i - 1]
end

# Classical right hand side (using plain OrdinaryDiffEq.jl, not PositiveIntegrators.jl)
# We use `rhs!` for non-Patanker schemes and `P!` for the Patanker methods (both describe the same ODE)
# We could also use `P!` for the non-Patanker schemes, but `rhs!` is more efficient
function rhs!(du, u, p, t)
    @unpack w_interfaces, w_cells, dw, D, Dprime, B, λ, C, δ, BB, u_tilde, fluxes = p
    N = length(w_interfaces) - 1
    # no flux boundary conditions
    fluxes[1] = 0.0
    fluxes[end] = 0.0
    for i in 2:N
        calculate_values!(BB, λ, C, δ, u_tilde, u, w_cells, w_interfaces, dw, D, Dprime, i)
        fluxes[i] = C[i] * u_tilde[i] + D[i] * (u[i] - u[i - 1])/dw
    end
    for i in 1:N
        du[i] = 1/dw * (fluxes[i + 1] - fluxes[i])
    end
end

# In-place implementation of the P matrix for the Fokker-Planck model
function P!(P, u, p, t)
    @unpack w_interfaces, w_cells, dw, D, Dprime, B, λ, C, δ, BB, u_tilde = p
    N = length(w_interfaces) - 1
    for i = 2:N
        calculate_values!(BB, λ, C, δ, u_tilde, u, w_cells, w_interfaces, dw, D, Dprime, i)
    end
    for i = 2:(N - 1)
        P[i, i + 1] = (max(0, C[i + 1]) * u_tilde[i + 1] + D[i + 1] * u[i + 1]/dw)/dw
        P[i, i - 1] = (-min(0, C[i]) * u_tilde[i] + D[i] * u[i - 1]/dw)/dw
    end
    P[1, 2] = (max(0, C[2]) * u_tilde[2] + D[2] * u[2]/dw)/dw
    P[N, N - 1] = (-min(0, C[N]) * u_tilde[N] + D[N] * u[N - 1]/dw)/dw
    return nothing
end

p_prototype = Tridiagonal(zeros(length(u0) - 1),
                          zeros(length(u0)),
                          zeros(length(u0) - 1))
