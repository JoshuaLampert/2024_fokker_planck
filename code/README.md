# Numerical Experiments

This directory contains all code required to reproduce the numerical
experiments. First, you need to install Julia, e.g., by downloading
the binaries from the [download page](https://julialang.org/downloads/).
The numerical experiments were performed using Julia v1.10.1.

The experiments are structured in different scripts creating the figures
and tables displayed in the article. The file `fokker_planck.jl` defines
the necessary structures and semidiscretisation, which is used in the other
files. The following list describes which script creates which figure(s) or
table and the names of the resulting .pdf files:

* Figure 2: `plot_solution_3d.jl` -> `solution_MPRK.pdf`
* Figures 3, 4, 5: `plot_error_vs_stationary.jl` -> `error_vs_stationary_0.0016.pdf`, `error_vs_stationary_0.025.pdf`, `error_vs_stationary_0.25.pdf`
* Figure 6: `plot_convergene_space.jl` -> `order_space.pdf`
* Figure 7: `plot_convergene_time.jl` -> `order_time.pdf`
* Figure 8: `plot_error_average_vs_computation_time.jl` -> `error_average_vs_computation_time.pdf`
* Figure 9: `plot_error_last_vs_computation_time.jl` -> `error_last_vs_computation_time.pdf`
* Table 1: `table_computation_time.jl`

The resulting figures are then saved as .pdf files in a new directory `out`
inside the folder of this `README.md`. The table (Table 1) is printed to the screen
as LaTeX code. Note that the Figures 8 and 9 as well as the Table 1 display computation times,
i.e. the results can differ depending on your machine. Note that Figures 6, 7, 8, and 9 need
a reference solution that has to be computed first. To do so, first execute the script
`save_ref_sol.jl`.

In order to execute a script, start Julia in this folder and execute

```julia
julia> include("path_to/file_name.jl")
```

in the Julia REPL. To execute the first script from the list above, e.g.,
execute

```julia
julia> include("plot_solution_3d.jl")
```
