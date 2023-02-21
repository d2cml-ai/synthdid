using DataFrames, Plots
using DataFrames, CSV
using Statistics, DataFrames, Distributions
using Statistics, CSV, DataFrames, Shuffle

include("utils.jl")
include("solver.jl")
include("synthdid.jl")
include("data.jl")
include("vcov.jl")
include("plots.jl")
include("summary.jl")

setup_data = panel_matrices(data("california_prop99"))

Y = setup_data.Y
N0 = setup_data.N0
T0 = setup_data.T0

tau_hat = synthdid_estimate(Y, N0, T0);