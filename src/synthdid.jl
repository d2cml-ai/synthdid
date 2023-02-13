include("utils.jl")
include("solver.jl")

# function sparsify_function(v::Vector{Number})
function sparsify_function(v::Vector)
  v[v.<=maximum(v)/4] .= 0
  return v ./ sum(v)
end


# isnothing(sparsify) ? max_iter_pre_sparsify : "sparsify_function"


data0 = random_low_rank()
Y = data0.Y
N0 = data0.n0
T0 = data0.t0

X = zeros(size(Y))

algo = reshape([1, 2, 4, 8, 0, 10, -1, 10, 4, 1, 2, 4, 8, 0, 10, -1, 10, 4, 0, 2], (4, 5))
algo1 = algo[1:3, 1:4]



std(diff(algo1, dims=2))


# function synthdid_estimate(
#   Y::Matrix, N0::Int, T0::Int, X::Matrix=zeros(size(Y));
noise_level::Float64 = std(diff(Y[1:N0, 1:T0], dims=2))#,
eta_omega::Float64 = ((size(Y, 1) - N0) * (size(Y, 2) - T0))^(1 / 4)#,
eta_lambda::Float64 = 1e-6#, 
zeta_omega::Float64 = eta_omega * noise_level#,
zeta_lambda::Float64 = eta_lambda * noise_level#,
omega_intercept::Bool = true#, 
lambda_intercept::Bool = true#,
weights5::Dict{String,Union{Nothing,Vector}} = Dict("omega" => nothing, "lambda" => nothing, "lambda" => [1, 2])
weights5
update_omega::Bool = isnothing(weights5["omega"])#,
update_lambda::Bool = isnothing(weights["lambda"])#,
min_decrease::Float64 = 1e-5 * noise_level
# ,
max_iter::Int = 10000#,
sparsify::Function = sparsify_function#, 
max_iter_pre_sparsify::Int = 100
# )

Dict("lambda" => [12], "omega" => nothing)

# sparsify_function(lambda_opt["lambda"])


N1 = size(Y, 1) - N0
T1 = size(Y, 2) - T0
X = Array{Any,3}(undef, size(Y, 1), size(Y, 2), 0)

if length(size(Y)) == 3
  weights5["vals"] = nothing
  weights5["lambda_vals"] = nothing
  weights5["omega_vals"] = nothing
  if (update_lambda)
    Yc = collapse_form(Y, N0, T0)
    lambda_opt = sc_weight_fw(
      Yc[1:N0, :], zeta_lambda, intercept=lambda_intercept, lambda=weights["lambda"], min_decrease=min_decrease, max_iter=max_iter_pre_sparsify)
    if (!isnothing(sparsify))
      lambda_lambda_opt = sparsify(lambda_opt["lambda"])
      lambda_opt = sc_weight_fw(Yc[1:N0, :], zeta_lambda, intercept=lambda_intercept, lambda=lambda_lambda_opt, min_decrease=min_decrease, max_iter=max_iter)
    end
    weights5["lambda"] = lambda_opt["lambda"]
    weights5["lambda_vals"] = lambda_opt["vals"]
    weights5["vals"] = lambda_opt["vals"]
  end
  if (update_omega)
    Yc = collapse_form(Y, N0, T0)
    matrix_yc = Matrix(transpose(Yc[:, 1:T0]))
    omega_opt = sc_weight_fw(
      matrix_yc, zeta_omega, intercept=omega_intercept, lambda=weights5["omega"], min_decrease=min_decrease, max_iter=max_iter_pre_sparsify)
    if (!isnothing(sparsify))
      omega_lambda_opt = sparsify(omega_opt["lambda"])
      omega_opt = sc_weight_fw(matrix_yc, zeta_omega, intercept=omega_intercept, lambda=omega_lambda_opt, min_decrease=min_decrease, max_iter=max_iter)
    end
    weights5["omega"] = omega_opt["lambda"]
    weights5["omega_vals"] = omega_opt["vals"]
    if (isnothing(weights5["vals"]))
      weights5["vals"] = omega_opt["vals"]
    else
      weights5["vals"] = pairwise_sum_decreasing(weights5["vals"], omega_opt["vals"])
    end
  end
else
  YC = collapse_form(Y, N0, T0)
  Xc = [collapse_form(Xi, N0, T0) for Xi in eachslice(X, dims=3)]
end

X = reshape(reverse(1:27), (3, 3, 3))


Y0 = collapse_form(Y, N0, T0)
# if omega_intercept
#   Y = hcat(ones(size(Y, 1)), Y)
# end
# if lambda_intercept
#   X = hcat(ones(size(X, 1)), X)
# end

# if update_omega
#   weights["omega"] = ones(size(Y, 2))
# end
# if update_lambda
#   print("ls")
# end


# weights["lambda"] = ones(
# end
