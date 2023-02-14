using Statistics
include("data.jl")

function contract3(X::Array{Float64,3}, v::Union{Vector,Nothing}=nothing)::Matrix{Float64}
  if !isnothing(v) && size(X, 3) != length(v)
    throw(ArgumentError("The length of `v` must match the size of the third dimension of `X`"))
  end
  out = zeros(Float64, size(X, 1), size(X, 2))
  if isnothing(v)
    return out
  end
  for ii in eachindex(v)
    out .+= v[ii] * X[:, :, ii]
  end
  return out
end
function contract3(X, v::Union{Vector,Nothing}=nothing)::Matrix{Float64}
  if !isnothing(v) && size(X, 3) != length(v)
    throw(ArgumentError("The length of `v` must match the size of the third dimension of `X`"))
  end
  out = zeros(Float64, size(X, 1), size(X, 2))
  if isnothing(v)
    return out
  end
  for ii in eachindex(v)
    out .+= v[ii] * X[:, :, ii]
  end
  return out
end

X = rand(3, 4, 5)
v = rand(5)
# size(X)
rand(3, 4, 12)
contract3(X, nothing)
X

# Frank_wolfe

function fw_step(A::Matrix, x::Vector{Float64}; b::Vector{Float64}, eta::Number, alpha::Union{Nothing,Float64}=nothing)::Vector{Float64}
  Ax = A * x
  half_grad = transpose(Ax .- b) * A + (eta * x)'
  i = findmin(half_grad)[2][2]
  if !isnothing(alpha)
    x *= (1 - alpha)
    x[i] += alpha
    return x
  else
    d_x = -x
    d_x[i] = 1 - x[i]
    if all(d_x .== 0)
      return x
    end
    d_err = A[:, i] - Ax
    step_upper = -half_grad * d_x
    step_bot = sum(d_err .^ 2) + eta * sum(d_x .^ 2)
    step = step_upper[1] / step_bot
    constrained_step = min(1, max(0, step))
    return x + constrained_step * d_x
  end
end

# A

# x
# A = reshape(1:9.0, 3, 3)
# x = [0.5, 0.5, 0.5]
# b = [1.0, 2.0, 3.0]
# eta = 0.1
# alpha = nothing
# alpha
# isnothing(alpha)
# fw_step1(A, x, b, eta)
# A = Matrix(A)
# fw_step1(A, x, b=b, eta=eta, alpha=nothing)

# Y = reshape(rand(100), 10, 10)
# zeta = 1 + 1im
# intercept = true
# lambda = nothing
# min_decrease = 1e-3
# max_iter = 1000
function sc_weight_fw(
  Y::Matrix, zeta::Number;
  intercept::Bool=true, lambda::Union{Vector,Nothing}=nothing,
  min_decrease::Number=1e-3, max_iter::Int64=1000)
  T0 = size(Y, 2) - 1
  N0 = size(Y, 1)
  if isnothing(lambda)
    lambda = fill(1 / T0, T0)
  end
  if intercept
    Y = Y .- mean(Y, dims=1)
  end

  t = 0
  vals = zeros(max_iter)
  # print("A")
  # Y
  A = Y[:, 1:T0]
  b = Y[:, T0+1]
  eta = N0 * real(zeta^2)
  while t < max_iter && (t < 2 || vals[t-1] - vals[t] > min_decrease^2)
    t += 1
    lambda_p = fw_step(A, lambda, b=b, eta=eta)
    lambda = lambda_p
    err = Y[1:N0, :] * [lambda; -1]
    vals[t] = real(zeta^2) * sum(lambda .^ 2) + sum(err .^ 2) / N0
  end
  Dict("lambda" => lambda, "vals" => vals)
end;

# sc_weight_fw(Matrix(reshape(rand(10000), 100, 100)), zeta)
# Y
# Y .- mean(Y, dims = 1)
# sc_weight_fw(Y, zeta, intercept, lambda)
# Y = reshape(rand(100), 10, 10)
# Y = [8 4 1; 10 1 8; 0 6 1]
# Y[:, 1]
# zeta = 1 + 1im
# intercept = true
# lambda = [1, 1, 1]
# min_decrease = 1e-3
# max_iter = 1000
# using Statistics
# Y[:, 1:2]

# TODO: sc.weigth.fw.covariates

# function sc_weight_fw_covariates(
#   Y, x, seta_lambda, zeta_ometa, lambda_intercept, omega_intercept,
#   min_decrease=1e-3, max_iter=1000,
#   lambda, omega, beta, update_lambda, update_omega
# )

# end;


