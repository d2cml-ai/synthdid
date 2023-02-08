using Statistics
include("data.jl")
function contract3(X, v=Nothing)
  # if(length(size(X)) == 3 && size(X, 3) == length(v)) return "error"
  out = zeros(size(X, 1), size(X, 2))
  if isnothing(v)
    return out
  end
  for ii in eachindex(v)
    out .+= v[ii] * X[:, :, ii]
  end
  return out
end;

# X = reshape(reverse(1:27), (3, 3, 3))
# v = 1:3

# contract3(X, v)

# Frank_wolfe
function fw_step(A, x, b, eta, alpha=nothing)
  Ax = A * x
  half_grad = transpose(Ax .- b) * A + eta * x
  i = findmin(half_grad)[2]
  if !isnothing(alpha)
    x *= (1 - alpha)
    x[i] += alpha
    return x
  else
    d_x = -x
    d_x[i] .= 1 - x[i]
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
end;
# A = reshape(1:9, 3, 3)
# x = [.5, .5, .5]
# b = [1, 2, 3]
# eta = 0.1
# alpha = nothing
# alpha
# isnothing(alpha)
# fw_step(A, x, b, eta, 0.5)

# function sc_weight_fw(Y, zeta, intercept=true, lambda=nothing, min_decrease=1e-3, max_iter=1000)
  T0 = size(Y, 2) - 1
  N0 = size(Y, 1)
  if isnothing(lambda)
    lambda = fill(1 / T0, T0)
  end
  if intercept
    Y = Y.- mean(Y, dims = 1)
  end

  t = 0
  vals = zeros(max_iter)
  # print("A")
  # Y
  A = Y[:, 1:T0]
  b = Y[:, T0+1]
  eta = N0 * real(zeta^2)
  # while t < max_iter && (t < 2 || vals[t-1] - vals[t] > min_decrease^2)
    t += 1
    lambda_p = fw_step(A, lambda, b, nothing)
    lambda = lambda_p
    err = Y[1:N0, :] * [lambda; -1]
    vals[t] = real(zeta^2) * sum(lambda .^ 2) + sum(err .^ 2) / N0
  # end
  # Dict("lambda" => lambda, "vals" => vals)
# end;
# Y
# Y .- mean(Y, dims = 1)
# sc_weight_fw(Y, zeta, intercept, lambda)
Y = reshape(rand(100), 10, 10)
Y = [8 4 1; 10 1 8; 0 6 1]
Y[:, 1]
zeta = 1+1im
intercept = true
lambda = [1, 1, 1]
min_decrease = 1e-3
max_iter = 1000
using Statistics
Y[:, 1:2]

# TODO: sc.weigth.fw.covariates

function sc_weight_fw_covariates(
  Y, x, seta_lambda, zeta_ometa, lambda_intercept, omega_intercept,
  min_decrease=1e-3, max_iter=1000,
  lambda, omega, beta, update_lambda, update_omega
)

end;

struct normal1
  x::Float64
  y::Float64
  function xy(self)
    return self.x + self.y
  end
end

normal1(1.0, 2.9)


struct MyClass4
  x
  y
  suma(self) = self.x + self.y

end

# Crear un objeto de la clase MyClass
# my_obj = MyClass4(1.0, 2.0)



# Llamar el método sum
println(my_obj.sum()) # Imprime 3.0
