function contract3(X, v = Nothing)
  if(length(size(X)) == 3 && size(X, 3) == length(v)) return "error"
  out = zeros(size(X, 1), size(X, 2))
  if isnothing(v)
    return out
  end
  for ii in eachindex(v)
    out .+= v[ii] * X[:, :, ii]
  end
  return out
end;

# Frank_wolfe
function fw_step(A, x, b, eta, alpha=Nothing)
  Ax = A * x
  half_grad = transpose(Ax .- b) * A + eta * x
  i = findmin(half_grad)[2]
  if !isnothing(alpha)
    x *= (1 - alpha)
    x[i] += alpha
    return x
  else
    d_x = x.copy() .* -1
    d_x[i] = 1 - x[i]
    if all(d_x .== 0)
      return x
    end
    d_err = A[:, i] - Ax
    step = -transpose(half_grad) * d_x / sum(d_err .^ 2) + eta * sum(d_x .^ 2)
    constrained_step = min(1, max(0, step))
    return x + constrained_step * d_x
  end
end;


function frank_wolfe(Y, zeta, max_iter, min_decrease, intercept, lambda)
  T0 = size(Y, 2) - 1
  N0 = size(Y, 1)
  if isnothing(lambda)
    lambda = fill(1 / T0, T0)
  end
  if intercept
    Y = [Y[:, i] .- mean(Y[:, i]) for i in 1:T0]
  end

  t = 0
  vals = zeros(max_iter)
  A = Y[:, 1:T0]
  b = Y[:, T0+1]
  eta = N0 * real(zeta^2)
  while t < max_iter && (t < 2 || vals[t-1] - vals[t] > min_decrease^2)
    t += 1
    lambda_p = fw_step(A, lambda, b, eta)
    lambda = lambda_p
    err = Y[1:N0, :] * [lambda; -1]
    vals[t] = real(zeta^2) * sum(lambda .^ 2) + sum(err .^ 2) / N0
  end
  Dict("lambda" => lambda, "vals" => vals)
end;

# TODO: sc.weigth.fw.covariates

function sc_weight_fw_covariates(
    Y, x, seta_lambda, zeta_ometa, lambda_intercept, omega_intercept,
    min_decrease=1e-3, max_iter=1000,
    lambda, omega, beta, update_lambda, update_omega
  )

end;