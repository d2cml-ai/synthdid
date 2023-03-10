using Statistics, Convex
function simplex_least_squares(; A, b, zeta=0, intercept=false)
  # n, t = size(A)
  # x = Variable(t)
  # constrainst = [sum(x) == 1, x .>= 0]
  # algo = zeta^2 * length(b) * sum(x .^ 2)
  # if intercept
  #   x0 = Variable(0)
  #   objetive = sum((A * x + x0 - b) .^ 2) + algo
  # else
  #   objetive = sum((A * x + B) .^ 2) + algo
  # end
  # problem = minimize(objetive, constrainst)
  # solve!(problem, ECOS.Optimizer)
  # return x.value
end




function sigma_defaul(Y, N0, T0)
  result = std(diff(Y[1:N0, 1:T0], dims=1), dims=1)
  return result
end

# Y = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20; 21 22 23 24 25]
# N0 = 3
# T0 = 3

# sigma_defaul(Y, N0, T0)

epsilon = 1e-6

function synthdid_reference(; Y, N0, T0, zeta_omega=((size(Y, 1) - N0) * (size(Y, 2) - T0)^(1 / 4) * sigma_defaul(Y, N0, T0)))
  N, T = size(Y)
  n1 = N - N0
  t1 = T - T0
  y_target = Y[1:n0, 1:T0]
  lambda = simplex_least_squares(y_target, mean(Y[1:N0, T0+1:T], dims=2), zeta=epsilon * sigma_default(Y, N0, T0), intercept=true)
  omega = simplex_least_squares(y_target', mean(Y[N0+1:N, 1:T0], dims=1), zeta=zeta_omega, intercept=true)
  estimate = transpose(hcat(-omega, fill(1 / N1, N1))) * Y * hcat(-lambda, fill(1 / T1, T1))
  return estimate
end

function sc_reference(Y, N0, T0, zeta_omega=epsilon * sigma_defaul(Y, N0, T0))
  n, t = size(Y)
  n1 = n - n0
  t1 = t - t0
  y_target = transpose(Y[1:N0, 1:T0], mean(Y[(N0+1):N, 1:T0], dims=1))
  omega = simplex_least_squares(y_target, zeta=zeta_omega, intercept=false)
  estimate = transpose(hcat(-omega, fill(1 / n1, n1))) * Y * hcat(fill(0, t0), fill(1 / t1, t1))
  return estimate

end

function did_reference(Y, N0, T0)
  n, t = size(Y)
  n1, t1 = n - N0, t - T0
  estimate = transpose(vcat(fill(-1 / N0, N0), fill(1 / n1, n1))) * Y * vcat(fill(1 / T0, T0), fill(1 / t1, t1))
  return estimate
end

Y = [1 3; 2 4]
did_reference(Y, 1, 1)

n, t = size(Y)

n1, t1 = n - 1, t - 1

vcat(fill(-1 / 1, 2), fill(1 / n1, n1))