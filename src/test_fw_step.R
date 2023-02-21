
#install.packages("synthdid")

#install.packages("remotes")
#remotes::install_github("synth-inference/synthdid")


library(synthdid)

contract3 = function(X, v) {
  stopifnot(length(dim(X)) == 3, dim(X)[3] == length(v))
  out = array(0, dim = dim(X)[1:2])
  if (length(v) == 0) { return(out) }
  for (ii in 1:length(v)) {
    out = out + v[ii] * X[, , ii]
  }
  return(out)
}


fw.step = function(A, x, b, eta, alpha = NULL) {
  Ax = A %*% x
  half.grad = t(Ax - b) %*% A + eta * x
  i = which.min(half.grad)
  if (!is.null(alpha)) {
    x = x * (1 - alpha)
    x[i] = x[i] + alpha
    return(x)
  } else {
    d.x = -x; d.x[i] = 1 - x[i]
    if (all(d.x == 0)) { return(x) }
    d.err = A[, i] - Ax
    step = -t(c(half.grad)) %*% d.x / (sum(d.err^2) + eta * sum(d.x^2))
    constrained.step = min(1, max(0, step))
    return(x + constrained.step * d.x)
  }
}


sc.weight.fw = function(Y, zeta, intercept = TRUE, lambda = NULL, min.decrease = 1e-3, max.iter = 1000) {
  T0 = ncol(Y) - 1
  N0 = nrow(Y)
  if (is.null(lambda)) { lambda = rep(1 / T0, T0) }
  if (intercept) {
    Y = apply(Y, 2, function(col) { col - mean(col) })
  }
  
  t = 0
  vals = rep(NA, max.iter)
  A = Y[, 1:T0]
  b = Y[, T0 + 1]
  eta = N0 * Re(zeta^2)
  while (t < max.iter && (t < 2 || vals[t - 1] - vals[t] > min.decrease^2)) {
    t = t + 1
    lambda.p = fw.step(A, lambda, b, eta)
    lambda = lambda.p
    err = Y[1:N0, ] %*% c(lambda, -1)
    vals[t] = Re(zeta^2) * sum(lambda^2) + sum(err^2) / N0
  }
  list(lambda = lambda, vals = vals)
}




california_prop99
setup = panel.matrices(california_prop99)
t(setup$Y)
Y = t(setup$Y)

#Y = read.csv(file = 'C:/Users/User/Documents/GitHub/synthdid/data/CPS.csv', sep = ";")
#head(Y)
#Y = Y[,c(3:5)]

dim(Y)
X = array(0, dim = c(dim(Y), 1))
zeta.lambda = 0
zeta.omega = 0
lambda.intercept = TRUE
omega.intercept = TRUE
min.decrease = 1e-3
max.iter = 1000
lambda = NULL
omega = NULL
beta = NULL
update.lambda = TRUE
update.omega = TRUE

#stopifnot(length(dim(Y)) == 2, length(dim(X)) == 3, all(dim(Y) == dim(X)[1:2]), all(is.finite(Y)), all(is.finite(X)))

T0 = ncol(Y) - 1
N0 = nrow(Y) - 1
if (length(dim(X)) == 2) { dim(X) = c(dim(X), 1) }

if (is.null(lambda)) {  lambda = rep(1 / T0, T0)   }
if (is.null(omega)) {  omega = rep(1 / N0, N0)    }
if (is.null(beta)) {  beta = rep(0, dim(X)[3]) }

Y.lambda = if (lambda.intercept) { apply(Y[1:N0, ], 2, function(row) { row - mean(row) }) } else { Y[1:N0, ] }

fw.step(Y.lambda[, 1:T0], lambda, Y.lambda[, T0 + 1], N0 * Re(zeta.lambda^2), 0.8)



#A = reshape(1:9.0, 3, 3)
x = c(0.5, 0.5, 0.5)
b = c(1.0, 2.0, 3.0)
eta = 0.1
alpha = NULL
alpha
is.null(alpha)
# fw_step1(A, x, b, eta)
A = matrix(1:9, ncol = 3, nrow = 3, byrow = FALSE)
fw.step(A, x, b, eta, alpha= NULL)




b = sc.weight.fw(Y, 4, lambda = c(1:38))

b$lambda
b$vals


























