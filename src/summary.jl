function synthdid_controls(estimates, sort_by=1, mass=0.9, weight_type="omega")
  # if (class(estimates) == 'synthdid_estimate') { estimates = list(estimates) }
  # if (is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
  # if (!weight.type %in% c('omega', 'lambda')) { stop('weight.type must be "omega" or "lambda"') } 
  # weights = do.call(cbind, lapply(estimates, function(est) { attr(est, 'weights')[[weight.type]] }))
  # if (is.null(dim(weights))) { dim(weights) = c(length(weights), 1) }

  # Y = attr(estimates[[1]], 'setup')$Y
  # o = rev(order(weights[, sort.by]))
  # tab = weights[o, , drop = FALSE]
  # rownames(tab) = if(weight.type == 'omega') { rownames(Y)[o] } else { colnames(Y)[o] }
  # colnames(tab) = names(estimates)

  # # truncate table to retain a weight sum of at least mass for each unit
  # tab.len = max(apply(tab, 2, function(col) { Position(function(x) { x >= mass }, cumsum(col), nomatch=nrow(tab)) }))
  # tab[1:tab.len, , drop=FALSE]
end

# summary.synthdid_estimate = function(object, weight.digits=3, fast=FALSE, ...) {
#   N0 = attr(object, 'setup')$N0
#   T0 = attr(object, 'setup')$T0
#   list(estimate = c(object),
#     se = sqrt(if(fast) { vcov(object, method = 'jackknife') } else { vcov(object) }),
#     controls = round(synthdid_controls(object, weight.type='omega'),  digits=weight.digits),
#     periods  = round(synthdid_controls(object, weight.type='lambda'), digits=weight.digits),
#     dimensions = c( N1 = nrow(Y(object))-N0, N0 = N0, N0.effective = round(1 / sum(omega(object)^2),  weight.digits),
# 		    T1 = ncol(Y(object))-T0, T0 = T0, T0.effective = round(1 / sum(lambda(object)^2), weight.digits)))
# }