# download(
#   "http://www.damianclarke.net/stata/quota_example.dta",
#   "./vignettes/data/quota.dta"
# )

pwd()
include("../src/setup.jl")
using Random

Random.seed!(1234)


cal = data("california_prop99")
quote_example = data("quota_example")
setup = panel_matrices(quote_example)

# setup = panel_matrices(cal);
Y, N0, T0 = setup.Y, setup.N0, setup.T0

tau_hat = synthdid_estimate(Y, N0, T0)

sc_estimate(Y, N0, T0)

point = tau_hat.estimate

p = synthdid_plot(tau_hat)
p["plot"]

vcov_synthdid_estimate(tau_hat, method="jackknife")


# The parallelogram

## The treatment Effect

tau_hat_plot = synthdid_plot(tau_hat, se_method="placebo", year_unit_trayectory=setup.time)
plot(tau_hat_plot["plot"])

synthdid_units_plot(tau_hat, x_ticks=setup.names, se_method="placebo", negligible_alpha=1.0, negligible_threshold=0.001)

tau_hat_overlay = synthdid_plot(tau_hat, overlay=1, se_method="placebo", year_unit_trayectory=setup.time)
plot(tau_hat_overlay["plot"])


# Checking for pre-treatment parallel trend
tau_hat_overlay = synthdid_plot(tau_hat, overlay=0.8, se_method="placebo", year_unit_trayectory=setup.time)
plot(tau_hat_overlay["plot"])


# compare to other estimators

tau_sc = sc_estimate(Y, N0, T0);
tau_did = did_estimate(Y, N0, T0);

Dict(
  "2 Diff-in-Diff" => tau_did.estimate,
  "2 Synthetic control" => tau_sc.estimate,
  "1 Synthetic Diff-in-Diff" => tau_hat.estimate
)

p1 = synthdid_plot(tau_did, year_unit_trayectory=setup.time)["plot"]
p2 = synthdid_plot(tau_sc, year_unit_trayectory=setup.time)["plot"]
p3 = synthdid_plot(tau_hat, year_unit_trayectory=setup.time)["plot"]
title!(p1, "Diff-in-Diff")
title!(p2, "Synthetic control")
title!(p3, "Synthetic Diff-in-Diff")



plot(
  p1, p2, p3
)

plot!(legend=false)


p1 = synthdid_units_plot(tau_did, x_ticks=setup.names, se_method="placebo")
ylims!(-120, 40)
xticks!([0])
p2 = synthdid_units_plot(tau_sc, x_ticks=setup.names, se_method="placebo")
ylims!(-120, 40)
xticks!([0])
p3 = synthdid_units_plot(tau_hat, x_ticks=setup.names, se_method="placebo")
ylims!(-120, 40)
xticks!([0])

title!(p1, "Diff-in-Diff")
title!(p2, "Synthetic control")
title!(p3, "Synthetic Diff-in-Diff")
plot(p1, p2, p3)


p0 = plot(p1, p2, p3, layout=(1, 3))
plot(p0[3])


