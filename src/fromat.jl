include("setup.jl")

function format(x::synthdid_est1)
  info = summary(x)
  d = info.dimensions
  summary_result = "synthdid: $(estimate) +- $(1.96 * info.se). Effective N0/N0 = $(d.N0_efective)/$(d.N0)~$(d.N0_efective / d.N0). Effective T0/T0 = $(d.T0_efective)/$(d.T0)~$(d.T0_efective/d.T0). N1,T1 = $(d.N1), $(d.T1)."
  return summary_result
end

