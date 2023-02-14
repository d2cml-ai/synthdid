function synthdid_plot(estimates::synthdid_est1, treated_name="treated", control_name="synthetic control",
  spaghetti_units=nothing, spaghetti_matrices=nothing,
  facet=nothing, facet_vertical=true, lambda_comparable=!isnothing(facet), overlay=0,
  lambda_plot_scale=3, trajectory_linetype=1, effect_curvature=3, line_width=5, guide_linetype=2, point_size=1,
  trajectory_alpha=5, diagram_alpha=95, effect_alpha=95, onset_alpha=3, ci_alpha=3,
  spaghetti_line_width=2, spaghetti_label_size=2,
  spaghetti_line_alpha=3, spaghetti_label_alpha=5,
  se_method="jackknife", alpha_multiplier=nothing)
end


function synthdid_placebo_plot(estimate::synthdid_est1)

end

function synthdid_placebo_plot(estimate::synthdid_est1)

end

function synthdid_rmse_plot(estimate)

end