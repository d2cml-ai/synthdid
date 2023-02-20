function synthdid_plot(estimates::synthdid_est1, treated_name="treated", control_name="synthetic control",
  spaghetti_units=nothing, spaghetti_matrices=nothing,
  facet=nothing, facet_vertical=true, lambda_comparable=!isnothing(facet), overlay=0,
  lambda_plot_scale=3, trajectory_linetype=1, effect_curvature=3, line_width=5, guide_linetype=2, point_size=1,
  trajectory_alpha=5, diagram_alpha=95, effect_alpha=95, onset_alpha=3, ci_alpha=3,
  spaghetti_line_width=2, spaghetti_label_size=2,
  spaghetti_line_alpha=3, spaghetti_label_alpha=5,
  se_method="jackknife", alpha_multiplier=nothing)

  multiple_frames = length(overlay)

  #-----
  treated = 1
  control = 1


  # crear la malla
  grid = DataFrame([repeat(x, outer=length(y)), vec(repeat(y, inner=length(x)))],
    [:x, :y])
  function description(row)
    est = estimates[1, row]
    # over = overlay
    se = if se_method == "none"
      NaN
    else
      sqrt(vcov(est, method=se_method))
    end
    setup = estimates.setup
    weights = estimates.weights
    Y = setup["Y"] .- contract3(setup["X"], weights["beta"])
    N0, N1 - setup$N0, size(Y, 1) - N0
    T0, T1 = setup$T0, size(Y, 2) - T0

    lambda_synth = hcat(weights["lambda"], fill(0, T1))
    lambda_target = hcat(fill(0, T0), fill(1 / T1, T1))
    omega_synth = hcat(weights["omega"], fill(0, N1))
    omega_target = hcat(fill(0, N0), fill(1 / N1, N1))

    if (!isnothing(estimates["overlay"]))
      over = estimates["overlay"]
    end
    is_sc = all(weights["lambda"] .== 0) || over .== 1

    intercept_offset - over * hcat(
      (omega_target - omega_synth) * Y * lambda_synth
    )


  end

  # plot_desctiptions = 


end


function synthdid_placebo_plot(estimate::synthdid_est1)

end

function synthdid_placebo_plot(estimate::synthdid_est1)

end

function synthdid_rmse_plot(estimate)

end