using Statistics, DataFrames
N0, T0 = 3, 4
Y = reshape(1:100, 10, 10)

function collapse_form(Y, N0::Float16, T0::Float16)
  N, T = size(Matrix(Y))
  head = Y[1:N0, 1:T0]
  head_row_mean = mean(Y[1:N0, (T0+1):T], dims=2)
  head_matrix = hcat(head, head_row_mean)
  bottom_col_mean = mean(Y[(N0+1):N, 1:T0], dims=1)
  bottom = mean(Y[(N0+1):N, (T0+1):T])
  bottom_matrix = hcat(bottom_col_mean, bottom)
  return vcat(head_matrix, bottom_matrix)
end

# collapse_form(Y, N0, T0)

function pairwise_sum_decreasing(x::Vector, y::Vector)
  na_x = isnan.(x)
  na_y = isnan.(y)
  x[na_x] .= minimum(x[.!na_x])
  y[na_y] .= minimum(y[.!na_y])
  pairwise_sum = x .+ y
  pairwise_sum[na_x.&na_y] .= NaN
  return pairwise_sum
end
# x = [1,NaN, 2 , 3, 5]
# y = [10,NaN , 30, 40, 50]

# na_x = isnan.(x)
# x[na_x] = 12
# minimum(x[.!na_x])
# pairwise_sum_decreasing(x, y)


mutable struct panelMatrix
  Y::Array
  W::Array
  names::Vector
  time::Vector
  N0::Int
  T0::Int
end
function panel_matrices(panel::DataFrame; unit=1, time=2, outcome=3, treatment=4, treated_last=true)
  index_to_name(x) = x in 1:size(panel, 2) ? names(panel)[x] : x
  if any(ismissing.(eachrow(panel)))
    error("Missing values in `panel`.")
  end
  keep = [unit, time, outcome, treatment]
  panel = panel[:, keep]
  panel = sort(panel, [unit, time])

  unique_years = unique(panel[:, time])
  unique_units = unique(panel[:, unit])

  num_years = length(unique(panel[:, time]))
  num_units = length(unique(panel[:, unit]))

  Y = reshape(panel[:, outcome], num_years, num_units)'
  W = reshape(panel[:, treatment], num_years, num_units)'
  # panel = data("california_prop99")

  w = sum(W, dims=2)
  w = [x > 0 for x in w]
  N0 = num_units - sum(w)
  treat_time = any(W, dims=1)
  T0 = ([i for i in 1:length(treat_time) if treat_time[i]] |> first) - 1

  # if !(all(W[.!w, :] .== 0) && all(W[:, 1:T0] .== 0) && all(W[w, (T0+1):size(Y, 2)] .== 1))
  #     throw("The package cannot use this data. Treatment adoption is not simultaneous.")
  # end

  unit_order = if treated_last
    sortperm(W[:, T0+1])
  else
    collect(1:size(Y, 1))
  end


  panel = panelMatrix(Y[unit_order, :], W[unit_order, :], unique_units, unique_years, N0, T0)
  return panel
end

# panel1 = panel_matrices(data("california_prop99"));
# panel1.Y