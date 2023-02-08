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


function panel_matrices(panel::DataFrame; unit=1, time=2, outcome=3, treatment=4, treated_last=true)
  index_to_name(x) = x in 1:size(panel, 2) ? names(panel)[x] : x
  if any(ismissing.(eachrow(panel)))
    error("Missing values in `panel`.")
  end
  # if length(unique(panel[:, treatment])) == 1
  #   error("There is no variation in treatment status.")
  # end
  # if !all(ismember(panel[:, treatment], [0, 1]))
  #   error("The treatment status should be in 0 or 1.")
  # end
  # panel = DataFrame(map(col ->
  #     begin
  #       if isa(col, CategoricalArray) || isa(col, Date)
  #         convert(Vector{String}, col)
  #       else
  #         convert(Vector, col)
  #       end
  #     end,
  #   panel)
  # )
  sort!(panel, [:x2, :x1])
  return panel

end

y_panel = DataFrame(reshape(rand(100), 10, 10), :auto)

panel_matrices(y_panel)

