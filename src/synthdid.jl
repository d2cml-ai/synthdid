function sparsify_function(v)
  v[v.<=maximum(v)/4] .= 0
  return v ./ sum(v)
end
# v = [1, 2, 4, 5, 7, 2, 3, 9]
# v[v .<= maximum(v) / 4] .= 0

# v