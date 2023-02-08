using DataFrames, CSV
function data(file)
  # path = joinpath(dirname(pathof()))
  # path = pwd() * "\$(file).csv"
  path = joinpath(pwd(), "data\\$(file).csv")
  return DataFrame(CSV.read(path, DataFrame))
end


