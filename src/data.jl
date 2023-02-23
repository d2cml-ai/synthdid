
function data(file::String)
  # path = joinpath(dirname(pathof()))
  # path = pwd() * "\$(file).csv"
  path = joinpath(pwd(), "data\\$(file).csv")
  return DataFrame(CSV.read(path, DataFrame))
end

function data1(file::String, jb=true)
  if jb
    path = joinpath("..\\data\\$(file).csv")
    print(path)
    return DataFrame(CSV.read(path, DataFrame))
  end
end


# data("california_prop99")