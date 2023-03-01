
include("setup.jl")

quota = data("quota_example_sample")

quota = sort(quota, [:year, :quota, :country])

quota[:, ["womparl", "country", "year", "quota"]]