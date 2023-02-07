function timesteps(Y)
  try
    names(Y)
  catch
    function error(e)
      names(Y)
    end
  end
end