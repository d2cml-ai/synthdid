using Statistics, CSV, DataFrames, Shuffle

function sum_normalize(x)
    if sum(x) != 0
        x/sum(x)
    else
        repeat([1/length(x)], length(x))
    end
end

function placebo_se(estimate, replications)
    setup = estimate.setup
    opts = estimate.opts
    weights = estimate.weight
    N1 = size(setup["Y"])[1] - setup["N0"]
    if setup["N0"] <= N1
        error("must have more controls than treated units to use the placebo se")
    end
    function theta(ind)
        N0 = length(ind)-N1
        weights1 = copy(weights)
        weights_boot = weights1
        weights_boot["omega"] = sum_normalize(weights["omega"][ind[1:N0]]);
        synthdid_estimate(setup["Y"][ind,:], N0, setup["T0"],  X=setup["X"][ind, :,:], weights=weights_boot,
                            zeta_omega = opts["zeta_omega"], zeta_lambda = opts["zeta_lambda"], omega_intercept = opts["omega_intercept"],
                            lambda_intercept = opts["lambda_intercept"], update_omega = opts["update_omega"], update_lambda = opts["update_lambda"],
                            min_decrease = opts["min_decrease"], max_iter = opts["max_iter"])
    end
    
    a = 0
    for i in 1:replications
        if i == 1
            a = theta(shuffle(1:setup["N0"])).estimate
        else
            a = vcat(a,theta(shuffle(1:setup["N0"])).estimate)
        end
    end
    sqrt((replications-1)/replications) * std(a)
end

function vcov_synthdid_estimate(object; method = ["bootstrap", "jackknife", "placebo"], replications = 200)
    if method == "placebo"
        se = placebo_se(object, replications)
    end
    se^2
end