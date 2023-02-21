---
title: "synthdid"
---



# synthdid: Synthetic Difference in Differences Estimation

This package implements the synthetic difference in difference estimator (SDID) for the average treatment effect in panel data,
as proposed in Arkhangelsky et al (2019). We observe matrices of outcomes Y and binary treatment indicators W
that we think of as satisfying Y<sub>ij</sub> = L<sub>ij</sub> + &tau;<sub>ij</sub> W<sub>ij</sub> + &epsilon;<sub>ij</sub>.
Here &tau;<sub>ij</sub> is the effect of treatment on the unit i at time j, and we estimate the average effect of
treatment when and where it happened: the average of &tau;<sub>ij</sub> over the observations with W<sub>ij</sub>=1.
All treated units must begin treatment simultaneously, so W is a block matrix: W<sub>ij</sub> = 1 for i > N<sub>0</sub> and j > T<sub>0</sub>
and zero otherwise, with N<sub>0</sub> denoting the number of control units and T<sub>0</sub> the number of observation times
before onset of treatment. This applies, in particular, to the case of a single treated unit or treated period.


This package is currently in beta and the functionality and interface is subject to change.

# Example

```julia
using Pkg
Pkg.add("synthdid")
using synthdid
```

```
Main.var"##WeaveSandBox#328".synthdid_est1(-15.603827872733847, Dict{String
, Union{Nothing, Vector{Float64}}}("omega_vals" => [13.18672026181369, 12.4
7377001377477, 12.106266678333789, 11.678082749493804, 11.558438925397132, 
11.246848929424237, 11.010734985812809, 10.979657913848131, 10.820174225624
315, 10.789142960356244  …  9.373927745928473, 9.373927263754181, 9.3739248
19153716, 9.373922880595307, 9.373922237848458, 9.37392037898267, 9.3739189
01210688, 9.373917235386077, 9.37391651650788, 9.373915017491015], "omega" 
=> [0.0, 0.003441829190052984, 0.057512787242020226, 0.0782872884996789, 0.
07036812272562899, 0.0015884105530842118, 0.03146821022946961, 0.0533878226
9887508, 0.0101354850429157, 0.02593862896097017  …  0.0, 0.004087093121374
664, 0.0, 0.009777529949684446, 0.04151765797797095, 0.0, 0.0, 0.0335691131
90658494, 0.03666708479079849, 0.0013864417998988702], "lambda_vals" => [77
.34951959277457, 77.34935942236041, 77.34935849346392, 77.34935848737783, 7
7.34935848733762, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0], "lambda" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.366470631936465, 0.20645305056017
524, 0.4270763175033598], "vals" => [90.53623985458826, 89.82312943613518, 
89.45562517179772, 89.02744123687162, 88.90779741273475, 11.246848929424237
, 11.010734985812809, 10.979657913848131, 10.820174225624315, 10.7891429603
56244  …  9.373927745928473, 9.373927263754181, 9.373924819153716, 9.373922
880595307, 9.373922237848458, 9.37392037898267, 9.373918901210688, 9.373917
235386077, 9.37391651650788, 9.373915017491015], "beta" => nothing), Dict{S
tring, Any}("Y" => [89.80000305 95.40000153 … 100.6999969 96.19999695; 100.
3000031 104.0999985 … 104.8000031 99.40000153; … ; 132.1999969 131.6999969 
… 104.8000031 90.5; 123.0 121.0 … 47.20000076 41.59999847], "X" => Array{An
y, 3}(undef, 39, 31, 0), "T0" => 19, "N0" => 38), Dict{String, Real}("updat
e_lambda" => true, "min_decrease" => 5.4944010185795114e-5, "max_iter" => 1
0000, "zeta_omega" => 10.226232571491236, "lambda_intercept" => true, "zeta
_lambda" => 5.494401018579511e-6, "omega_intercept" => true, "update_omega"
 => true), 38, 19)
```



```julia
setup_data = panel_matrices(data("california_prop99"))
tau_hat = synthdid_estimate(setup_data.Y, setup_data.N0, setup_data.T0)
summary_synth(tau_hat, panel = setup_data);
```

```
synthdid: -15.604 +- NaN. Effective N0/N0 = 16.388/38~0.431. Effective T0/T
0 = 2.783/19 ~ 0.146. N1,T1 = 1, 12.
```



```julia
p = synthdid_plot(tau_hat)
plot(p["plot"])
```

![](figures/README_4_1.png)

```julia
synthdid_units_plot(tau_hat, x_ticks = setup_data.names)
```

![](figures/README_5_1.png)
