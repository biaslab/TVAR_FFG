using ForneyLab
using ForneyLab: unsafeMean, unsafeCov, unsafePrecision, unsafeMeanCov
using LAR
using LAR.Node
using ProgressMeter
using JLD
using Plots
using PGFPlotsX
pgfplotsx()
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepgfplotslibrary{fillbetween}")

# # NOTE Verification

d = load("data/verification/svmp.jld")
sfe = d["sfe"]
smx_set = d["smx_set"]
svx_set = d["svx_set"]
smθ_set = d["smθ_set"]
svθ_set = d["svθ_set"]
sγx_set = d["sγx_set"];

ds = load("data/verification/sdata.jld")
observations_set = ds["observations_set"]
process_set = ds["process_set"]
coefs_set = ds["coefs_set"]

index = 99 # dataset index

using LaTeXStrings
n_samples = length(process_set[index])
axis1 = @pgf Axis({xlabel=L"t",
           ylabel=L"x_t",
           legend_pos = "south east",
           legend_cell_align="{left}",
           grid = "major",
           title="Dataset $(index)",
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(1:n_samples), y = observations_set[index]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(1:n_samples), process_set[index])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(1:n_samples), smx_set[index])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(1:n_samples), smx_set[index] .+  sqrt.(svx_set[index]))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(1:n_samples), smx_set[index] .-  sqrt.(svx_set[index]))),
    Plot({ thick, color = "blue", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/stat_verification_estimation.tikz", axis1)

# Verification estimation TVAR

d = load("data/verification/tvvmp.jld")
tvfe = d["tvfe"]
tvmx_set = d["tvmx_set"]
tvvx_set = d["tvvx_set"]
tvmθ_set = d["tvmθ_set"]
tvvθ_set = d["tvvθ_set"]
tvγx_set = d["tvγx_set"];

ds = load("data/verification/tvdata.jld")
observations_set = ds["observations_set"]
process_set = ds["process_set"]
coefs_set = ds["coefs_set"]
coefs_set = ds["coefs_set"]

index = 10
n_samples = length(process_set[index])
axis2 = @pgf Axis({xlabel=L"t",
           ylabel=L"x_t",
           legend_pos = "north east",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="Dataset $(index)",
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(1:n_samples), y = observations_set[index]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(1:n_samples), process_set[index])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(1:n_samples), tvmx_set[index])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(1:n_samples), tvmx_set[index] .+  sqrt.(tvvx_set[index]))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(1:n_samples), tvmx_set[index] .-  sqrt.(tvvx_set[index]))),
    Plot({ thick, color = "blue", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/tv_verification_estimation.tikz", axis2)

index = 10
n_samples = length(process_set[index])
mθ = collect(Iterators.flatten(tvmθ_set[index]))
vθ = collect(Iterators.flatten(tvvθ_set[index]))
real_coefs = collect(Iterators.flatten(coefs_set[index]))
axis3 = @pgf Axis({xlabel=L"t",
           ylabel=L"\theta_t",
           legend_pos = "south east",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="Dataset $(index)",
    },
    Plot({no_marks,color="orange"}, Coordinates(collect(1:n_samples), real_coefs)), LegendEntry("real"),
    Plot({no_marks,color="blue", style ="{dashed}"}, Coordinates(collect(1:n_samples), mθ)),
    Plot({"name path=f", no_marks,color="blue",opacity=0.2 }, Coordinates(collect(1:n_samples), mθ .+  sqrt.(vθ))),
    Plot({"name path=g", no_marks, color="blue",opacity=0.2}, Coordinates(collect(1:n_samples), mθ .-  sqrt.(vθ))),
    Plot({ thick, color = "blue", fill = "blue", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/tv_coefficient_estimation.tikz", axis3)

# FE
n_dataset = size(sfe, 1)
vmp_iter = size(sfe, 2)
sFE = sum(sfe, dims=1) ./ n_dataset
sFE = reshape(FE, (10, ))
tvFE = sum(tvfe, dims=1) ./ n_dataset
tvFE = reshape(tvFE, (10, ))
n_samples = length(process_set[index])
axis4 = @pgf Axis({xlabel="iteration",
                  ylabel="Free energy [nats]",
                  legend_pos = "north east",
                  legend_cell_align="{left}",
                  scale = 1.0,
                  grid = "major",
    },
    Plot({mark = "square*", "blue"}, Coordinates(collect(1:vmp_iter), sFE)), LegendEntry("stationary AR"),
    Plot({mark = "o", "red"}, Coordinates(collect(1:vmp_iter), tvFE)), LegendEntry("time-varying AR"))

pgfsave("figures/tikz/verification_fe.tikz", axis4)


# NOTE: WEATHER
using Random
Random.seed!(42)
d = load("data/weather/weather.jld")
params_x_1, params_θ_1, params_η_1, params_γ_1, params_ξ_1, predictions_1, FE_1 = d["params_x_1"], d["params_θ_1"], d["params_η_1"], d["params_γ_1"], d["params_ξ_1"], d["predictions_1"], d["FE_1"]
params_x_2, params_θ_2, params_η_2, params_γ_2, params_ξ_2, predictions_2, FE_2 = d["params_x_2"], d["params_θ_2"], d["params_η_2"], d["params_γ_2"], d["params_ξ_2"], d["predictions_2"], d["FE_2"]
params_x_3, params_θ_3, params_η_3, params_γ_3, params_ξ_3, predictions_3, FE_3 = d["params_x_3"], d["params_θ_3"], d["params_η_3"], d["params_γ_3"], d["params_ξ_3"], d["predictions_3"], d["FE_3"]
params_x_4, params_θ_4, params_η_4, params_γ_4, params_ξ_4, predictions_4, FE_4 = d["params_x_4"], d["params_θ_4"], d["params_η_4"], d["params_γ_4"], d["params_ξ_4"], d["predictions_4"], d["FE_4"]
x = d["real"]
y = d["observed"]

n_samples = length(d["real"])

xmin=2000; xmax=2200
ar1 = @pgf Axis({xlabel="day",
           ylabel="temperature",
           legend_pos = "south east",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="AR(1)",
           xmin=xmin, xmax=xmax,
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(xmin:xmax), y = y[xmin:xmax]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(xmin:xmax), x[xmin:xmax])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(xmin:xmax), params_x_1[1][1, :][2:end][xmin:xmax])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(xmin:xmax), params_x_1[1][1, :][2:end][xmin:xmax] .+  sqrt.(inv.(params_x_1[2][1, 1, :][2:end][xmin:xmax])))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(xmin:xmax), params_x_1[1][1, :][2:end][xmin:xmax] .-  sqrt.(inv.(params_x_1[2][1, 1, :][2:end][xmin:xmax])))),
    Plot({ thick, color = "black", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/weather_ar1.tikz", ar1)


xmin=2000; xmax=2200
ar2 = @pgf Axis({xlabel="day",
           ylabel="temperature",
           legend_pos = "south east",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="AR(2)",
           xmin=xmin, xmax=xmax,
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(xmin:xmax), y = y[xmin:xmax]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(xmin:xmax), x[xmin:xmax])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(xmin:xmax), params_x_2[1][1, :][2:end][xmin:xmax])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(xmin:xmax), params_x_2[1][1, :][2:end][xmin:xmax] .+  sqrt.(inv.(params_x_2[2][1, 1, :][2:end][xmin:xmax])))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(xmin:xmax), params_x_2[1][1, :][2:end][xmin:xmax] .-  sqrt.(inv.(params_x_2[2][1, 1, :][2:end][xmin:xmax])))),
    Plot({ thick, color = "black", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/weather_ar2.tikz", ar2)


xmin=2000; xmax=2200
ar3 = @pgf Axis({xlabel="day",
           ylabel="temperature",
           legend_pos = "south east",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="AR(3)",
           xmin=xmin, xmax=xmax,
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(xmin:xmax), y = y[xmin:xmax]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(xmin:xmax), x[xmin:xmax])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(xmin:xmax), params_x_3[1][1, :][2:end][xmin:xmax])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(xmin:xmax), params_x_3[1][1, :][2:end][xmin:xmax] .+  sqrt.(inv.(params_x_3[2][1, 1, :][2:end][xmin:xmax])))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(xmin:xmax), params_x_3[1][1, :][2:end][xmin:xmax] .-  sqrt.(inv.(params_x_3[2][1, 1, :][2:end][xmin:xmax])))),
    Plot({ thick, color = "black", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/weather_ar3.tikz", ar3)


xmin=2000; xmax=2200
ar3 = @pgf Axis({xlabel="day",
           ylabel="temperature",
           legend_pos = "south east",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="AR(4)",
           xmin=xmin, xmax=xmax,
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(xmin:xmax), y = y[xmin:xmax]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(xmin:xmax), x[xmin:xmax])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(xmin:xmax), params_x_4[1][1, :][2:end][xmin:xmax])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(xmin:xmax), params_x_4[1][1, :][2:end][xmin:xmax] .+  sqrt.(inv.(params_x_4[2][1, 1, :][2:end][xmin:xmax])))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(xmin:xmax), params_x_4[1][1, :][2:end][xmin:xmax] .-  sqrt.(inv.(params_x_4[2][1, 1, :][2:end][xmin:xmax])))),
    Plot({ thick, color = "black", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/weather_ar4.tikz", ar4)


fe_1 = reshape(sum(FE_1, dims=1) ./ length(y), (10, ))
fe_2 = reshape(sum(FE_2, dims=1) ./ length(y), (10, ))
fe_3 = reshape(sum(FE_3, dims=1) ./ length(y), (10, ))
fe_4 = reshape(sum(FE_4, dims=1) ./ length(y), (10, ))

n_samples = length(process_set[index])
weather_vmp_iterations = @pgf Axis({xlabel="iteration",
                  ylabel="Free energy [nats]",
                  legend_pos = "outer north east",
                  legend_cell_align="{left}",
                  scale = 1.0,
                  grid = "major",
    },
    Plot({mark = "square*", "blue"}, Coordinates(collect(1:vmp_iter), fe_1)), LegendEntry("TVAR(1)"),
    Plot({mark = "triangle*", "green"}, Coordinates(collect(1:vmp_iter), fe_2)), LegendEntry("TVAR(2)"),
    Plot({mark = "*", "black"}, Coordinates(collect(1:vmp_iter), fe_3)), LegendEntry("TVAR(3)"),
    Plot({mark = "o", "red"}, Coordinates(collect(1:vmp_iter), fe_4)), LegendEntry("TVAR(4)"))

pgfsave("figures/tikz/weather_vmp_iterations.tikz", weather_vmp_iterations)

xmin, xmax = 591, 600
len = xmax - xmin + 1
xt = collect(xmin:xmax)
weather_vmp_time_frame = @pgf Axis({xlabel="last VMP iteration of a day",
           ylabel="Free Energy [nats]",
           legend_pos = "north east",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
    },
    Plot({mark = "square*", "blue"}, Coordinates(xt, FE_1[xmin:xmax, end])), LegendEntry("TVAR(1)"),
    Plot({mark = "triangle*", "green"}, Coordinates(xt, FE_2[xmin:xmax, end])), LegendEntry("TVAR(2)"),
    Plot({mark = "*", "black"}, Coordinates(xt, FE_3[xmin:xmax, end])), LegendEntry("TVAR(3)"),
    Plot({mark = "o", "red"}, Coordinates(xt, FE_4[xmin:xmax, end])), LegendEntry("TVAR(4)"))

pgfsave("figures/tikz/weather_vmp_time_frame.tikz", weather_vmp_time_frame)

len = xmax - xmin + 1
xt = collect(1:size(FE_1, 1))
weather_vmp_over_time = @pgf Axis({xlabel="last VMP iteration of a day",
           ylabel="Free Energy [nats]",
           legend_pos = "north east",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
    },
    Plot({no_marks, "blue"}, Coordinates(xt, FE_1[:, end])), LegendEntry("TVAR(1)"),
    Plot({no_marks, "green"}, Coordinates(xt, FE_2[:, end])), LegendEntry("TVAR(2)"),
    Plot({no_marks, "black"}, Coordinates(xt, FE_3[:, end])), LegendEntry("TVAR(3)"),
    Plot({no_marks, "red"}, Coordinates(xt, FE_4[:, end])), LegendEntry("TVAR(4)"))

pgfsave("figures/tikz/weather_vmp_over_time.tikz", weather_vmp_over_time)


# NOTE: Speech
d = load("data/speech/signal.jld")
x, y = d["real"], d["observed"]
x, y = reshape(x, (size(x, 1),)), reshape(y, (size(y, 1),))
cl, ns = x, y
# dividing into 10ms frames with 2.5ms overlap
fs = 8000
start = 1
l = Int(0.01*fs)
overlap = Int(0.0025*fs)
totseg = Int(ceil(length(ns)/(l-overlap)))
segment = zeros(totseg, l)
zseg = zeros(totseg, l)
for i in 1:totseg - 1
    global start
    segment[i,1:l]=ns[start:start+l-1]
    zseg[i, 1:l] = cl[start:start+l-1]
    start = (l-overlap)*i+1
end
segment[totseg, 1:length(ns)-start+1] = ns[start:length(ns)]
zseg[totseg, 1:length(cl)-start+1] = cl[start:length(cl)];

d = load("data/speech/inference.jld")
rmx_rw, rvx_rw, rγ_rw, fe_rw  = d["rmx_rw"], d["rvx_rw"], d["rγ_rw"], d["fe_rw"]
rmx_1, rvx_1, rmθ_1, rvθ_1, rγ_1, fe_1 = d["rmx_1"], d["rvx_1"], d["rmθ_1"], d["rvθ_1"], d["rγ_1"], d["fe_1"]
rmx_2, rvx_2, rmθ_2, rvθ_2, rγ_2, fe_2 = d["rmx_2"], d["rvx_2"], d["rmθ_2"], d["rvθ_2"], d["rγ_2"], d["fe_2"]
rmx_tv_1, rvx_tv_1, rmθ_tv_1, rvθ_tv_1, rγ_tv_1, fe_tv_1  = d["rmx_tv_1"], d["rvx_tv_1"], d["rmθ_tv_1"], d["rvθ_tv_1"], d["rγ_tv_1"], d["fe_tv_1"]
rmx_tv_2, rvx_tv_2, rmθ_tv_2, rvθ_tv_2, rγ_tv_2, fe_tv_2 = d["rmx_tv_2"], d["rvx_tv_2"], d["rmθ_tv_2"], d["rvθ_tv_2"], d["rγ_tv_2"], d["fe_tv_2"]

n_samples = length(process_set[index])
speech_free_energy = @pgf Axis({xlabel="segment number",
                  ylabel="Free energy [nats]",
                  legend_pos = "outer north east",
                  legend_cell_align="{left}",
                  scale = 1.0,
                  grid = "major",
    },
    Plot({mark = "square*", "blue"}, Coordinates(collect(1:totseg), fe_rw[:,end])), LegendEntry("RW"),
    Plot({mark = "triangle*", "green"}, Coordinates(collect(1:totseg), fe_1[:,end])), LegendEntry("AR(1)"),
    Plot({mark = "*", "black"}, Coordinates(collect(1:totseg), fe_2[:,end])), LegendEntry("AR(2)"),
    Plot({mark = "o", "red"}, Coordinates(collect(1:totseg), fe_tv_1[:,end])), LegendEntry("TVAR(1)"),
    Plot({mark = "-", "cyan"}, Coordinates(collect(1:totseg), fe_tv_2[:,end])), LegendEntry("TVAR(2)"))

pgfsave("figures/tikz/speech_free_energy.tikz", speech_free_energy)

# Advantage of TVAR(2) over RW
segnum = 293
n_samples = length(segment[segnum, :])
tvar2vsrw_a = @pgf Axis({xlabel="sample",
           ylabel="amplitude",
           legend_pos = "north west",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="TVAR(2)",
           ymax=0.15, ymin=-0.1,
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(1:n_samples), y = segment[segnum, :]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(1:n_samples), zseg[segnum, :])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(1:n_samples), rmx_tv_2[segnum, :])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(1:n_samples), rmx_tv_2[segnum, :] .+  sqrt.(rvx_tv_2[segnum, :]))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(1:n_samples), rmx_tv_2[segnum, :] .-  sqrt.(rvx_tv_2[segnum, :]))),
    Plot({ thick, color = "black", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/tvar2vsrw_a.tikz", tvar2vsrw_a)


tvar2vsrw_b = @pgf Axis({xlabel="sample",
           ylabel="amplitude",
           legend_pos = "north west",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="RW",
           ymax=0.15, ymin=-0.1,
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(1:n_samples), y = segment[segnum, :]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(1:n_samples), zseg[segnum, :])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(1:n_samples), rmx_rw[segnum, :])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(1:n_samples), rmx_rw[segnum, :] .+  sqrt.(rvx_rw[segnum, :]))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(1:n_samples), rmx_rw[segnum, :] .-  sqrt.(rvx_rw[segnum, :]))),
    Plot({ thick, color = "black", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/tvar2vsrw_b.tikz", tvar2vsrw_b)


n_samples = length(process_set[index])
speech_free_energy_tvar2vsrw = @pgf Axis({xlabel="segment number",
                  ylabel="Free energy [nats]",
                  legend_pos = "outer north east",
                  legend_cell_align="{left}",
                  scale = 1.0,
                  grid = "major",
                  xmin=segnum - 2, xmax=segnum+2,
    },
    Plot({mark = "square*", "blue"}, Coordinates(collect(segnum-2:segnum+2), fe_rw[segnum-2:segnum+2,end])), LegendEntry("RW"),
    Plot({mark = "triangle*", "green"}, Coordinates(collect(segnum-2:segnum+2), fe_1[segnum-2:segnum+2,end])), LegendEntry("AR(1)"),
    Plot({mark = "*", "black"}, Coordinates(collect(segnum-2:segnum+2), fe_2[segnum-2:segnum+2,end])), LegendEntry("AR(2)"),
    Plot({mark = "o", "red"}, Coordinates(collect(segnum-2:segnum+2), fe_tv_1[segnum-2:segnum+2,end])), LegendEntry("TVAR(1)"),
    Plot({mark = "x", "cyan"}, Coordinates(collect(segnum-2:segnum+2), fe_tv_2[segnum-2:segnum+2,end])), LegendEntry("TVAR(2)"))

pgfsave("figures/tikz/speech_free_energy_tvar2vsrw.tikz", speech_free_energy_tvar2vsrw)

# Advantage of AR(2) over RW
segnum = 208
n_samples = length(segment[segnum, :])
ar2vsrw_a = @pgf Axis({xlabel="sample",
           ylabel="amplitude",
           legend_pos = "north west",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="AR(2)",
           ymax=0.16, ymin=-0.15,
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(1:n_samples), y = segment[segnum, :]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(1:n_samples), zseg[segnum, :])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(1:n_samples), rmx_tv_2[segnum, :])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(1:n_samples), rmx_tv_2[segnum, :] .+  sqrt.(rvx_tv_2[segnum, :]))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(1:n_samples), rmx_tv_2[segnum, :] .-  sqrt.(rvx_tv_2[segnum, :]))),
    Plot({ thick, color = "black", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/ar2vsrw_a.tikz", ar2vsrw_a)


ar2vsrw_b = @pgf Axis({xlabel="sample",
           ylabel="amplitude",
           legend_pos = "north west",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="AR(2)",
           ymax=0.16, ymin=-0.15,
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(1:n_samples), y = segment[segnum, :]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(1:n_samples), zseg[segnum, :])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(1:n_samples), rmx_rw[segnum, :])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(1:n_samples), rmx_rw[segnum, :] .+  sqrt.(rvx_rw[segnum, :]))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(1:n_samples), rmx_rw[segnum, :] .-  sqrt.(rvx_rw[segnum, :]))),
    Plot({ thick, color = "black", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/ar2vsrw_b.tikz", ar2vsrw_b)


n_samples = length(process_set[index])
speech_free_energy_ar2vsrw = @pgf Axis({xlabel="segment number",
                  ylabel="Free energy [nats]",
                  legend_pos = "outer north east",
                  legend_cell_align="{left}",
                  scale = 1.0,
                  grid = "major",
                  xmin=segnum - 2, xmax=segnum+2,
    },
    Plot({mark = "square*", "blue"}, Coordinates(collect(segnum-2:segnum+2), fe_rw[segnum-2:segnum+2,end])), LegendEntry("RW"),
    Plot({mark = "triangle*", "green"}, Coordinates(collect(segnum-2:segnum+2), fe_1[segnum-2:segnum+2,end])), LegendEntry("AR(1)"),
    Plot({mark = "*", "black"}, Coordinates(collect(segnum-2:segnum+2), fe_2[segnum-2:segnum+2,end])), LegendEntry("AR(2)"),
    Plot({mark = "o", "red"}, Coordinates(collect(segnum-2:segnum+2), fe_tv_1[segnum-2:segnum+2,end])), LegendEntry("TVAR(1)"),
    Plot({mark = "x", "cyan"}, Coordinates(collect(segnum-2:segnum+2), fe_tv_2[segnum-2:segnum+2,end])), LegendEntry("TVAR(2)"))

pgfsave("figures/tikz/speech_free_energy_ar2vsrw.tikz", speech_free_energy_ar2vsrw)


# Advantage of RW over RW
segnum = 62
n_samples = length(segment[segnum, :])
rwvstvar2_a = @pgf Axis({xlabel="sample",
           ylabel="amplitude",
           legend_pos = "north west",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="TVAR(2)",
           ymax=0.15, ymin=-0.1,
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(1:n_samples), y = segment[segnum, :]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(1:n_samples), zseg[segnum, :])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(1:n_samples), rmx_tv_2[segnum, :])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(1:n_samples), rmx_tv_2[segnum, :] .+  sqrt.(rvx_tv_2[segnum, :]))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(1:n_samples), rmx_tv_2[segnum, :] .-  sqrt.(rvx_tv_2[segnum, :]))),
    Plot({ thick, color = "black", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/rwvstvar2_a.tikz", rwvstvar2_a)


rwvstvar2_b = @pgf Axis({xlabel="sample",
           ylabel="amplitude",
           legend_pos = "north west",
           legend_cell_align="{left}",
           scale = 1.0,
           grid = "major",
           title="RW",
           ymax=0.15, ymin=-0.1,
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = collect(1:n_samples), y = segment[segnum, :]
        ),
    ), LegendEntry("observations"),
    Plot({no_marks,color="magenta"}, Coordinates(collect(1:n_samples), zseg[segnum, :])), LegendEntry("real"),
    Plot({no_marks,color="black", style ="{dashed}"}, Coordinates(collect(1:n_samples), rmx_rw[segnum, :])),
    Plot({"name path=f", no_marks,color="black",opacity=0.2 }, Coordinates(collect(1:n_samples), rmx_rw[segnum, :] .+  sqrt.(rvx_rw[segnum, :]))),
    Plot({"name path=g", no_marks, color="black",opacity=0.2}, Coordinates(collect(1:n_samples), rmx_rw[segnum, :] .-  sqrt.(rvx_rw[segnum, :]))),
    Plot({ thick, color = "black", fill = "black", opacity = 0.2 },
            raw"fill between [of=f and g]"), LegendEntry("inferred"))

pgfsave("figures/tikz/rwvstvar2_b.tikz", rwvstvar2_b)

n_samples = length(process_set[index])
speech_free_energy_rwvstvar2 = @pgf Axis({xlabel="segment number",
                  ylabel="Free energy [nats]",
                  legend_pos = "outer north east",
                  legend_cell_align="{left}",
                  scale = 1.0,
                  grid = "major",
                  xmin=segnum - 2, xmax=segnum+2,
    },
    Plot({mark = "square*", "blue"}, Coordinates(collect(segnum-2:segnum+2), fe_rw[segnum-2:segnum+2,end])), LegendEntry("RW"),
    Plot({mark = "triangle*", "green"}, Coordinates(collect(segnum-2:segnum+2), fe_1[segnum-2:segnum+2,end])), LegendEntry("AR(1)"),
    Plot({mark = "*", "black"}, Coordinates(collect(segnum-2:segnum+2), fe_2[segnum-2:segnum+2,end])), LegendEntry("AR(2)"),
    Plot({mark = "o", "red"}, Coordinates(collect(segnum-2:segnum+2), fe_tv_1[segnum-2:segnum+2,end])), LegendEntry("TVAR(1)"),
    Plot({mark = "x", "cyan"}, Coordinates(collect(segnum-2:segnum+2), fe_tv_2[segnum-2:segnum+2,end])), LegendEntry("TVAR(2)"))

pgfsave("figures/tikz/speech_free_energy_rwvstvar2.tikz", speech_free_energy_rwvstvar2)

# Reconstructing the signal
FEs = [fe_rw[:, end], fe_1[:, end], fe_2[:, end], fe_tv_1[:, end], fe_tv_2[:, end]]
FEs_sorted = collect.(eachrow(Matrix(transpose(hcat(map(d -> sortperm(collect(d)), zip(FEs...))...)))))
findall(map(s -> s != [1, 2, 3, 4, 5], FEs_sorted));
rmx = [rmx_rw, rmx_1, rmx_2, rmx_tv_1, rmx_tv_2]
cleanSpeech = zeros(length(ns))
cleanSpeech[1:l] = rmx[FEs_sorted[1][1]][1, 1:l]

start = l + 1
for i in 2:totseg - 1
    global start
    cleanSpeech[start:start+(l-overlap)] = rmx[FEs_sorted[i][1]][i,overlap:end]
    start = start + l - overlap - 1
end
cleanSpeech[start:start+l-1] = rmx[FEs_sorted[totseg][1]][totseg,1:l]
cleanSpeech = cleanSpeech[1:length(ns)];
x = collect((length(ns)/fs)/length(ns):(length(ns)/fs)/length(ns):length(ns)/fs)
#cl = circshift(cl, 1)
xmin, xmax = 17050, 17200
reconstructed = @pgf Axis({xlabel=L"$t$ sec",
                  ylabel="Amplitude",
                  legend_pos = "outer north east",
                  legend_cell_align="{left}",
                  scale = 1.0,
                  grid="major", xmin=x[xmin], xmax=x[xmax],
    },
    Plot(
        {only_marks,color="black",opacity=0.6, mark="x"},
        Table(
            {x = "x", y = "y"},
             x = x[xmin:xmax], y = ns[xmin:xmax]
        ),
    ), LegendEntry("Noisy signal"),
    Plot({no_marks, "magenta"}, Coordinates(x[xmin:xmax], cl[xmin:xmax])), LegendEntry("Clean"),
    Plot({no_marks, "black", style ="{dashed}"}, Coordinates(x[xmin:xmax], cleanSpeech[xmin:xmax])), LegendEntry("Filtered"))

pgfsave("figures/tikz/reconstructed.tikz", reconstructed)
