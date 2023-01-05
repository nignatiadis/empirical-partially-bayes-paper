using FileIO
using CategoricalArrays
using DataFrames
using Plots
using StatsPlots
using StatsBase
using LaTeXStrings

pgfplotsx()
theme(
    :default;
    background_color_legend = :transparent,
    foreground_color_legend = :transparent,
    grid = nothing,
    frame = :box,
    legendfonthalign = :left,
    thickness_scaling = 1.3,
    size = (440, 330),
)




method_res1 = load("simulations/method_res_x.jld2", "method_res")
method_res2 = load("simulations/method_res_twopoint.jld2", "method_res")

method_res = vcat(method_res1, method_res2)

summary_tbl =
    method_res |>
    x ->
        groupby(x, [:simulation_setting, :ν, :adversarial_ordering, :method_name]) |>
        x -> DataFrames.combine(x, nrow, [:FDP, :Power] .=> mean, [:FDP, :Power] .=> std)

unique(summary_tbl.simulation_setting)

cleaned_method_names = CategoricalArray(string.(summary_tbl.method_name))
recode!(
    cleaned_method_names,
    "t_npmle_BH" => "NPMLE (BH)",
    "t_limma_BH" => "Limma (BH)",
    "t_standard_BH" => "t-test (BH)",
    "t_oracle_BH" => "Oracle (BH)",
    "t_oracle_storey" => "Oracle (Storey)",
)

levels!(
    cleaned_method_names,
    ["t-test (BH)", "Limma (BH)", "NPMLE (BH)", "Oracle (BH)", "Oracle (Storey)"],
)
summary_tbl.cleaned_method_names = cleaned_method_names
summary_tbl.code = Int.(summary_tbl.cleaned_method_names.refs)
summary_tbl.ν_factor = CategoricalArray(string.(summary_tbl.ν))
levels!(summary_tbl.ν_factor, ["2", "4", "8", "16", "32", "64"])


#method_shapes = [:circle  :dtriangle :rtriangle :ltriangle :utriangle]
method_colors =
    [colorant"#4C413F" colorant"#278B9A" colorant"#E75B64" colorant"#D8AF39" colorant"#E8C4A2"]
#method_linestyles = [:solid :dash :dot :dashdot :dashdotdot]

dirac_nonadv = filter(
    row -> (row.simulation_setting == :Dirac) & !row.adversarial_ordering,
    summary_tbl,
)

dirac_fdr = @df dirac_nonadv groupedbar(
    :ν_factor,
    :FDP_mean,
    group = :cleaned_method_names,
    ylim = (0, 0.148),
    yticks = 0:0.02:0.12,
    color = method_colors,
    linestyle = :solid,
    legend = :top,
    alpha = 0.7,
    legend_column = 3,
)
hline!(
    dirac_fdr,
    [0.1],
    linestyle = :dash,
    color = :grey,
    label = "",
    xlabel = L"\nu",
    ylabel = "FDR",
    title = "Dirac",
)

dirac_power = @df dirac_nonadv groupedbar(
    :ν_factor,
    :Power_mean,
    group = :cleaned_method_names,
    ylim = (0, 1),
    color = method_colors,
    linestyle = :solid,
    legend = :top,
    alpha = 0.7,
    legend_column = 3,
    xlabel = L"\nu",
    ylabel = "Power",
    title = "Dirac",
)

inversegamma_nonadv = filter(
    row -> (row.simulation_setting == :InverseGamma) & !row.adversarial_ordering,
    summary_tbl,
)


inversegamma_nonadv_fdr = @df inversegamma_nonadv groupedbar(
    :ν_factor,
    :FDP_mean,
    group = :cleaned_method_names,
    ylim = (0, 0.148),
    yticks = 0:0.02:0.12,
    color = method_colors,
    linestyle = :solid,
    legend = :top,
    alpha = 0.7,
    legend_column = 3,
)

hline!(
    inversegamma_nonadv_fdr,
    [0.1],
    linestyle = :dash,
    color = :grey,
    label = "",
    xlabel = L"\nu",
    ylabel = "FDR",
    title = L"\textrm{Scaled inverse } \chi^2",
)

inversegamma_nonadv_power = @df inversegamma_nonadv groupedbar(
    :ν_factor,
    :Power_mean,
    group = :cleaned_method_names,
    ylim = (0, 1),
    color = method_colors,
    linestyle = :solid,
    legend = :top,
    alpha = 0.7,
    legend_column = 3,
    xlabel = L"\nu",
    ylabel = "Power",
    title = L"\textrm{Scaled inverse } \chi^2",
)


twopoint_nonadv = filter(
    row -> (row.simulation_setting == :TwoPointPrior3) & !row.adversarial_ordering,
    summary_tbl,
)

twopoint_nonadv_fdr = @df twopoint_nonadv groupedbar(
    :ν_factor,
    :FDP_mean,
    group = :cleaned_method_names,
    ylim = (0, 0.148),
    yticks = 0:0.02:0.12,
    color = method_colors,
    linestyle = :solid,
    legend = :top,
    alpha = 0.7,
    legend_column = 3,
)
hline!(
    twopoint_nonadv_fdr,
    [0.1],
    linestyle = :dash,
    color = :grey,
    label = "",
    xlabel = L"\nu",
    ylabel = "FDR",
    title = "Two Point",
)

twopoint_nonadv_power = @df twopoint_nonadv groupedbar(
    :ν_factor,
    :Power_mean,
    group = :cleaned_method_names,
    ylim = (0, 1),
    color = method_colors,
    linestyle = :solid,
    legend = :top,
    alpha = 0.7,
    legend_column = 3,
    xlabel = L"\nu",
    ylabel = "Power",
    title = "Two Point",
)

nonadv_plot = plot(
    dirac_fdr,
    dirac_power,
    inversegamma_nonadv_fdr,
    inversegamma_nonadv_power,
    twopoint_nonadv_fdr,
    twopoint_nonadv_power,
    layout = (3, 2),
    size = (1200, 900),
    margin = 3 * Plots.mm,
)

savefig(nonadv_plot, "main_simulations.pdf")


inversegamma_adv = filter(
    row -> (row.simulation_setting == :InverseGamma) & row.adversarial_ordering,
    summary_tbl,
)
inversegamma_adv_fdr = @df inversegamma_adv groupedbar(
    :ν_factor,
    :FDP_mean,
    group = :cleaned_method_names,
    ylim = (0, 0.158),
    yticks = 0:0.02:0.12,
    color = method_colors,
    linestyle = :solid,
    legend = :top,
    alpha = 0.7,
    legend_column = 3,
)
hline!(
    inversegamma_adv_fdr,
    [0.1],
    linestyle = :dash,
    color = :grey,
    label = "",
    xlabel = L"\nu",
    ylabel = "FDR",
    title = L"\textrm{Scaled inverse } \chi^2",
)

inversegamma_adv_power = @df inversegamma_adv groupedbar(
    :ν_factor,
    :Power_mean,
    group = :cleaned_method_names,
    ylim = (0, 1),
    color = method_colors,
    linestyle = :solid,
    legend = :top,
    alpha = 0.7,
    legend_column = 3,
    xlabel = L"\nu",
    ylabel = "Power",
    title = L"\textrm{Scaled inverse } \chi^2",
)


twopoint_adv = filter(
    row -> (row.simulation_setting == :TwoPointPrior3) & row.adversarial_ordering,
    summary_tbl,
)

twopoint_adv_fdr = @df twopoint_adv groupedbar(
    :ν_factor,
    :FDP_mean,
    group = :cleaned_method_names,
    ylim = (0, 0.158),
    yticks = 0:0.02:0.12,
    color = method_colors,
    linestyle = :solid,
    legend = :top,
    alpha = 0.7,
    legend_column = 3,
)
hline!(
    twopoint_adv_fdr,
    [0.1],
    linestyle = :dash,
    color = :grey,
    label = "",
    xlabel = L"\nu",
    ylabel = "FDR",
    title = "Two Point",
)

twopoint_adv_power = @df twopoint_adv groupedbar(
    :ν_factor,
    :Power_mean,
    group = :cleaned_method_names,
    ylim = (0, 1),
    color = method_colors,
    linestyle = :solid,
    legend = :top,
    alpha = 0.7,
    legend_column = 3,
    xlabel = L"\nu",
    ylabel = "Power",
    title = "Two Point",
)


adv_plot = plot(
    inversegamma_adv_fdr,
    inversegamma_adv_power,
    twopoint_adv_fdr,
    twopoint_adv_power,
    layout = (2, 2),
    size = (1200, 650),
    margin = 3 * Plots.mm,
)

savefig(adv_plot, "adv_simulations.pdf")
