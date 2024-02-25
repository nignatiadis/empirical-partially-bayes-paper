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


dir = @__DIR__
method_res = load(joinpath(dir, "method_res_3k.jld2"), "method_res")
names(method_res)


function sanitize_data(method_res)
    summary_tbl =
        method_res |>
        x -> groupby(x, [:simulation_setting, :ν, :adversarial_ordering, :method_name]) |>
        x -> DataFrames.combine(x, nrow, [:FDP, :Power] .=> mean, [:FDP, :Power] .=> std)

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
 
    return summary_tbl
end


function plot_groupedbar(summary_tbl::DataFrame, setting::Symbol, adv_ordering::Bool)

    simulation_name_map =   Dict(
        :Dirac => L"\textrm{Dirac}",
        :InverseGamma => L"\textrm{Scaled inverse } \chi^2",
        :TwoPointPrior => L"\textrm{Two-point}",
    )

    method_colors =
    [colorant"#4C413F" colorant"#278B9A" colorant"#E75B64" colorant"#D8AF39" colorant"#E8C4A2"]

    filtered_data = filter(
        row -> (row.simulation_setting == setting) & (row.adversarial_ordering == adv_ordering),
        summary_tbl,
    )

    fdr_plot = @df filtered_data groupedbar(
        :ν_factor,
        :FDP_mean,
        group = :cleaned_method_names,
        ylim = (0, 0.158), # Adjust based on your data
        yticks = 0:0.02:0.12,
        color = method_colors,
        linestyle = :solid,
        legend = :top,
        alpha = 0.7,
        legend_column = 3,
        xlabel = L"\nu",
        ylabel = "FDR",
        title = simulation_name_map[setting],
    )

    hline!(fdr_plot, [0.1], linestyle = :dash, color = :grey, label = "")

    power_plot = @df filtered_data groupedbar(
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
        title =  simulation_name_map[setting],
    )

    return fdr_plot, power_plot
end

summary_tbl = sanitize_data(method_res)

simulation_settings = unique(summary_tbl.simulation_setting)


nonadversarial_plots = []

for setting in simulation_settings
    push!(nonadversarial_plots, plot_groupedbar(summary_tbl, setting, false)...)
end


nonadv_plot = plot(nonadversarial_plots..., 
        layout = @layout([a b; c d; e f]), 
        margin = 3 * Plots.mm,
        size = (1200, 900))
        

savefig(nonadv_plot, "main_simulations.pdf")

adversarial_plots = []

for setting in [:InverseGamma, :TwoPointPrior]
    push!(adversarial_plots, plot_groupedbar(summary_tbl, setting, true)...)
end

adv_plot = plot(adversarial_plots..., 
        layout = @layout([a b; c d]), 
        margin = 3 * Plots.mm,
        size = (1200, 650))


savefig(adv_plot, "adv_simulations.pdf")