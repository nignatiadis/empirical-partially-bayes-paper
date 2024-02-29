using FileIO
using CategoricalArrays
using DataFrames
using Plots
using StatsPlots
using StatsBase
using LaTeXStrings
using CSV


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
data_dir = joinpath(dir, "data")

# List all .jld2 files in the directory
files = filter(x -> occursin(r".jld2$", x), readdir(data_dir))

# Initialize an empty DataFrame to concatenate all results
method_res = DataFrame()

# Loop through each file, load its DataFrame, and concatenate
for file in files
    filepath = joinpath(data_dir, file)
    df = load(filepath, "method_res")  # Adjust "method_res" if the variable name differs
    method_res = vcat(method_res, df)
end

names(method_res)


function sanitize_data(method_res; keep_npmle_variants = false)
    summary_tbl =
        method_res |>
        x -> groupby(x, [:signal_setting, :variance_setting, :ν, :ordering, :method_name]) |>
        x -> DataFrames.combine(x, nrow, [:FDP, :Power, :cond_cov] .=> mean, [:FDP, :Power] .=> std)

    if keep_npmle_variants
        npmle_methods = [:t_npmle_BH_0; :t_npmle_BH; :t_npmle_BH_02]
        summary_tbl = filter(row -> row.method_name ∈ npmle_methods, summary_tbl)
        summary_tbl.cleaned_method_names = CategoricalArray(string.(summary_tbl.method_name))
        recode!(
            summary_tbl.cleaned_method_names,
            "t_npmle_BH" => "NPMLE(0.01)",
            "t_npmle_BH_02" => "NPMLE(0.2)",
            "t_npmle_BH_0" => "NPMLE(min)",
        )
        levels!(
            summary_tbl.cleaned_method_names,
            ["NPMLE(min)"; "NPMLE(0.01)"; "NPMLE(0.2)"]
        )

    else
        cleaned_method_names = CategoricalArray(string.(summary_tbl.method_name))
        recode!(
            cleaned_method_names,
            "t_npmle_BH" => "NPMLE (BH)",
            "t_limma_BH" => "Limma (BH)",
            "t_standard_BH" => "t-test (BH)",
            "t_oracle_BH" => "Oracle (BH)",
            "t_oracle_storey" => "Oracle (Storey)",
        )
        primary_method_list = ["t-test (BH)", "Limma (BH)", "NPMLE (BH)", "Oracle (BH)", "Oracle (Storey)"]
        summary_tbl.cleaned_method_names = cleaned_method_names
        summary_tbl = filter(row -> row.cleaned_method_names in primary_method_list, summary_tbl)
        levels!(
            summary_tbl.cleaned_method_names,
            primary_method_list
        )
        summary_tbl.code = Int.(summary_tbl.cleaned_method_names.refs)
    end
    summary_tbl.ν_factor = CategoricalArray(string.(summary_tbl.ν))
    levels!(summary_tbl.ν_factor, ["2", "4", "8", "16", "32", "64"])
 
    return summary_tbl
end

_default_colors = [colorant"#4C413F" colorant"#278B9A" colorant"#E75B64" colorant"#D8AF39" colorant"#E8C4A2"]

function plot_groupedbar(summary_tbl::DataFrame,  
    variance_setting::Symbol,
    signal_setting::Symbol,
    ordering_setting::Symbol;
    legend = true,
    cond_cov_legend = false,
    method_colors = _default_colors
    )

    variance_setting_name_map =   Dict(
        :Dirac => L"\textrm{Dirac}",
        :InverseGamma => L"\textrm{Scaled inverse } \chi^2",
        :TwoPointPrior => L"\textrm{Two-point}",
    )


    filtered_data = filter(
        row -> (row.signal_setting == signal_setting) &
              (row.variance_setting == variance_setting) &
              (row.ordering == ordering_setting),
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
        legend = legend ? :top : nothing,
        alpha = 0.7,
        legend_column = 3,
        xlabel = L"\nu",
        ylabel = "FDR",
        title = variance_setting_name_map[variance_setting],
    )

    hline!(fdr_plot, [0.1], linestyle = :dash, color = :grey, label = "")

    power_plot = @df filtered_data groupedbar(
        :ν_factor,
        :Power_mean,
        group = :cleaned_method_names,
        ylim = (0, 1.01),
        yticks = 0:0.2:1,
        color = method_colors,
        linestyle = :solid,
        legend = nothing,
        alpha = 0.7,
        legend_column = 3,
        xlabel = L"\nu",
        ylabel = "Power",
        title =  variance_setting_name_map[variance_setting],
    )

    cond_cov_plot = @df filtered_data groupedbar(
        :ν_factor,
        :cond_cov_mean,
        group = :cleaned_method_names,
        ylim = (0, 1.01),
        yticks = 0:0.2:1,
        color = method_colors,
        linestyle = :solid,
        legend = cond_cov_legend ? :top : nothing,
        alpha = 0.7,
        legend_column = 3,
        xlabel = L"\nu",
        ylabel = L"\textrm{MinSVarFP}(\leq 0.2)",
        title =  variance_setting_name_map[variance_setting],
    )

    hline!(cond_cov_plot, [0.2], linestyle = :dash, color = :grey, label = "")


    return fdr_plot, power_plot, cond_cov_plot
end

summary_tbl = sanitize_data(method_res)

nonadversarial_conjugate_plots = []

for variance_setting in [ :Dirac; :InverseGamma; :TwoPointPrior]
    push!(nonadversarial_conjugate_plots, 
        plot_groupedbar(summary_tbl, variance_setting,
        :conjugate,
        :random)...)
end


nonadversarial_conjugate_plot = plot(nonadversarial_conjugate_plots..., 
        layout =  @layout([a b c ; c d e; e f g]),    #@layout([a b ; c d; e f]), 
        margin = 3 * Plots.mm,
        size = (1400, 900))
        

savefig(nonadversarial_conjugate_plot, "main_simulations.pdf")

adversarial_conjugate_plots = []

for variance_setting in [:InverseGamma; :TwoPointPrior]
    push!(adversarial_conjugate_plots, 
    plot_groupedbar(summary_tbl, variance_setting,
    :conjugate,
    :decreasing)...)
end

adversarial_conjugate_plot = plot(adversarial_conjugate_plots..., 
        layout = @layout([a b c; d e f]), 
        margin = 3 * Plots.mm,
        size = (1400, 650))


savefig(adversarial_conjugate_plot, "adv_simulations.pdf")


nonadversarial_normal_plots = []

for variance_setting in [:InverseGamma; :TwoPointPrior]
    push!(nonadversarial_normal_plots, 
        plot_groupedbar(summary_tbl, variance_setting,
        :normal,
        :random)...)
end

nonadversarial_normal_plot = plot(nonadversarial_normal_plots..., 
    layout = @layout([a b c; d e f]), 
    margin = 3 * Plots.mm,
    size = (1400, 650))

savefig(nonadversarial_normal_plot, "suppl_simulations_normal_signal_strength.pdf")


nonadversarial_dirac_plots = []

for variance_setting in [:Dirac; :InverseGamma; :TwoPointPrior]
    push!(nonadversarial_dirac_plots, 
        plot_groupedbar(summary_tbl, variance_setting,
        :point_mass,
        :random)...)
end


nonadversarial_dirac_plot = plot(nonadversarial_dirac_plots..., 
        layout =  @layout([a b c ; c d e; e f g]),    #@layout([a b ; c d; e f]), 
        margin = 3 * Plots.mm,
        size = (1400, 900))
        
savefig(nonadversarial_dirac_plot, "suppl_simulations_dirac_signal_strength.pdf")

    

adversarial_normal_plots = []

for variance_setting in [:InverseGamma; :TwoPointPrior]
    push!(adversarial_normal_plots, 
        plot_groupedbar(summary_tbl, variance_setting,
        :normal,
        :decreasing)...)
end

adversarial_normal_plot = plot(nonadversarial_normal_plots..., 
    layout = @layout([a b c; d e f]), 
    margin = 3 * Plots.mm,
    size = (1400, 650))


adversarial_dirac_plots = []

for variance_setting in [:InverseGamma; :TwoPointPrior]
    push!(adversarial_dirac_plots, 
        plot_groupedbar(summary_tbl, variance_setting,
        :point_mass,
        :decreasing)...)
end


adversarial_dirac_plot = plot(adversarial_dirac_plots..., 
        layout =  @layout([a b c ; c d e; e f g]),    #@layout([a b ; c d; e f]), 
        margin = 3 * Plots.mm,
        size = (1400, 900))

summary_tbl_npmle = sanitize_data(method_res; keep_npmle_variants = true)
summary_tbl_npmle = filter(row->row.signal_setting == :conjugate, summary_tbl_npmle)

summary_tbl_npmle_FDR = unstack(summary_tbl_npmle, 
    [:signal_setting, :variance_setting, :ν, :ordering],
    :cleaned_method_names,
    :FDP_mean,
    renamecols=x->Symbol(:FDR, :_, x))


maximum(abs.(summary_tbl_npmle_FDR.FDR_t_npmle_BH .- summary_tbl_npmle_FDR.FDR_t_npmle_BH_0))
maximum(abs.(summary_tbl_npmle_FDR.FDR_t_npmle_BH .- summary_tbl_npmle_FDR.FDR_t_npmle_BH_02))
maximum(abs.(summary_tbl_npmle_FDR.FDR_t_npmle_BH_0 .- summary_tbl_npmle_FDR.FDR_t_npmle_BH_02))



summary_tbl_npmle_power = unstack(summary_tbl_npmle, 
    [:signal_setting, :variance_setting, :ν, :ordering],
    :method_name,
    :Power_mean,
    renamecols=x->"Power_$(x)")

summary_tbl_npmle_power.power_diff = max.( abs.(summary_tbl_npmle_power.Power_t_npmle_BH .- summary_tbl_npmle_power.Power_t_npmle_BH_02),
                                           abs.(summary_tbl_npmle_power.Power_t_npmle_BH .- summary_tbl_npmle_power.Power_t_npmle_BH_0))

maximum(abs.(summary_tbl_npmle_power.Power_t_npmle_BH .- summary_tbl_npmle_power.Power_t_npmle_BH_0))
maximum(abs.(summary_tbl_npmle_power.Power_t_npmle_BH .- summary_tbl_npmle_power.Power_t_npmle_BH_02))
maximum(abs.(summary_tbl_npmle_power.Power_t_npmle_BH_0 .- summary_tbl_npmle_power.Power_t_npmle_BH_02))

CSV.write("summary_tbl_npmle_power.csv", sort(summary_tbl_npmle_power, order(:power_diff, rev=true))[1:4,:])

summary_tbl_npmle_cond = unstack(summary_tbl_npmle, 
    [:signal_setting, :variance_setting, :ν, :ordering],
    :method_name,
    :cond_cov_mean,
    renamecols=x->Symbol(:Cond, :_, x))

sensitivity_plot = plot(
 getindex.(plot_groupedbar.(Ref(summary_tbl_npmle),
   [:Dirac; :InverseGamma; :TwoPointPrior],
   :conjugate,
   :random; 
   cond_cov_legend=true,
   method_colors = [:orange _default_colors[3] :purple]), 
   3)...,
   layout = (1,3),
   margin = 3 * Plots.mm,
   size = (1400, 400))

savefig(sensitivity_plot, "sensitivity_plot_conditionality.pdf")

adversarial_sensitivity_plot = plot(
    plot(legend=false, grid=false, axis=false, ticks=nothing, border=:none, margin=100 * Plots.mm),
 getindex.(plot_groupedbar.(Ref(summary_tbl_npmle),
   [:InverseGamma; :TwoPointPrior],
   :conjugate,
   :decreasing; 
   cond_cov_legend=true,
   method_colors = [:orange _default_colors[3] :purple]), 
   3)...,
   layout = (1,3),
   margin = 3 * Plots.mm,
   size = (1400, 400))

savefig(adversarial_sensitivity_plot, "adversarial_sensitivity_plot_conditionality.pdf")


maximum(abs.(summary_tbl_npmle_cond.Cond_t_npmle_BH .- summary_tbl_npmle_cond.Cond_t_npmle_BH_0))
maximum(abs.(summary_tbl_npmle_power.Power_t_npmle_BH .- summary_tbl_npmle_power.Power_t_npmle_BH_02))
maximum(abs.(summary_tbl_npmle_power.Power_t_npmle_BH_0 .- summary_tbl_npmle_power.Power_t_npmle_BH_02))
