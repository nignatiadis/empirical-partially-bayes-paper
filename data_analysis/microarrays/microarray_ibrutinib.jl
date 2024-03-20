using Pkg

dir = @__DIR__
Pkg.activate(joinpath(dir, "..","..","simulations"))

using Empirikos
using CSV
using LaTeXStrings
using Statistics
using Plots
using MosekTools
using DataFrames

include(joinpath(dir, "..", "..", "helpers.jl"))

pgfplotsx()
theme(
    :default;
    background_color_legend = :transparent,
    foreground_color_legend = :transparent,
    grid = nothing,
    frame = :box,
    legendfonthalign = :left,
    thickness_scaling = 1.3,
    size = (420, 330),
)


microarray = CSV.File(joinpath(dir, "microarray_ibrutinib.csv")) |> DataFrame
microarray = filter(row -> row.residual_dof == 11, microarray)


# Code to be added to Empirikos.jl tests
# limma_ν_prior = 10.42875
# limma_prior_var_naive_vs_act = 0.007079249

microarray_Ss =
    ScaledChiSquareSample.(abs2.(microarray.se_hat), microarray.residual_dof)
microarray_mu_hat = microarray.mu_hat
microarray_samples = NormalChiSquareSample.(microarray.mu_hat, microarray_Ss)

microarray_npmle = fit(
    Empirikos.EmpiricalPartiallyBayesTTest(
        prior = DiscretePriorClass(),
        solver = Mosek.Optimizer,
        α = 0.05,
        discretize_marginal = false,
    ),
    microarray_samples,
)

microarray_limma =
    fit(Empirikos.EmpiricalPartiallyBayesTTest(
        prior = Empirikos.Limma(),
        solver = nothing,
        α = 0.05),
        microarray_samples)

# Comparison to limma results from R package
# limma_prior_var_naive_vs_act = 0.007078982
# limma_ν_prior = 10.43969
microarray_limma.prior

microarray_t = fit(SimultaneousTTest(α=0.05), microarray_samples)

# Rejections by method
(
    microarray_npmle.total_rejections,
    microarray_limma.total_rejections,
    microarray_t.total_rejections,
)
#---------------------------------------------------------
# Panel A: Marginal density of variances
#---------------------------------------------------------

extrema(support(microarray_npmle.prior))

var_grid_refine = 0.0001:0.0001:0.05
var_grid_refine_samples = ScaledChiSquareSample.(var_grid_refine, 11)

marginal_pdf_limma = pdf.(microarray_limma.prior, var_grid_refine_samples)
marginal_pdf_npmle = pdf.(microarray_npmle.prior, var_grid_refine_samples)

histogram_plot = Plots.histogram(
    abs2.(microarray.se_hat)[abs2.(microarray.se_hat).<0.1],
    fill = "lightgrey",
    normalize = true,
    fillalpha = 0.5,
    linewidth = 0.3,
    linealpha = 0.5,
    bins = 100,
    xlim = (0, 0.05),
    ylim = (0, 140),
    label = "Histogram",
    xlabel = L"S_i^2",
    ylabel = "Density",
    legend = :topright,
)

Plots.plot!(
    histogram_plot,
    var_grid_refine,
    [marginal_pdf_npmle marginal_pdf_limma],
    label = ["NPMLE" "Limma"],
    color = [:purple :darkorange],
    linewidth = 0.8,
)

#---------------------------------------------------------
# Panel B: Prior
#---------------------------------------------------------

prior_plot = Plots.plot(
    support(microarray_npmle.prior),
    300 * probs(microarray_npmle.prior),
    seriestype = :sticks,
    xlim = (0, 0.05),
    legend = :topright,
    label = "NPMLE",
    color = :purple,
    xlabel = L"\sigma^2",
    ylabel = L"g(\sigma^2)",
)

Plots.plot!(
    prior_plot,
    u -> pdf(microarray_limma.prior, u),
    label = "Limma",
    color = :darkorange,
)


#---------------------------------------------------------
# Panel C: 2D-rejection regions
#---------------------------------------------------------

extrema(log10.(abs2.(microarray.se_hat)))


log_grid = [0.001; 0.01; 0.1]

twod_histogram_plot = histogram2d(
    log.(abs2.(microarray.se_hat)),
    microarray_mu_hat,
    bins = 50,
    c = cgrad(:algae, rev = false, scale = :exp),
    xlabel = L"S_i^2",
    ylabel = L"Z_i",
    xticks = (log.(log_grid), string.(log_grid)),
)

equidistant_grid = 0.0:0.0002:1
threshold_grid_Ss =
    ScaledChiSquareSample.(quantile(response.(microarray_Ss), equidistant_grid), 11)


limma_z_cutoffs_bh =
    invert_limma_pvalue.(
        microarray_limma.cutoff,
        threshold_grid_Ss,
        Ref(microarray_limma.prior),
    )
npmle_z_cutoffs_bh =
    invert_limma_pvalue.(
        microarray_npmle.cutoff,
        threshold_grid_Ss,
        Ref(microarray_npmle.prior),
    )
ttest_z_cutoffs_bh = invert_ttest_pvalue.(microarray_t.cutoff, threshold_grid_Ss)


limma_z_cutoffs_005 =
    invert_limma_pvalue.(0.05, threshold_grid_Ss, Ref(microarray_limma.prior))
npmle_z_cutoffs_005 =
    invert_limma_pvalue.(0.05, threshold_grid_Ss, Ref(microarray_npmle.prior))
ttest_z_cutoffs_005 = invert_ttest_pvalue.(0.05, threshold_grid_Ss)

cutoff_matrix =
    [npmle_z_cutoffs_bh limma_z_cutoffs_bh ttest_z_cutoffs_bh npmle_z_cutoffs_005 limma_z_cutoffs_005 ttest_z_cutoffs_005]
cutoff_colors = [:purple :darkorange :grey :purple :darkorange :grey]
cutoff_linestyles = [:solid :solid :solid :dash :dash :dash]

plot!(
    twod_histogram_plot,
    log.(response.(threshold_grid_Ss)),
    cutoff_matrix,
    color = cutoff_colors,
    ylim = (-2.2, 2.2),
    linestyle = cutoff_linestyles,
    label = ["NPMLE (BH)" "Limma (BH)" "t-test (BH)" "NPMLE (unadj.)" "Limma (unadj)" "t-test (unadj.)"],
)

plot!(
    twod_histogram_plot,
    log.(response.(threshold_grid_Ss)),
    -cutoff_matrix,
    color = cutoff_colors,
    linestyle = cutoff_linestyles,
    label = "",
)

#---------------------------------------------------------
# Put all panels together
#---------------------------------------------------------



Plots.plot(histogram_plot, prior_plot, size = (800, 300), layout = (1, 2))

savefig("first_two_microarray_panels.pdf")

Plots.plot(twod_histogram_plot, size = (400, 300))

savefig("microarray_rejection_regions.pdf")
