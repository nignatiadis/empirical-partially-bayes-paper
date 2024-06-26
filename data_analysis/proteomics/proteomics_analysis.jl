using Pkg

dir = @__DIR__
Pkg.activate(joinpath(dir, "..", "..", "simulations"))

using Empirikos
using CSV
using LaTeXStrings
using Statistics
using Plots

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


proteomics = CSV.File(joinpath(dir, "molecular_oncology_proteomics.csv"))


proteomics_Ss =
    ScaledChiSquareSample.(abs2.(proteomics.se_hat), proteomics.residual_dof)
proteomics_mu_hat = proteomics.mu_hat
proteomics_samples = NormalChiSquareSample.(proteomics_mu_hat, proteomics_Ss)

proteomics_npmle = fit(
    Empirikos.EmpiricalPartiallyBayesTTest(
        prior = DiscretePriorClass(),
        solver = Mosek.Optimizer,
        α = 0.05,
        discretize_marginal = false,
    ),
    proteomics_samples,
)

proteomics_limma =
    fit(Empirikos.EmpiricalPartiallyBayesTTest(
        prior = Empirikos.Limma(),
        solver = nothing,
        α = 0.05),
        proteomics_samples)

# Comparison to limma results from R package
#limma_ν_prior = 5.174972897 limma_prior_var = 0.07665211404
proteomics_limma.prior

proteomics_t = fit(SimultaneousTTest(α=0.05), proteomics_samples)

# Rejections by method
        
(
    proteomics_npmle.total_rejections,
    proteomics_limma.total_rejections,
    proteomics_t.total_rejections,
)

#---------------------------------------------------------
# Panel A: Marginal density of variances
#---------------------------------------------------------

extrema(support(proteomics_npmle.prior))

var_grid_refine = 0.01:0.001:1.0
var_grid_refine_samples = ScaledChiSquareSample.(var_grid_refine, 28)

marginal_pdf_limma = pdf.(proteomics_limma.prior, var_grid_refine_samples)
marginal_pdf_npmle = pdf.(proteomics_npmle.prior, var_grid_refine_samples)

histogram_plot = Plots.histogram(
    abs2.(proteomics.se_hat)[abs2.(proteomics.se_hat).<0.7],
    fill = "lightgrey",
    normalize = true,
    fillalpha = 0.5,
    linewidth = 0.3,
    linealpha = 0.5,
    bins = 100,
    xlim = (0, 0.7),
    ylim = (0, 10),
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
    support(proteomics_npmle.prior),
    50 * probs(proteomics_npmle.prior),
    seriestype = :sticks,
    xlim = (0, 0.7),
    legend = :topright,
    label = "NPMLE",
    color = :purple,
    xlabel = L"\sigma^2",
    ylabel = L"g(\sigma^2)",
)

Plots.plot!(
    prior_plot,
    u -> pdf(proteomics_limma.prior, u),
    label = "Limma",
    color = :darkorange,
)


#---------------------------------------------------------
# Panel C: 2D-rejection regions
#---------------------------------------------------------

extrema(log10.(abs2.(proteomics.se_hat)))

log_grid = [0.01; 0.1; 1.0]

twod_histogram_plot = histogram2d(
    log.(abs2.(proteomics.se_hat)),
    proteomics.mu_hat,
    bins = 50,
    c = cgrad(:algae, rev = false, scale = :exp),
    xlabel = L"S_i^2",
    ylabel = L"Z_i",
    ylim = (-3.4, 3.4),
    xticks = (log.(log_grid), string.(log_grid)),
)


extrema(response.(proteomics_Ss))
threshold_grid_Ss = ScaledChiSquareSample.(0.0058:0.002:1.52, 28)


limma_z_cutoffs_bh =
    invert_limma_pvalue.(
        proteomics_limma.cutoff,
        threshold_grid_Ss,
        Ref(proteomics_limma.prior),
    )
npmle_z_cutoffs_bh =
    invert_limma_pvalue.(
        proteomics_npmle.cutoff,
        threshold_grid_Ss,
        Ref(proteomics_npmle.prior),
    )
ttest_z_cutoffs_bh = invert_ttest_pvalue.(proteomics_t.cutoff, threshold_grid_Ss)


limma_z_cutoffs_005 =
    invert_limma_pvalue.(0.05, threshold_grid_Ss, Ref(proteomics_limma.prior))
npmle_z_cutoffs_005 =
    invert_limma_pvalue.(0.05, threshold_grid_Ss, Ref(proteomics_npmle.prior))
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

savefig("first_two_proteomics_panels.pdf")

Plots.plot(twod_histogram_plot, size = (400, 300))

savefig("proteomics_rejection_regions.pdf")
