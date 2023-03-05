using Empirikos
using CSV
using LaTeXStrings
using Statistics
using Plots

dir = @__DIR__
include(joinpath(dir, "..","..", "helpers.jl"))


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

methylation = CSV.File(joinpath(dir, "methylation_naive_vs_act.csv"))

# Code to be added to Empirikos.jl tests
# limma_ν_prior = 3.957179
# limma_prior_var_naive_vs_act = 0.03664045646

methylation_Ss =
    Empirikos.ScaledChiSquareSample.(abs2.(methylation.sigma_hat), 4)
methylation_mu_hat = methylation.mu_hat

methylation_npmle = fit(PartiallyBayesTest(prior=NPMLE, discretize_marginal=true),
    methylation_Ss, methylation_mu_hat)
methylation_limma = fit(PartiallyBayesTest(prior=Limma), methylation_Ss, methylation_mu_hat)
methylation_limma.prior

methylation_t = fit(SimultaneousTTest(), methylation_Ss, methylation_mu_hat)



# rejections by method

(methylation_npmle.total_rejections, methylation_limma.total_rejections, methylation_t.total_rejections)
(methylation_npmle.cutoff, methylation_limma.cutoff, methylation_t.cutoff)

# Identify rejection that was only done by t-test
# and check that it is the one with 2nd smallest sample variance
t_test_idx = findfirst(methylation_t.rj_idx)
methylation_Ss[t_test_idx]
methylation_mu_hat[t_test_idx]
sort(response.(methylation_Ss))[1:2]

#---------------------------------------------------------
# Panel A: Marginal density of variances
#---------------------------------------------------------

extrema(support(methylation_npmle.prior))

var_grid_refine = 0.001:0.001:1.0
var_grid_refine_samples = Empirikos.ScaledChiSquareSample.(var_grid_refine, 4)

marginal_pdf_limma = pdf.(methylation_limma.prior, var_grid_refine_samples)
marginal_pdf_npmle = pdf.(methylation_npmle.prior, var_grid_refine_samples)

histogram_plot = Plots.histogram(
    abs2.(methylation.sigma_hat)[abs2.(methylation.sigma_hat).<0.4],
    fill = "lightgrey",
    normalize = true,
    fillalpha = 0.5,
    linewidth = 0.3,
    linealpha = 0.5,
    bins = 100,
    xlim = (0, 0.4),
    ylim = (0, 18),
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
    support(methylation_npmle.prior),
    100 * probs(methylation_npmle.prior),
    seriestype = :sticks,
    xlim = (0, 0.4),
    legend = :topright,
    label = "NPMLE",
    color = :purple,
    xlabel = L"\sigma^2",
    ylabel = L"g(\sigma^2)",
)

Plots.plot!(
    prior_plot,
    u -> pdf(methylation_limma.prior, u),
    label = "Limma",
    color = :darkorange,
)



#---------------------------------------------------------
# Panel C: 2D-rejection regions
#---------------------------------------------------------

extrema(log.(abs2.(methylation.sigma_hat)))


log_grid = [0.0001; 0.001; 0.01; 0.1; 1; 10]

twod_histogram_plot = histogram2d(
    log.(abs2.(methylation.sigma_hat)),
    methylation_mu_hat,
    bins = 50,
    c = cgrad(:algae, rev = false, scale = :exp),
    xlabel = L"S_i^2",
    ylabel = L"Z_i",
    xticks = (log.(log_grid), string.(log_grid)),
)

equidistant_grid = 0.0:0.001:1
extrema(response.(methylation_Ss))
threshold_grid_Ss = ScaledChiSquareSample.(quantile(response.(methylation_Ss), equidistant_grid), 4)


limma_z_cutoffs_bh = invert_limma_pvalue.(methylation_limma.cutoff, threshold_grid_Ss, Ref(methylation_limma.prior))
npmle_z_cutoffs_bh = invert_limma_pvalue.(methylation_npmle.cutoff, threshold_grid_Ss, Ref(methylation_npmle.prior))
ttest_z_cutoffs_bh = invert_ttest_pvalue.(methylation_t.cutoff, threshold_grid_Ss)


limma_z_cutoffs_005 = invert_limma_pvalue.(0.05, threshold_grid_Ss, Ref(methylation_limma.prior))
npmle_z_cutoffs_005 = invert_limma_pvalue.(0.05, threshold_grid_Ss, Ref(methylation_npmle.prior))
ttest_z_cutoffs_005 = invert_ttest_pvalue.(0.05, threshold_grid_Ss)

cutoff_matrix = [npmle_z_cutoffs_bh limma_z_cutoffs_bh ttest_z_cutoffs_bh npmle_z_cutoffs_005  limma_z_cutoffs_005 ttest_z_cutoffs_005]
cutoff_colors = [:purple :darkorange :grey :purple :darkorange :grey]
cutoff_linestyles = [:solid :solid :solid :dash :dash :dash]

plot!(twod_histogram_plot,  log.(response.(threshold_grid_Ss)),cutoff_matrix,
    color = cutoff_colors,
    ylim = (-6, 6),
    linestyle= cutoff_linestyles,
    label = ["NPMLE (BH)" "Limma (BH)" "t-test (BH)" "NPMLE (unadj.)" "Limma (unadj)" "t-test (unadj.)"])

plot!(twod_histogram_plot, log.(response.(threshold_grid_Ss)), -cutoff_matrix,
    color = cutoff_colors,
    linestyle= cutoff_linestyles,
    label = "")

#---------------------------------------------------------
# Put all panels together
#---------------------------------------------------------


Plots.plot(
    histogram_plot,
    prior_plot,
    size = (800, 300),
    layout = (1,2)
)

savefig("first_two_methylation_panels.pdf")

Plots.plot(
    twod_histogram_plot,
    size = (400, 300),
)

savefig("methylation_rejection_regions.pdf")

# Stratified p-value histograms


ten_percent_quantile_low = quantile(response.(methylation_Ss), 0.10)
ten_percent_quantile_high = quantile(response.(methylation_Ss), 0.90)

keep_idx_low = response.(methylation_Ss) .<= ten_percent_quantile
keep_idx_high = response.(methylation_Ss) .>= ten_percent_quantile_high
keep_idx_all = fill(true, length(keep_idx_low))

hist_params = (
normalize = true,
fillalpha = 0.2,
linewidth = 0.3,
bins = 20,
xlim = (0, 1),
ylim = (0.001,7),
xlabel = L"P_i",
ylabel = "Density",
legend = :topright,
titlefontsize=10
)


ttest_histograms = plot(
    histogram(methylation_t.pvalue[keep_idx_all];
        fillcolor=:grey,
        title = L"\textrm{ all }\, S_i^2",
        label = "t-test",
        hist_params...),
    histogram(methylation_t.pvalue[keep_idx_low];
        fillcolor=:grey,
        title = L"S_i^2\, \textrm{ in  bottom }\, 10\%",
        label = "t-test",
        hist_params...),
    histogram(methylation_t.pvalue[keep_idx_high];
        fillcolor=:grey,
        title = L"S_i^2\, \textrm{ in  top }\, 10\%",
        label = "t-test",
        hist_params...),
    histogram(methylation_npmle.pvalue[keep_idx_all];
        fillcolor=:purple,
        title = L"\textrm{ all }\, S_i^2",
        label = "NPMLE",
        hist_params...),
    histogram(methylation_npmle.pvalue[keep_idx_low];
        fillcolor=:purple,
        title = L"S_i^2\, \textrm{ in  bottom }\, 10\%",
        label = "NPMLE",
        hist_params...),
    histogram(methylation_npmle.pvalue[keep_idx_high];
        fillcolor=:purple,
        title = L"S_i^2\, \textrm{ in  top }\, 10\%",
        label = "NPMLE",
        hist_params...),
    bottom_margin = 9*Plots.mm,
    layout = (2,3),
    size = (1000, 500),
    thickness_scaling=1.6
)

savefig("methylation_conditional_pvalue_histograms.pdf")
