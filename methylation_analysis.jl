using Empirikos
using CSV
using MultipleTesting
using MosekTools
using LaTeXStrings
using Statistics
using Plots
using Test


pgfplotsx()
theme(
    :default;
    background_color_legend = :transparent,
    foreground_color_legend = :transparent,
    grid = nothing,
    frame = :box,
    legendfonthalign = :left,
    thickness_scaling = 1.2,
    size = (420, 330),
)

#Plots.gr_cbar_width[] = 0.005

methylation = CSV.File("methylation_naive_vs_act.csv")
limma_ν_prior = 3.957179
limma_prior_var_naive_vs_rtreg = 0.0492029
limma_prior_var_naive_vs_act = 0.03664045646

@test sum(adjust(methylation.limma_pvalue, BenjaminiHochberg()) .<= 0.05) == 627


limma_prior = Empirikos.InverseScaledChiSquare(limma_prior_var_naive_vs_act, limma_ν_prior)

methylation_Ss = Empirikos.ScaledChiSquareSample.(abs2.(methylation.sigma_hat), 4)
mu_hat = methylation.mu_hat
limma_pvalues =   Empirikos.limma_pvalue.(mu_hat, methylation_Ss, limma_prior)
@test limma_pvalues ≈ methylation.limma_pvalue

bh_limma_adjp = adjust(limma_pvalues, BenjaminiHochberg())
α = 0.05
bh_limma_rj_idx = bh_limma_adjp .<= α
@test sum(bh_limma_rj_idx) == 627

quantile(log.( abs2.(methylation.sigma_hat) ), [0.01;1])

extrema(log.( abs2.(methylation.sigma_hat) ))
sigma_sq_grid = exp.(-6.2:0.05:2.6)
sigma_disc_grid = exp.(-6.2:0.05:2.6)
disc = interval_discretizer(sigma_disc_grid)
_prior = DiscretePriorClass(sigma_sq_grid)
Ss_summary = summarize(disc.(methylation_Ss))
npmle_prior = fit(NPMLE(_prior, Mosek.Optimizer), Ss_summary)
npmle_pvalues =   Empirikos.limma_pvalue.(mu_hat, methylation_Ss, npmle_prior.prior)
bh_npmle_adjp = adjust(npmle_pvalues, BenjaminiHochberg())
bh_npmle_rj_idx = bh_npmle_adjp .<= α

bh_cutoff_limma = maximum(limma_pvalues[bh_limma_rj_idx])
bh_cutoff_npmle = maximum(npmle_pvalues[bh_npmle_rj_idx])

using Roots

S = Empirikos.ScaledChiSquareSample(0.1, 4)

function invert_limma_pvalue(threshold, S, prior)
    pval_fun(z) = Empirikos.limma_pvalue(z, S, prior) - threshold
    largest_se = 100* sqrt(response(S))
    find_zero(pval_fun, (0, largest_se))
end

t_test_cutoff = quantile(TDist(4), 1-0.025) * sqrt.(response.(SS_grid))

Ss_grid = Empirikos.ScaledChiSquareSample.(sigma_sq_grid, 4)
limma_z_cutoffs = invert_limma_pvalue.(0.05, Ss_grid, Ref(limma_prior))

npmle_z_cutoffs = invert_limma_pvalue.(0.05, Ss_grid, Ref(npmle_prior.prior))

var_grid_refine_samples = Empirikos.ScaledChiSquareSample.(var_grid_refine, 4)

plot(log.(response.(Ss_grid)), [limma_z_cutoffs npmle_z_cutoffs t_test_cutoff], label=["limma" "npmle" "t-test"])


invert_limma_pvalue()
gr()
plot(0:0.01:10, x-> Empirikos.limma_pvalue(x,  maximum(methylation_Ss), npmle_prior.prior))

plot!(0:0.01:10, x-> Empirikos.limma_pvalue(x,  Empirikos.ScaledChiSquareSample(0.22, 4), npmle_prior.prior),
lable="quant")


sum(bh_npmle_rj_idx)
#1090, 1020
sum(bh_limma_rj_idx .* bh_npmle_rj_idx)

ttest_pvalues = 2*ccdf(TDist(4), abs.(mu_hat ./ methylation.sigma_hat) )

findall(adjust(ttest_pvalues, BenjaminiHochberg()) .<= α)

abs2(methylation.sigma_hat[335269])


pgfplotsx()
twod_histogram_plot = histogram2d(log.(methylation.sigma_hat),methylation.mu_hat,
 bins=50, c=cgrad(:algae, rev=false, scale=:exp),
 xlabel=L"S_i^2", ylabel=L"Z_i",
 tex_output_standalone = true)

Plots.plot!(twod_histogram_plot,
 log.(methylation.sigma_hat)[bh_limma_rj_idx .* bh_npmle_rj_idx],
 methylation.mu_hat[bh_limma_rj_idx .* bh_npmle_rj_idx], seriestype=:scatter,
 markershape = :xcross,
 xticks =( -5:1:1, string.(round.(exp.(-5:1:1);digits=3))),
 label="Limma & NPMLE",
 markerstrokecolor=:orange, legend=:topleft)

Plots.plot!(twod_histogram_plot,
 log.(methylation.sigma_hat)[(.!bh_limma_rj_idx) .* bh_npmle_rj_idx],
 methylation.mu_hat[(.!bh_limma_rj_idx) .* bh_npmle_rj_idx], seriestype=:scatter,
 markershape = :cross, legendfontsize=10,
 markerstrokecolor=:purple,
 label="NPMLE only")

Plots.plot!(twod_histogram_plot,
 [log(methylation.sigma_hat[335269])],
 [methylation.mu_hat[335269]],
 seriestype = :scatter,
 markershape=:utriangle,markercolor=:grey,
 label="t-test")

savefig("2dhistogram.pdf")

var_grid_refine =  1e-5:0.001:0.5
var_grid_refine_samples = Empirikos.ScaledChiSquareSample.(var_grid_refine, 4)

marginal_pdf_limma= pdf.(limma_prior, var_grid_refine_samples)
marginal_pdf_npmle = pdf.(npmle_prior.prior, var_grid_refine_samples)

histogram_plot = Plots.histogram(
        abs2.(methylation.sigma_hat)[methylation.sigma_hat .< sqrt(0.4) ],
        fill="lightgrey",
        normalize=true, fillalpha=0.5, linewidth=0.3,
        bins=100, xlim=(0,0.4), ylim=(0,18), label="Histogram",
        xlabel=L"S_i^2", ylabel="Density", legend=:topright)

Plots.plot!(histogram_plot, var_grid_refine,
        [marginal_pdf_npmle marginal_pdf_limma], label=["NPMLE" "Limma"],
        color = [:purple :darkorange],
        linewidth=0.7)



prior_plot = Plots.plot(support(npmle_prior.prior), 100 * probs(npmle_prior.prior),
     seriestype=:sticks,
     xlim=(0, 0.5),
     legend = :topright, label="NPMLE",
     color = :purple,
     xlabel=L"\sigma^2",
     ylabel=L"g(\sigma^2)")
#limma_f = u->pdf(limma_prior,u^2 )*2*u
Plots.plot!(prior_plot,
    u->pdf(limma_prior,u ), label="Inverse Gamma", color=:darkorange)
savefig("methylation_estimated_priors.pdf")

pvalue_hist = histogram(npmle_pvalues, fill=:blue,
    normalize=true, fillalpha=0.2, linewidth=0.3,
    bins=20, xlim=(0,1), ylim=(0,2.7), label="NPMLE",
    xlabel=L"P_i", ylabel="Density", legend=:topright)
histogram!(pvalue_hist, limma_pvalues, fill=:darkorange,
    normalize=true, fillalpha=0.2, linewidth=0.3,
    bins=20, xlim=(0,1),  label="Limma",
    xlabel=L"P_i", ylabel="Density", legend=:topright)

savefig("methylation_pvalue_histograms.pdf")


Plots.plot(histogram_plot, prior_plot, pvalue_hist, twod_histogram_plot,
             size=(1000,850), margin = 5Plots.mm,
             layout = (2,2), title=["a)" "b)" "c)" "d)"])

savefig("methylation_plot.pdf")



















sum( probs(npmle_prior.prior) .> 1e-4)

loglikelihood(methylation_Ss, npmle_prior.prior)

loglikelihood(methylation_Ss, limma_prior)



loglikelihood(methylation_Ss, npmle_prior.prior)






































using StatsBase
ecdf_limma_ps = ecdf(limma_pvalues)
ecdf_npmle_ps = ecdf(npmle_pvalues)
ps = 0.00:0.01:1.0

plot(ps, [ecdf_npmle_ps.(ps) ecdf_limma_ps.(ps)],
    linestyle=[:dot :solid], linewidth=[1.3 1.1],
    linealpha=0.7, color=[:blue :darkorange],
    label=["NPMLE" "Inverse Gamma"], legend=:bottomright,
    xlabel=L"P_i", ylabel="Distribution")

plot(ps, ecdf_npmle_ps.(ps),
    color=:blue, markershape=:circle,
    label="NPMLE", legend=:bottomright, markersize=0.3,
    xlabel=L"P_i", ylabel="Distribution")

savefig("methylation_cdfs.pdf")

plot(1:10,1:10)

filter_idx =  methylation.sigma_hat .< 1
histogram(abs2.(methylation.sigma_hat[filter_idx]), normalize=true, fill="white",
    bins=100)







mean(limma_prior)
std(limma_prior)
mean(npmle_prior.prior)
median(InverseGamma(limma_prior))

median(npmle_prior.prior)
#std(npmle_prior.prior)



plot(support(npmle_prior.prior), 1000 * probs(npmle_prior.prior), seriestype=:scatter,
     xlim=(0, 1), markersize = 2, markerstrokewidth = 0, markershape=:x)

plot(support(npmle_prior.prior), 1000 * probs(npmle_prior.prior), seriestype=:line,
     xlim=(0, 1), markersize = 2, markerstrokewidth = 0)
plot!(u->pdf(limma_prior,u ))

plot(sigma_sq_grid , u->pdf(npmle_prior.prior,u ))


plot(sqrt.(support(npmle_prior.prior)), 50 * probs(npmle_prior.prior), seriestype=:sticks,
     xlim=(0, 1))
#limma_f ./ sum(limma_f, sqrt.(support(npmle_prior.prior)))

var_grid = 1e-5:0.005:7
var_grid_samples = Empirikos.ScaledChiSquareSample.(var_grid, 4)
postmeans = PosteriorMean.(var_grid_samples)
npmle_postmeans = postmeans.(npmle_prior.prior)
limma_postmeans = postmeans.(limma_prior)

plot(var_grid, [npmle_postmeans limma_postmeans], label=["npmle" "limma"])

pvals_npmle_fun = Empirikos.limma_pvalue.(1.0, var_grid_samples, npmle_prior.prior)
pvals_limma_fun = Empirikos.limma_pvalue.(1.0, var_grid_samples, limma_prior)


var_grid = 0.001:0.001:0.5
var_grid_samples = Empirikos.ScaledChiSquareSample.(var_grid, 4)

pvals_npmle_fun = Empirikos.limma_pvalue.(1.0 .* sqrt.(var_grid), var_grid_samples, npmle_prior.prior)
pvals_limma_fun = Empirikos.limma_pvalue.(1.0 .* sqrt.(var_grid), var_grid_samples, limma_prior)

pvals_npmle_fun2 = Empirikos.limma_pvalue.(2.0 .* sqrt.(var_grid), var_grid_samples, npmle_prior.prior)
pvals_limma_fun2 = Empirikos.limma_pvalue.(2.0 .* sqrt.(var_grid), var_grid_samples, limma_prior)
pvals_limma_fun3 = Empirikos.limma_pvalue.(3.0 .* sqrt.(var_grid), var_grid_samples, limma_prior)
pvals_npmle_fun3 = Empirikos.limma_pvalue.(3.0 .* sqrt.(var_grid), var_grid_samples, npmle_prior.prior)

plot(var_grid, [pvals_npmle_fun pvals_limma_fun], label=["npmle 1" "limma 1"])
plot(var_grid, [pvals_npmle_fun2 pvals_limma_fun2], label=["npmle 2" "limma 2"])
plot!(var_grid, )
plot(var_grid, [pvals_npmle_fun3 pvals_limma_fun3], label=["npmle 3" "limma 3"])


effects = 0.0:0.1:3.0

pvals_npmle_var005 = Empirikos.limma_pvalue.(effects.* sqrt.(0.3), Empirikos.ScaledChiSquareSample(0.3, 4), npmle_prior.prior)
pvals_limma_var005 = Empirikos.limma_pvalue.(effects.* sqrt.(0.3), Empirikos.ScaledChiSquareSample(0.3, 4), limma_prior)

plot(effects, [pvals_npmle_var005 pvals_limma_var005], label=["npmle 1" "limma 1"])



var_grid_refine =  1e-5:0.001:1
var_grid_refine_samples = Empirikos.ScaledChiSquareSample.(var_grid_refine, 4)

marginal_pdf_limma= pdf.(limma_prior, var_grid_refine_samples)
marginal_pdf_npmle = pdf.(npmle_prior.prior, var_grid_refine_samples)


filter_idx =  methylation.sigma_hat .< 1
histogram(abs2.(methylation.sigma_hat[filter_idx]), normalize=true, fill="white",
    bins=100)

plot!(var_grid_refine, [marginal_pdf_npmle marginal_pdf_limma], label=["npmle" "limma"],
linesize=3)

plot(1:10,1:10)

savefig("methylation_pvalue_vs_pvalue.pdf")

using StatsPlots





findmin(npmle_pvalues)
findmin(limma_pvalues)

methylation[202561]

cor( limma_pvalues, npmle_pvalues)
marginalhist(limma_pvalues, npmle_pvalues)
plot!([0;1]; [0,1])

plot(-log10.(limma_pvalues), -log10.(npmle_pvalues), markersize=1,
    markerstrokewidth = 0, seriestype=:scatter,
    aspect_ratio=:equal)
plot!([-8;8]; [-8,8])
σ²_prior = 4.0
ν_lik = 4
ν_prior = 3
n = 10_000
π1 = 0.05
v0 = 2


Z = Empirikos.ScaledChiSquareSample(2.0, ν_lik)
G = Empirikos.InverseScaledChiSquare(σ²_prior, ν_prior)


pvalue = Empirikos.limma_pvalue(2*1.96, Z, G)
_post = Empirikos.posterior(Z, G)


pvalue2_fun(β) = quadgk(u-> 2*ccdf(Normal(), β/sqrt(u))*pdf(_post, u), 0 , Inf)[1]
pvalue2_fun(2*1.96)

Zplus1 = Empirikos.ScaledChiSquareSample(2.0, ν_lik+1)
marginal = Empirikos.marginalize(Z, G)

marginal_plus1 = Empirikos.marginalize(Zplus1, G)
ssq = Zs.Z
int_lb = (abs2(2*1.96) + ν_lik * ssq)/(1+ν_lik)
pvalue3 = quadgk(t -> _const * t^(-(ν_lik-1)/2)*pdf(marginal_plus1, t) / (2*sqrt((ν_lik+1)t - ν_lik*ssq)), int_lb  , Inf)[1]

quadgk(x->pdf(marginal, x), 0, Inf)[1]
quadgk(x->pdf(marginal_plus1, x), 0, Inf)[1]

using SpecialFunctions
_const = ((1+1/ν_lik)^(-ν_lik/2)*sqrt(ν_lik +1)/ sqrt(pi)/pdf(marginal, ssq)) * ssq^(ν_lik/2-1) * gamma((ν_lik+1)/2) / gamma(ν_lik/2)

function pvalue3_fun(β)
    int_lb = (abs2(β) + ν_lik * ssq)/(1+ν_lik)
    quadgk(t -> _const * t^(-(ν_lik-1)/2)*pdf(marginal_plus1, t) / (sqrt((ν_lik+1)t - ν_lik*ssq)), int_lb  , Inf)[1]
end
pvalue3_fun(2*1.96)/pvalue2_fun(2*1.96)
pvalue3_fun(3*1.96)/pvalue2_fun(3*1.96)
pvalue3_fun(0.5)/pvalue2_fun(0.5)
pvalue3_fun(0.1)/pvalue2_fun(0.)


pvalue3/pvalue2
pvalue3_fun(0.5),pvalue2_fun(0.5)
pvalue3_fun(0.1),pvalue2_fun(0.1)


using QuadGK



loglikelihood(Zs, Empirikos.InverseScaledChiSquare(σ²_prior, ν_prior))




function single_sim()
    #σs_squared = rand(InverseGamma(Empirikos.InverseScaledChiSquare(σ²_prior, ν_prior)), n)
    σs_squared = fill(σ²_prior, n)
    β = fill(0.0, n)
    alt_idxs = Base.OneTo(floor(Int, π1*n))
    β[alt_idxs] .= 6.0
    #β[alt_idxs] = rand(Normal(), length(alt_idxs)) .* sqrt(v0) .* sqrt.(σs_squared[alt_idxs])
    Ss_squared = σs_squared ./ ν_lik .* rand(Chisq(ν_lik), n)
    Ss_squared = Empirikos.ScaledChiSquareSample.(Ss_squared, ν_lik)
    βs_hat = rand(Normal(), n) .* sqrt.(σs_squared) .+  β
    (βs_hat = βs_hat, Ss_squared = Ss_squared, βs = β, σs_squared = σs_squared)
end

Random.seed!(1)

function single_sim_eval()
    sim_res = single_sim()

    true_prior = Empirikos.InverseScaledChiSquare(σ²_prior, ν_prior)
    sigma_sq_grid = abs2.(0.0001:0.1:20)
    _prior = DiscretePriorClass(sigma_sq_grid)
    npmle_prior = fit(NPMLE(_prior, Mosek.Optimizer), sim_res.Ss_squared)

    npmle_pvalues = Empirikos.limma_pvalue.(sim_res.βs_hat, sim_res.Ss_squared, npmle_prior.prior)
    npmle_adjusted_pvalues = adjust(npmle_pvalues, BenjaminiHochberg())
    npmle_rjs = npmle_adjusted_pvalues .<= 0.1

    npmle_discoveries = max(sum( npmle_rjs), 1)
    npmle_true_discoveries = sum( npmle_rjs .&  (sim_res.βs .!= 0 ) )
    npmle_FDP = sum( npmle_rjs .*  (sim_res.βs .== 0 ) ) / npmle_discoveries

    oracle_pvalues = pvalue.(OneSampleZTest.( sim_res.βs_hat, sqrt(σ²_prior), 1))
    oracle_rjs = adjust(oracle_pvalues, BenjaminiHochberg()).<= 0.1
    oracle_discoveries = max(sum( oracle_rjs), 1)
    oracle_true_discoveries = sum( oracle_rjs .&  (sim_res.βs .!= 0 ) )
    oracle_FDP = sum( oracle_rjs .*  (sim_res.βs .== 0 ) ) / oracle_discoveries

    t_pvalues = pvalue.(OneSampleTTest.( sim_res.βs_hat,
    sqrt(ν_lik + 1) .* sqrt.(response.(sim_res.Ss_squared)), ν_lik + 1))
    t_rjs = adjust(t_pvalues, BenjaminiHochberg()).<= 0.1
    t_discoveries = max(sum( t_rjs), 1)
    t_true_discoveries = sum( t_rjs .&  (sim_res.βs .!= 0 ) )
    t_FDP = sum( t_rjs .*  (sim_res.βs .== 0 ) ) / t_discoveries

    return(
            t_true_discoveries = t_true_discoveries,
            t_discoveries = t_discoveries,
            t_FDP = t_FDP,
            oracle_true_discoveries = oracle_true_discoveries,
            oracle_discoveries = oracle_discoveries,
            oracle_FDP = oracle_FDP,
            npmle_true_discoveries = npmle_true_discoveries,
            npmle_discoveries = npmle_discoveries,
            npmle_FDP = npmle_FDP
    )
end

res = [ single_sim_eval() for i in 1:25]

npmle_FDR = mean(getfield.(res, :npmle_FDP))
#t_FDR = mean(getfield.(res, :t_FDP))
oracle_FDR = mean(getfield.(res, :oracle_FDP))

npmle_power = mean(getfield.(res, :npmle_true_discoveries))
#t_power = mean.(getfield.(res, :t_true_discoveries))
oracle_power = mean(getfield.(res, :oracle_true_discoveries))


histogram(tmp.t_pvalues[501:n], bins=0:0.02:1.0)
histogram(tmp.oracle_pvalues[501:n], bins=0:0.02:1.0)

histogram(tmp.npmle_pvalues[501:n], bins=0:0.02:1.0)

using Plots
using StatsPlots



hat_β = β .+ rand(Normal(), n) .* sqrt.(σs_squared)

Zs = Empirikos.ScaledChiSquareSample.(Zs, ν_lik)
using Optim




loglikelihood(Zs, Empirikos.InverseScaledChiSquare(σ²_prior, ν_prior))

ll_fun(v) = -loglikelihood(Zs, Empirikos.InverseScaledChiSquare(v[1], v[2]))
dfc = TwiceDifferentiableConstraints([0.0, 0.0], [Inf, Inf])


res = optimize(ll_fun, dfc, [1.0; 1.0], IPNewton())

estimated_prior = Empirikos.InverseScaledChiSquare(Optim.minimizer(res)...)


using StatsPlots
parametric_pvalues = Empirikos.limma_pvalue.(hat_β, Zs, estimated_prior)

true_pvalues = Empirikos.limma_pvalue.(hat_β, Zs, true_prior)

histogram(parametric_pvalues, bins=0:0.05:1.0)
histogram(true_pvalues, bins=0:0.05:1.0)


plot( -log10.(true_pvalues), -log10.(parametric_pvalues))

Empirikos.posterior(Zs[1], estimated_prior)

quantile(InverseGamma(true_prior), 0.99999)



using Hypatia


npmle_prior = fit(NPMLE(_prior, Hypatia.Optimizer), Zs)

maximum( support(npmle_prior.prior)[ probs(npmle_prior.prior) .> 1e-6])
maximum(σs_squared)


plot( support(npmle_prior.prior),  probs(npmle_prior.prior) , seriestype=:sticks, xlim=(0,100))
plot!( 0:0.01:100, estimated_prior )
loglikelihood(Zs, npmle_prior.prior)
loglikelihood(Zs, estimated_prior)

loglikelihood(Zs, Empirikos.InverseScaledChiSquare(σ²_prior, ν_prior))



quantile(σs_squared)


using Plots

npmle_pvalues = Empirikos.limma_pvalue.(hat_β, Zs, npmle_prior.prior)
histogram(npmle_pvalues, bins=0:0.02:1.0)
histogram(oracle_pvalues, bins=0:0.02:1.0)
histogram(t_pvalues, bins=0:0.02:1.0)



plot( -log10.(true_pvalues), -log10.(npmle_pvalues), seriestype=:scatter, label="")
plot!(0:4, 0:4, linestyle=:dash, label="")
sum(npmle_pvalues .<= 0.05/n * 300)
sum(true_pvalues .<= 0.05/n *300)
