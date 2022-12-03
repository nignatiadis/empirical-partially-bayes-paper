using Empirikos
using RangeHelpers
using MultipleTesting
using MosekTools
using Roots

function npmle_limma_analysis(Ss, mu_hat;
    multiple_test = BenjaminiHochberg(),
    α = 0.05,
    prior_grid=500,
    discretize_marginal = false)

    a_min = quantile(response.(Ss), 0.01)
    a_max = maximum(response.(Ss))

    grid = exp.(range(start=log(a_min), stop=log(a_max), length=prior_grid))
    _prior = DiscretePriorClass(grid)
    _npmle = NPMLE(_prior, Mosek.Optimizer)

    if discretize_marginal
        disc = interval_discretizer(grid)
        Ss_summary = summarize(disc.(methylation_Ss))
        npmle_prior = fit(_npmle, Ss_summary)
    else
        npmle_prior = fit(_npmle, Ss)
    end

    npmle_pvalues =   Empirikos.limma_pvalue.(mu_hat, Ss, npmle_prior.prior)
    npmle_adjp = adjust(npmle_pvalues, multiple_test)
    rj_idx = npmle_adjp .<= α
    total_rejections = sum(rj_idx)
    if iszero(total_rejections)
        cutoff = zero(Float64)
    else
        cutoff = maximum(npmle_pvalues[rj_idx])
    end

    (
    prior = npmle_prior.prior,
    pvalue = npmle_pvalues,
    cutoff = cutoff,
    adjp = npmle_adjp,
    rj_idx = rj_idx,
    total_rejections = total_rejections
    )

end


function parametric_limma_analysis(Ss, mu_hat;
    multiple_test = BenjaminiHochberg(),
    α = 0.05)

    prior = Empirikos.fit_limma(Ss)
    pvalues =   Empirikos.limma_pvalue.(mu_hat, Ss, prior)
    adjp = adjust(pvalues, multiple_test)
    rj_idx = adjp .<= α
    total_rejections = sum(rj_idx)
    if iszero(total_rejections)
        cutoff = zero(Float64)
    else
        cutoff = maximum(pvalues[rj_idx])
    end

    (
    prior = prior,
    pvalue = pvalues,
    cutoff = cutoff,
    adjp = adjp,
    rj_idx = rj_idx,
    total_rejections = total_rejections
    )

end

function ttest_analysis(Ss, mu_hat;
    multiple_test = BenjaminiHochberg(),
    α = 0.05)
    pvalues = 2*ccdf.(TDist.(dof.(Ss)), abs.(mu_hat) ./ sqrt.(response.(Ss)))
    adjp = adjust(pvalues, multiple_test)
    rj_idx = adjp .<= α
    total_rejections = sum(rj_idx)
    if iszero(total_rejections)
        cutoff = zero(Float64)
    else
        cutoff = maximum(pvalues[rj_idx])
    end

    (
    pvalue = pvalues,
    cutoff = cutoff,
    adjp = adjp,
    rj_idx = rj_idx,
    total_rejections = total_rejections
    )

end

function invert_limma_pvalue(threshold, S, prior)
    pval_fun(z) = Empirikos.limma_pvalue(z, S, prior) - threshold
    largest_se = 100* sqrt(response(S))
    find_zero(pval_fun, (0, largest_se))
end


function invert_ttest_pvalue(threshold, S)
    quantile(TDist(dof(S)), 1 - threshold/2) * sqrt(response(S))
end
