using Empirikos
using RangeHelpers
using MultipleTesting
using MosekTools
using Roots
import StatsBase:fit

struct Limma end

Base.@kwdef struct PartiallyBayesTest
    prior
    multiple_test =  BenjaminiHochberg()
    α::Float64 = 0.05
    prior_grid = 300
    discretize_marginal = false
    solver = Mosek.Optimizer
end

function fit(test::PartiallyBayesTest, Ss, mu_hat)

    prior = fit_prior(test, test.prior, Ss)
    multiple_test = test.multiple_test
    α = test.α

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
    method = test,
    prior = prior,
    pvalue = pvalues,
    cutoff = cutoff,
    adjp = adjp,
    rj_idx = rj_idx,
    total_rejections = total_rejections
    )
end

function fit_prior(test::PartiallyBayesTest, prior::Type{NPMLE}, Ss)
    discretize_marginal = test.discretize_marginal
    prior_grid = test.prior_grid
    solver = test.solver

    a_min = quantile(response.(Ss), 0.01)
    a_max = maximum(response.(Ss))

    grid = exp.(range(start=log(a_min), stop=log(a_max), length=prior_grid))
    _prior = DiscretePriorClass(grid)
    _npmle = NPMLE(_prior, solver)

    if discretize_marginal
        disc = interval_discretizer(grid)
        Ss_summary = summarize(disc.(Ss))
        npmle_prior = fit(_npmle, Ss_summary)
    else
        npmle_prior = fit(_npmle, Ss)
    end
    npmle_prior.prior
end


function fit_prior(test::PartiallyBayesTest, prior::Type{Limma}, Ss)
    Empirikos.fit_limma(Ss)
end

function fit_prior(test::PartiallyBayesTest, prior::Distribution, Ss)
    prior
end

Base.@kwdef struct SimultaneousTTest
    multiple_test =  BenjaminiHochberg()
    α::Float64 = 0.05
end
function fit(test::SimultaneousTTest, Ss, mu_hat)
    multiple_test = test.multiple_test
    α = test.α
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
    largest_se = 100* max( sqrt(response(S)), quantile(prior, 0.999))
    find_zero(pval_fun, (0, largest_se))
end


function invert_ttest_pvalue(threshold, S)
    quantile(TDist(dof(S)), 1 - threshold/2) * sqrt(response(S))
end
