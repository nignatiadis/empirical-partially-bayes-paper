using Empirikos
using MultipleTesting
using MosekTools
using Roots
import StatsBase:fit

Base.@kwdef struct SimultaneousTTest
    multiple_test =  BenjaminiHochberg()
    α::Float64 = 0.05
end

function fit(test::SimultaneousTTest, samples)
    mu_hat = getproperty.(samples, :Z)
    Ss = ScaledChiSquareSample.(samples)

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
