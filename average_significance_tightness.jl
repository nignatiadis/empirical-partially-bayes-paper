using Cubature
using Empirikos
using LinearAlgebra
using QuadGK
using Expectations

Gmix = DiscreteNonParametric([1.0;0.5], [1/2; 1/2])
Gsingle = DiscreteNonParametric([1.0], [1.0])


pval(z, G) = dot(probs(G), 2*ccdf.( Normal.(0, support(G)), abs(z)))

# sanity check: should be 0.05
pval(1.96, DiscreteNonParametric([1.0], [1.0]))

# sanity check: should be 0.01
quadgk(z-> (pval(z, Gsingle) <= 0.01)*pdf(Normal(),z), -7, +7)

# assume the null has variance 1 (throughout this file)
function limma_type_pval_inner(params; t=0.01, ν=3, G=Gsingle)
    z = params[1]
    s_squared = params[2]
    S_squared = Empirikos.ScaledChiSquareSample(s_squared, ν)
    Gmix = Empirikos.posterior(S_squared, G)
    pdf(Normal(),z)*(pval(z, Gmix) <= t)
end


function limma_level(G, ν, t)
    E = expectation(
    likelihood_distribution(Empirikos.ScaledChiSquareSample(1.0, ν), 1)
    )
    E(s_sq -> quadgk(z->limma_type_pval_inner([z;s_sq]; G=G, ν=ν, t=t), -7, 7)[1])
end


limma_level(Gsingle, 3, 0.01)
limma_level(Gmix, 3, 0.01)
