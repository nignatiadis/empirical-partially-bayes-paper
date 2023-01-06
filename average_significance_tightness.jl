using Cubature
using Empirikos
using LinearAlgebra
using QuadGK
using Expectations

Gmix = DiscreteNonParametric([1.0;0.5], [1/2; 1/2])
Gsingle = DiscreteNonParametric([1.0], [1.0])


pval(z, G) = dot(probs(G), 2*ccdf.( Normal.(0, support(G)), abs(z)))

pval(1.96, DiscreteNonParametric([1.0], [1.0]))

quadgk(z-> (pval(z, Gsingle) <= 0.01)*pdf(Normal(),z), -7, +7)

Empirikos.limma_pvalue

# assume the null has variance 1
function limma_type_pval_inner(params; t=0.01, ν=3, G=Gsingle)
    z = params[1]
    s_squared = params[2]
    S_squared = Empirikos.ScaledChiSquareSample(s_squared, ν)
    Gmix = Empirikos.posterior(S_squared, G)
    pdf(Normal(),z)*(pval(z, Gmix) <= t)
end







function limma_level(G, ν, t)
    E = expectation(
    likelihood_distribution(Empirikos.ScaledChiSquareSample(1.0, ν),1)
    )
    E(s_sq -> quadgk(z->limma_type_pval_inner([z;s_sq]; G=G, ν=ν, t=t), -7, 7)[1])
end



νs = Base.OneTo(60)
_level = limma_level.(Ref(Gmix), νs, 0.01)



using Plots
pgfplotsx()

theme(
    :default;
    background_color_legend = :transparent,
    foreground_color_legend = :transparent,
    grid = nothing,
    frame = :box,
    legendfonthalign = :left,
    thickness_scaling = 1.1
)

plot(νs, _level, seriestype=:scatter, label="",  xguide="ν")
plot!([1; 60], [0.01; 0.01], linestyle=:dash,
    label="",
    ylabel="Average significance", xlim=(1,60),
    xticks=1:6:60, ylim=(0,0.0205))
savefig("avg_significance_violation.pdf")

function limma_type_pval(params; ν=3, G=Gsingle)
    z = params[1]
    s_squared = params[2]
    S_squared = Empirikos.ScaledChiSquareSample(s_squared, ν)
    lik = likelihood_distribution(S_squared, 1)
    Gmix = Empirikos.posterior(S_squared, G)
    pdf(Normal(),z)*pdf(lik, s_squared)*(pval(z, Gmix) <= 0.01)
end

using Cubature

mean(Empirikos.posterior(Empirikos.ScaledChiSquareSample(1, 3), Gmix))
limma_type_pval([1.96;1.0];G=Gsingle)



hcubature(z->limma_type_pval(z;G=Gmix), [-7; 1e-5], [+7;20])
