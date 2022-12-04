include("../helpers.jl")
include("simulation_helpers.jl")

using DataFrames
#σ²_prior = 4.0
#ν_lik = 4
#ν_prior = 3
#n = 15_000
#π1 = 0.02
#v0 = 2*3

#νs = [2; 4; 8; 16; 32; 64; 128]


#ground_truth_prior = Empirikos.InverseScaledChiSquare(σ²_prior, ν_prior)
#ground_truth_prior = Dirac(1.0)
#ground_truth_prior = DiscreteNonParametric([1.0; 2.0], [0.5; 0.5])
#ground_truth_prior = DiscreteNonParametric([1.0; 10.0], [0.5; 0.5])

Random.seed!(1)


method_list = (t_standard_BH =  SimultaneousTTest(α=0.1),
    t_limma_BH = PartiallyBayesTest(prior=Limma, α=0.1),
    t_npmle_BH = PartiallyBayesTest(prior=NPMLE, prior_grid=300, discretize_marginal=true, α = 0.1),
    t_oracle_BH = PartiallyBayesTest(prior=ground_truth_prior, α = 0.1),
    t_oracle_storey = PartiallyBayesTest(prior=ground_truth_prior, α = 0.1,
     multiple_test = BenjaminiHochbergAdaptive(Storey(0.5)))
)

# BH methods

limma_sim = LimmaSimulation(prior =  ground_truth_prior,  ν= 2, adversarial_ordering=true)

method_res = DataFrame(simulation_setting = Symbol[],
                       simulation_number = Int64[],
                       method_name = Symbol[],
                       FDP = Float64[],
                       power = Float64[],
                       rejections = Int64[],
                       α = Float64[]
                       )


simulation_setting = limma_sim

_sim = run_simulation(limma_sim)

key = :t_standard_BH
_method = getproperty(method_list, key)

_apply_method = fit(t_limma_BH, sim.Ss, sim.Zs)




collect( method_list )



res_t_standard_BH = fit(t_standard_BH, sim.Ss, sim.Zs)
res_t_limma_BH = fit(t_limma_BH, sim.Ss, sim.Zs)
res_t_npmle_BH = fit(t_npmle_BH, sim.Ss, sim.Zs)
res_t_oracle_BH = fit(t_oracle_BH, sim.Ss, sim.Zs)
res_t_storey_BH = fit(t_oracle_storey, sim.Ss, sim.Zs)


evaluate_method(sim,  res_t_standard_BH)
evaluate_method(sim,  res_t_limma_BH)
evaluate_method(sim,  res_t_npmle_BH)
evaluate_method(sim,  res_t_oracle_BH)
evaluate_method(sim,  res_t_storey_BH)



# dg = 4, s20 = 4, vg = 1/3 and v0 = 2.
