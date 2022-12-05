using Pkg
Pkg.activate(".")

include("../helpers.jl")
include("simulation_helpers.jl")

using DataFrames
using Random
using JLD2
Random.seed!(1)

νs = [2; 4; 8; 16; 32; 64; 128]

settings = (
    Dirac = Dirac(1.0),
    InverseGamma = Empirikos.InverseScaledChiSquare(1.0, 6.0),
    TwoPointPrior = DiscreteNonParametric([1.0; 0.5], [0.9; 0.1]),
    TwoPointPriorSeparated = DiscreteNonParametric([10.0; 1.0], [0.9; 0.1]),
)


nreps = 2


method_res = DataFrame(
    simulation_setting = Symbol[],
    ν = Int64[],
    adversarial_ordering = Bool[],
    simulation_number = Int64[],
    method_name = Symbol[],
    FDP = Float64[],
    Power = Float64[],
    discoveries = Int64[],
    α = Float64[],
)

for adversarial_ordering in [false; true]
    for simulation_name in keys(settings)
        ground_truth_prior = getproperty(settings, simulation_name)
        for ν in νs

            for k in Base.OneTo(nreps)
                limma_sim = LimmaSimulation(
                    prior = ground_truth_prior,
                    ν = ν,
                    adversarial_ordering = adversarial_ordering,
                )

                method_list = (
                    t_standard_BH = SimultaneousTTest(α = 0.1),
                    t_limma_BH = PartiallyBayesTest(prior = Limma, α = 0.1),
                    t_npmle_BH = PartiallyBayesTest(
                        prior = NPMLE,
                        prior_grid = 300,
                        discretize_marginal = true,
                        α = 0.1,
                    ),
                    t_oracle_BH = PartiallyBayesTest(prior = ground_truth_prior, α = 0.1),
                    t_oracle_storey = PartiallyBayesTest(
                        prior = ground_truth_prior,
                        α = 0.1,
                        multiple_test = BenjaminiHochbergAdaptive(Storey(0.5)),
                    ),
                )

                _sim = run_simulation(limma_sim)
                for key in keys(method_list)
                    _method = getproperty(method_list, key)
                    _apply_method = fit(_method, _sim.Ss, _sim.Zs)
                    _method_eval = evaluate_method(_sim, _apply_method)

                    push!(
                        method_res,
                        (
                            simulation_setting = simulation_name,
                            ν = ν,
                            adversarial_ordering = adversarial_ordering,
                            simulation_number = k,
                            method_name = key,
                            _method_eval...,
                            α = _method.α,
                        ),
                    )
                end
            end
        end
    end
end

jldsave("method_res_x.jld2"; method_res=method_res)
