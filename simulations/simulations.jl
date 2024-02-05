task_id = parse(Int, ARGS[1])  # SLURM_ARRAY_TASK_ID

using Pkg
Pkg.activate(".")


dir = @__DIR__
include(joinpath(dir, "..", "helpers.jl"))
include(joinpath(dir, "simulation_helpers.jl"))


using DataFrames
using Random
using JLD2
using MosekTools
using Setfield

# Monte Carlo replicates
nreps = 1_000


Random.seed!(1)


# settings
νs = [2; 4; 8; 16; 32; 64]

variance_settings = (
    Dirac = Dirac(1.0),
    InverseGamma = Empirikos.InverseScaledChiSquare(1.0, 6.0),
    TwoPointPrior = DiscreteNonParametric([10.0; 1.0], [0.5; 0.5])
)

signal_settings = (conjugate = :conjugate,
                   point_mass = Dirac(4), 
                   normal = Normal(0, 4))

ordering_settings = (
                    random = :random, 
                    decreasing = :decreasing
                    )

variance_keys = keys(variance_settings)
signal_keys = keys(signal_settings)
ordering_keys = keys(ordering_settings)

key_combinations = collect(Iterators.product(variance_keys, signal_keys, ordering_keys))
(variance_key, signal_key, ordering_key) = key_combinations[task_id] #[task_id]


ground_truth_variance_prior = getproperty(variance_settings, variance_key)
signal_setting = getproperty(signal_settings, signal_key)
ordering_setting = getproperty(ordering_settings, ordering_key)

general_limma_sim = LimmaSimulation(
    variance_prior = ground_truth_variance_prior,
    alternative_effects = signal_setting,
    ordering = ordering_setting,
)



basic_npmle = Empirikos.EmpiricalPartiallyBayesTTest(
                prior = DiscretePriorClass(),
                prior_grid_size = 300,
                discretize_marginal = true,
                α = 0.1,
                solver = Mosek.Optimizer,
            )


method_res = DataFrame(
    variance_setting = Symbol[],
    signal_setting = Symbol[],
    ν = Int64[],
    ordering = Symbol[],
    simulation_number = Int64[],
    method_name = Symbol[],
    FDP = Float64[],
    Power = Float64[],
    discoveries = Int64[],
    cond_cov = Bool[],
    α = Float64[],
)

for ν in νs
    for k in Base.OneTo(nreps)
        limma_sim = @set general_limma_sim.ν = ν
        method_list = (
            t_standard_BH = SimultaneousTTest(α = 0.1),
            t_limma_BH = Empirikos.EmpiricalPartiallyBayesTTest(
                prior = Empirikos.Limma(),
                solver = nothing,
                α = 0.1,
            ),
            t_npmle_BH = basic_npmle,
            t_npmle_BH_0 = (@set basic_npmle.lower_quantile = 0.0),
            t_npmle_BH_02 = (@set basic_npmle.lower_quantile = 0.2),
            t_oracle_BH = Empirikos.EmpiricalPartiallyBayesTTest(
                prior = ground_truth_variance_prior,
                solver = nothing,
                α = 0.1,
            ),
            t_oracle_storey = Empirikos.EmpiricalPartiallyBayesTTest(
                prior = ground_truth_variance_prior,
                α = 0.1,
                solver = nothing,
                multiple_test = BenjaminiHochbergAdaptive(Storey(0.5)),
            ),
        )

        _sim = run_simulation(limma_sim)
        for key in keys(method_list)
            _method = getproperty(method_list, key)
            _apply_method = fit(_method, _sim.samples)
            _method_eval = evaluate_method(_sim, _apply_method)

            push!(
                method_res,
                (
                    variance_setting = variance_key,
                    signal_setting = signal_key,
                    ν = ν,
                    ordering = ordering_key,
                    simulation_number = k,
                    method_name = key,
                    _method_eval...,
                    α = _method.α,
                ),
            )
        end
    end
end



jldsave("method_res_$(variance_key)_$(signal_key)_$(ordering_key).jld2"; method_res=method_res)
