
summary_tbl =
method_res |>
x ->
    groupby(x, [:simulation_setting, :ν, :adversarial_ordering, :method_name]) |>
    x -> DataFrames.combine(x, nrow, [:FDP, :Power] .=> mean)

using StatsPlots

dirac_nonadv = filter(row -> (row.simulation_setting == :Dirac) & !row.adversarial_ordering, summary_tbl)
inversegamma_nonadv = filter(row -> (row.simulation_setting == :InverseGamma) & !row.adversarial_ordering, summary_tbl)
twopoint_nonadv = filter(row -> (row.simulation_setting == :TwoPointPrior) & !row.adversarial_ordering, summary_tbl)
twopointseparated_nonadv = filter(row -> (row.simulation_setting == :TwoPointPriorSeparated) & !row.adversarial_ordering, summary_tbl)

twopoint_adv = filter(row -> (row.simulation_setting == :TwoPointPrior) & row.adversarial_ordering, summary_tbl)


@df twopointseparated_nonadv plot(:ν, :FDP_mean, group = :method_name, ylim=(0,1), xscale=:log2)
@df twopointseparated_nonadv plot(:ν, :Power_mean, group = :method_name, ylim=(0,1), xscale=:log2)


ground_truth_prior = DiscreteNonParametric([1.0; 0.5], [0.5; 0.5])



ground_truth_prior = Empirikos.InverseScaledChiSquare(4.0, 3.0)

_method =  t_oracle_storey = PartiallyBayesTest(
prior = ground_truth_prior,
α = 0.1,
multiple_test = BenjaminiHochbergAdaptive(Storey(0.5)))

limma_sim = LimmaSimulation(
prior = ground_truth_prior,
ν = 2,
π1 = 0.1,
adversarial_ordering = true,
)

_sim = run_simulation(limma_sim)
ps = Empirikos.limma_pvalue.( _sim.Zs, _sim.Ss, Ref(ground_truth_prior))
histogram(ps)
histogram(ps[_sim.Hs .== 0])

_apply_method = fit(_method, _sim.Ss, _sim.Zs)

evaluate_method(_sim, _apply_method)

estimate(ps, Storey(0.5))

rand(ground_truth_prior, 2)
histogram(ps)
