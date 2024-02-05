using Empirikos

Base.@kwdef struct LimmaSimulation{D,S}
    variance_prior::D
    alternative_effects::S = :conjugate
    v0::Float64 = 16.0
    ν::Int64 = 4
    n::Int64 = 10_000
    π1::Float64 = 0.1
    ordering::Symbol = :random
end


function run_simulation(sim::LimmaSimulation)
    variance_prior = sim.variance_prior
    v0 = sim.v0
    ν = sim.ν
    n = sim.n
    π1 = sim.π1
    ordering = sim.ordering
    alternative_effects = sim.alternative_effects

    n1 = floor(Int, π1*n); n0 = n - n1
    Hs = [ones(Bool, n1); zeros(Bool, n0)]

    σs_squared = rand(variance_prior, n)

    if ordering == :decreasing # old "adversarial" --- means that alternatives have small variance
        σs_squared = sort(σs_squared; rev=false)
    elseif ordering == :increasing #  means that nulls have small variance
        σs_squared = sort(σs_squared; rev=true)
    end

    βs = fill(0.0, n)
    if alternative_effects == :conjugate
        βs[Hs] = rand(Normal(), sum(Hs)) .* sqrt(v0) .* sqrt.(σs_squared[Hs])
    elseif isa(alternative_effects, Distribution)
        βs[Hs] = rand(alternative_effects, sum(Hs)) 
    end
    Ss = ScaledChiSquareSample.(σs_squared ./ ν .* rand(Chisq(ν), n), ν)
    Zs = rand(Normal(), n) .* sqrt.(σs_squared) .+  βs
    samples = NormalChiSquareSample.(Zs, Ss)
    (samples = samples, σs_squared = σs_squared, Hs = Hs)
end



Base.@kwdef struct ExtremeNullSimulation
    ν::Int64 = 4
    n::Int64 = 100
end

function run_simulation(sim::ExtremeNullSimulation)
    ν = sim.ν
    n = sim.n

    Hs = zeros(Bool, n)
    #σs_squared = [1.0; fill(40_000, n-1)]
    σs_squared = rand(DiscreteNonParametric([1.0; 40_000], [0.01; 0.99]), n)
    Ss = ScaledChiSquareSample.(σs_squared ./ ν .* rand(Chisq(ν), n), ν)
    Zs = rand(Normal(), n) .* sqrt.(σs_squared) 
    samples = NormalChiSquareSample.(Zs, Ss)
    (samples = samples, σs_squared = σs_squared, Hs = Hs)
end

function evaluate_method(simulation, method_res)
    discoveries = method_res.total_rejections
    true_discoveries = sum(method_res.rj_idx .&  simulation.Hs)
    false_discoveries = discoveries - true_discoveries
    FDP = false_discoveries / max(discoveries, 1)
    Power = true_discoveries / max( sum(simulation.Hs), 1)
    idx_min_ss = findmin(
        getproperty.(simulation.samples[.! simulation.Hs], :S²))[2]
    cond_cov = method_res.pvalue[.! simulation.Hs][idx_min_ss] ≤ 0.2
    (FDP = FDP, Power=Power, discoveries = discoveries, 
    cond_cov = cond_cov)
end
