using Empirikos

Base.@kwdef struct LimmaSimulation{D}
    prior::D
    v0::Float64 = 9.0
    ν::Int64 = 4
    n::Int64 = 10_000
    π1::Float64 = 0.1
    adversarial_ordering::Bool = false
end


function run_simulation(sim::LimmaSimulation)
    prior = sim.prior
    v0 = sim.v0
    ν = sim.ν
    n = sim.n
    π1 = sim.π1
    adversarial_ordering = sim.adversarial_ordering

    n1 = floor(Int, π1*n); n0 = n - n1
    Hs = [ones(Bool, n1); zeros(Bool, n0)]

    σs_squared = rand(prior, n)

    if adversarial_ordering
        σs_squared = sort(σs_squared; rev=true)
    end

    βs = fill(0.0, n)
    βs[Hs] = rand(Normal(), sum(Hs)) .* sqrt(v0) .* sqrt.(σs_squared[Hs])
    Ss = ScaledChiSquareSample.(σs_squared ./ ν .* rand(Chisq(ν), n), ν)
    Zs = rand(Normal(), n) .* sqrt.(σs_squared) .+  βs
    (Zs = Zs, Ss = Ss, βs = βs, σs_squared = σs_squared, Hs = Hs)
end

function evaluate_method(simulation, method_res)
    discoveries = method_res.total_rejections
    true_discoveries = sum(method_res.rj_idx .&  simulation.Hs)
    false_discoveries = discoveries - true_discoveries
    FDP = false_discoveries / max(discoveries, 1)
    Power = true_discoveries / max( sum(simulation.Hs), 1)
    (FDP = FDP, Power=Power, discoveries = discoveries)
end
