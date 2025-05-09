using BenchmarkTools
using Plots
include("src/formatteddata.jl")
include("src/renewableenergysimulation.jl")
include("src/showboundrycase.jl")
include("src/readdatafromexcel.jl")
include("src/cuccommitmentmodel.jl")
include("src/tuccommitmentmodel.jl")
include("src/casesploting.jl")
include("src/creatfrequencyconstraints.jl")
include("src/cluster_units.jl")
include("src/recognizingcriticalscenarios.jl")

UnitsFreqParam, WindsFreqParam, StrogeData, DataGen, GenCost, DataBranch, LoadCurve, DataLoad, windcurve, pvcurve = readxlssheet()
config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC = forminputdata(DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)
windspmax, pvpmax = 10.0, 20.0
winds, NW = genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 1)

rampingup_critical_scenario, frequency_critical_scenario = recognizing_critical_scenarios(winds, loads, NT)
# rampingup_critical_scenario, frequency_critical_scenario = reshape(rampingup_critical_scenario[1, :], Int64(7), Int64(24)), reshape(frequency_critical_scenario[1, :], Int64(7), Int64(24))
p51 = Plots.heatmap(rampingup_critical_scenario, c=cgrad([:white, :red]), title="critical scenarios for flexility-check", ylabel="scenarios", xlabel="time", fa=0.75)
p52 = Plots.heatmap(frequency_critical_scenario, c=cgrad([:white, :blue]), title="critical scenarios for frequency-dynamics", ylabel="scenarios", xlabel="time", fa=0.75)
p5 = Plots.plot(p1, p2; size=(800, 300), titlefontsize=8, layout=(1, 2))

rampingup_critical_scenario

p61 = Plots.heatmap(reshape(rampingup_critical_scenario[1, 61:84], 24, 1)', c=cgrad([:white, :red]), title="critical scenarios for flexility-check", ylabel="scenarios", xlabel="time", fa=0.75)
p62 = Plots.heatmap(reshape(frequency_critical_scenario[1, 61:84], 24, 1)', c=cgrad([:white, :blue]), title="critical scenarios for frequency-dynamics", ylabel="scenarios", xlabel="time", fa=0.75)
p6 = Plots.plot(p1, p2; size=(800, 300), titlefontsize=8, layout=(1, 2))

p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻ = scucmodel(
    NT, NB, NG, ND, NC, units, loads, winds, lines, config_param, rampingup_critical_scenario, frequency_critical_scenario)

cunits, cNG, cluster_cunitsset, cluster_featurematrix = calculating_clustered_units(units, DataGen, GenCost, UnitsFreqParam)
# sum(cluster_featurematrix[:,2])
p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻ = refined_cscucmodel(
    NT, NB, NG, cNG, ND, NC, units, cunits, loads, winds, lines, config_param, cluster_cunitsset, cluster_featurematrix, rampingup_critical_scenario, frequency_critical_scenario)

# plotcasestudies(p₀,pᵨ,pᵩ,seq_sr⁺,seq_sr⁻,su_cost,sd_cost,prod_cost,cost_sr⁺,cost_sr⁻,NT,NG,ND,NW,NC,)
