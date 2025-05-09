using BenchmarkTools

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
windspmax, pvpmax = 0.5, 1.0
winds, NW = genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 1)
# using Plots

# netloadcurve, loadcurve = zeros(size(winds.scenarios_curve, 2), size(loads.load_curve', 1)), zeros(size(winds.scenarios_curve, 2), size(loads.load_curve', 1))
# for t in 1:size(loads.load_curve', 1)
#     for s in 1:size(winds.scenarios_curve, 2)
#         loadcurve[s, t], netloadcurve[s, t] =
#             sum(loads.load_curve[:, t]), sum(loads.load_curve[:, t]) - winds.scenarios_curve[:, t][1, 1]
#     end
# end

# Plots.plot(loadcurve[1, :], ylims=(0, 3.0))
# Plots.plot(netloadcurve[1, :])
# # Plots.plot(reshape(loadcurve[1, :], 24, 7))
# # Plots.plot(reshape(netloadcurve[1, :], 24, 7))
# Plots.plot(winds.scenarios_curve')

rampingup_critical_scenario, frequency_critical_scenario = recognizing_critical_scenarios(winds, loads, NT)
# rampingup_critical_scenario, frequency_critical_scenario = reshape(rampingup_critical_scenario[1, :], Int64(7), Int64(24)), reshape(frequency_critical_scenario[1, :], Int64(7), Int64(24))
# p1 = Plots.heatmap(rampingup_critical_scenario, c=cgrad([:white, :red]), title="critical scenarios for flexility-check", ylabel="scenarios", xlabel="time")
# p2 = Plots.heatmap(frequency_critical_scenario, c=cgrad([:white, :blue]), title="critical scenarios for frequency-dynamics", ylabel="scenarios", xlabel="time")
# Plots.plot(p1, p2; size=(800, 300), titlefontsize=8, layout=(1, 2))

cunits, cNG, cluster_cunitsset, cluster_featurematrix = calculating_clustered_units(units, DataGen, GenCost, UnitsFreqParam)
# sum(cluster_featurematrix[:,2])
p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻ = refined_cscucmodel(
    NT, NB, NG, cNG, ND, NC, units, cunits, loads, winds, lines, config_param, cluster_cunitsset, cluster_featurematrix, rampingup_critical_scenario, frequency_critical_scenario)

# plotcasestudies(p₀,pᵨ,pᵩ,seq_sr⁺,seq_sr⁻,su_cost,sd_cost,prod_cost,cost_sr⁺,cost_sr⁻,NT,NG,ND,NW,NC,)
