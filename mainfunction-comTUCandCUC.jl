ENV["GKS_ENCODING"] = "utf-8"
using TimerOutputs
using BenchmarkTools, Plots, JLD2, StatsPlots, LaTeXStrings, PlotThemes
filepath = ""
using Random
Random.seed!(1234)

include("src/formatteddata.jl")
include("src/renewableenergysimulation.jl")
include("src/showboundrycase.jl")
include("src/readdatafromexcel.jl")
# include("src/casesploting.jl")
include("src/creatfrequencyconstraints.jl")
include("src/cluster_units.jl")
include("src/recognizingcriticalscenarios.jl")
include("com/cuccommitmentmodel_bench.jl")
# include("com/cuccommitmentmodel_pro.jl")
include("src/tuccommitmentmodel.jl")

UnitsFreqParam, WindsFreqParam, StrogeData, DataGen, GenCost, DataBranch, LoadCurve, DataLoad, windcurve, pvcurve = readxlssheet()
config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC = forminputdata(DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)
windspmax, pvpmax = 0.250, 2.750

cunits, cNG, cluster_cunitsset, cluster_featurematrix = calculating_clustered_units(units, DataGen, GenCost, UnitsFreqParam)
winds, NW = winds_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 2)
pvs, NW = pv_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 2)
rampingup_critical_scenario, frequency_critical_scenario, seqential_rampingrate, netloadcurves = recognizing_critical_scenarios(winds, pvs, loads, NT)

to = TimerOutput()
@timeit to "tuc" traditionalscucmodel(NT, NB, NG, ND, NC, units, loads, winds, pvs, lines, config_param, rampingup_critical_scenario, frequency_critical_scenario)
@timeit to "cuc" refined_cscucmodel_withoutFreqandFlex(
	NT, NB, NG, cNG, ND, NC, units, cunits, loads, winds, pvs, lines,
	config_param, cluster_cunitsset, cluster_featurematrix, rampingup_critical_scenario, frequency_critical_scenario)
show(to; allocations = true)
# @time p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻ =
#     traditionalscucmodel(NT, NB, NG, ND, NC, units, loads, winds, pvs, lines, config_param, rampingup_critical_scenario, frequency_critical_scenario)

# @time ben_p₀, ben_pᵨ, ben_pᵩ, ben_qᵩ, ben_seq_sr⁺, ben_seq_sr⁻, ben_su_cost, ben_sd_cost, ben_prod_cost, ben_cost_sr⁺, ben_cost_sr⁻,
# ben_x, ben_z, ben_each_pg₀ = refined_cscucmodel_withoutFreqandFlex(
#     NT, NB, NG, cNG, ND, NC, units, cunits, loads, winds, pvs, lines,
#     config_param, cluster_cunitsset, cluster_featurematrix, rampingup_critical_scenario, frequency_critical_scenario)
