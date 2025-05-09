ENV["GKS_ENCODING"] = "utf-8"

using BenchmarkTools, Plots, JLD2, StatsPlots, LaTeXStrings, Plots, PlotThemes
filepath = "/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/"

include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/formatteddata.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/renewableenergysimulation.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/showboundrycase.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/readdatafromexcel.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/casesploting.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/creatfrequencyconstraints.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/cluster_units.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/recognizingcriticalscenarios.jl")
# include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/com/cuccommitmentmodel_bench.jl")
# include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/com/cuccommitmentmodel_pro.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/tuccommitmentmodel.jl")

UnitsFreqParam, WindsFreqParam, StrogeData, DataGen, GenCost, DataBranch, LoadCurve, DataLoad, windcurve, pvcurve = readxlssheet()
config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC = forminputdata(DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)

windspmax, pvpmax = 0.250, 2.750
# winds, NW = winds_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 1)
# pvs, NW = pv_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 1)

cunits, cNG, cluster_cunitsset, cluster_featurematrix = calculating_clustered_units(units, DataGen, GenCost, UnitsFreqParam)
winds, NW = winds_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 1)
pvs, NW = pv_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 1)
# pvs.scenarios_curve
rampingup_critical_scenario, frequency_critical_scenario, seqential_rampingrate, netloadcurves = recognizing_critical_scenarios(winds, pvs, loads, NT)
@time p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻ =
    traditionalscucmodel(NT, NB, NG, ND, NC, units, loads, winds, pvs, lines, config_param, rampingup_critical_scenario, frequency_critical_scenario)

# cunits, cNG, cluster_cunitsset, cluster_featurematrix = calculating_clustered_units(units, DataGen, GenCost, UnitsFreqParam)

# @btime p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻ = refined_cscucmodel(
#     NT, NB, NG, cNG, ND, NC, units, cunits, loads, winds, lines, config_param, cluster_cunitsset, cluster_featurematrix)

# plotcasestudies(p₀,pᵨ,pᵩ,seq_sr⁺,seq_sr⁻,su_cost,sd_cost,prod_cost,cost_sr⁺,cost_sr⁻,NT,NG,ND,NW,NC,)
