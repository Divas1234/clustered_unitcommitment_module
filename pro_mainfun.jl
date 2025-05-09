ENV["GKS_ENCODING"] = "utf-8"

using BenchmarkTools, Plots, JLD2, StatsPlots, LaTeXStrings, Plots, PlotThemes
filepath = "D:/GithubClonefiles/clustered_unitcommitment_module/"

include("D:/GithubClonefiles/clustered_unitcommitment_module/src/formatteddata.jl")
include("D:/GithubClonefiles/clustered_unitcommitment_module/src/renewableenergysimulation.jl")
include("D:/GithubClonefiles/clustered_unitcommitment_module/src/showboundrycase.jl")
include("D:/GithubClonefiles/clustered_unitcommitment_module/src/readdatafromexcel.jl")
include("D:/GithubClonefiles/clustered_unitcommitment_module/src/casesploting.jl")
include("D:/GithubClonefiles/clustered_unitcommitment_module/src/creatfrequencyconstraints.jl")
include("D:/GithubClonefiles/clustered_unitcommitment_module/src/cluster_units.jl")
include("D:/GithubClonefiles/clustered_unitcommitment_module/src/recognizingcriticalscenarios.jl")
include("D:/GithubClonefiles/clustered_unitcommitment_module/com/cuccommitmentmodel_bench.jl")
# include("D:/GithubClonefiles/clustered_unitcommitment_module/com/cuccommitmentmodel_pro.jl")

UnitsFreqParam, WindsFreqParam, StrogeData, DataGen, GenCost, DataBranch, LoadCurve, DataLoad, windcurve, pvcurve = readxlssheet()
config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC = forminputdata(DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)

windspmax, pvpmax = 0.250, 2.750
# winds, NW = winds_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 1)
# pvs, NW = pv_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 1)

cunits, cNG, cluster_cunitsset, cluster_featurematrix = calculating_clustered_units(units, DataGen, GenCost, UnitsFreqParam)
winds, NW = winds_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 2)
pvs, NW = pv_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 2)
# pvs.scenarios_curve

# NOTE uc model
rampingup_critical_scenario, frequency_critical_scenario, seqential_rampingrate, netloadcurves = recognizing_critical_scenarios(winds, pvs, loads, NT)
@time ben_p₀, ben_pᵨ, ben_pᵩ, ben_qᵩ, ben_seq_sr⁺, ben_seq_sr⁻, ben_su_cost, ben_sd_cost, ben_prod_cost, ben_cost_sr⁺, ben_cost_sr⁻,
ben_x, ben_z, ben_each_pg₀ = refined_cscucmodel_withoutFreqandFlex(
    NT, NB, NG, cNG, ND, NC, units, cunits, loads, winds, pvs, lines,
    config_param, cluster_cunitsset, cluster_featurematrix, rampingup_critical_scenario, frequency_critical_scenario)
