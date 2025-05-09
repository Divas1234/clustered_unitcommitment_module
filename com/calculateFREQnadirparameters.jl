# test frequency dynamics

using BenchmarkTools, Plots, JLD2
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/src/formatteddata.jl")
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/src/renewableenergysimulation.jl")
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/src/showboundrycase.jl")
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/src/readdatafromexcel.jl")
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/src/casesploting.jl")
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/src/creatfrequencyconstraints.jl")
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/src/cluster_units.jl")
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/src/recognizingcriticalscenarios.jl")
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/src/creatfrequencyconstraints.jl")
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/com/cuccommitmentmodel_bench.jl")
include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/com/cuccommitmentmodel_pro.jl")

UnitsFreqParam, WindsFreqParam, StrogeData, DataGen, GenCost, DataBranch, LoadCurve, DataLoad, windcurve, pvcurve = readxlssheet()
config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC = forminputdata(DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)

# NT = 24
windspmax, pvpmax = 5.0, 30.0
winds, NW = genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 2)
loadcurve, netloadcurve = zeros(1, size(loads.load_curve', 1)), zeros(size(winds.scenarios_curve, 1), size(loads.load_curve', 1))

# xdata = collect(1:1:168)
# p3 = Plots.plot(xdata, zeros(1, 168)[1, :], size=(300, 200), fillrange=loadcurve[1, :], fillalpha=0.35, c=1, fc=:skyblue, label="load curve", ylims=(0, 100), legend=:topleft)
# p3 = Plots.plot!(xdata, zeros(1, 168)[1, :], fillrange=netloadcurve[1, :], fillalpha=0.35, c=1, fc=:orange, label="netload curve", legend=:topleft)
cunits, cNG, cluster_cunitsset, cluster_featurematrix = calculating_clustered_units(units, DataGen, GenCost, UnitsFreqParam)

rampingup_critical_scenario, frequency_critical_scenario, seqential_rampingrate, netloadcurves = recognizing_critical_scenarios(winds, loads, NT)

@time p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻, x, z = refined_cscucmodel_withFreqandFlex(
    NT, NB, NG, cNG, ND, NC, units, cunits, loads, winds, lines,
    config_param, cluster_cunitsset, cluster_featurematrix,
    rampingup_critical_scenario, frequency_critical_scenario, seqential_rampingrate, netloadcurves)

A, b = creatfrequencyfittingfunction(cunits, winds, cNG, NW, cluster_cunitsset)

# Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Rg, unitsamplestatues = montecalrosimulation(cunits, winds, cNG, NW, cluster_cunitsset)

# Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Kg, set_Rg = zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1)

# for i in 1:sampleNumber
#     Sampling_Statue = x[:, i]
#     M, H, D, T, R, F, K, δp, Hg, Dg = calculate_aggregatedfrequencyparameters(winds, cunits, NW, cNG, cluster_cunitsset, Sampling_Statue)
#     f_nadir = calculate_frequencynadir(M, H, D, T, R, F, K, δp)
#     Set_f_nadir[i, 1], set_H[i, 1], set_δp[i, 1], set_Dg[i, 1], set_Fg[i, 1], set_Kg[i, 1], set_Rg[i, 1] = f_nadir, Hg, δp, Dg, F, K, R
# end

# y = Set_f_nadir
# xx = zeros(size(set_H, 1), 4)
# str = size(set_H, 1)
# for i in 1:str
#     xx[i, 1:4] = [set_H[i, 1], set_Dg[i, 1], set_Fg[i, 1], set_Rg[i, 1]]
# end
# x₀ = hcat((filter(!isnan, xx[:, i]) for i in 1:size(xx, 2))...)
# y₀ = hcat((filter(!isnan, y[:, i]) for i in 1:size(y, 2))...)
# # println(x), println(y)
# sol = llsq(x₀, y₀)

# A, b = sol[1:end-1, :], sol[end, :]

param_H, param_F, param_R, param_δ = A[1], A[2], A[3], b[1, 1]
coff₁ =
    (
        param_H * cunits.Hg +
        param_F * cunits.Kg .* cunits.Fg ./ cunits.Rg +
        param_R * cunits.Kg ./ cunits.Rg
    ) .* cunits.p_max
coff₂ = (param_δ)

tem = zeros(NT, 3)
for t in 1:168
    Sampling_Statue = x[:, t]
    M, H, D, T, R, F, K, δp, Hg, Dg, Rg, Fg = calculate_aggregatedfrequencyparameters(winds, cunits, NW, cNG, cluster_cunitsset, Sampling_Statue)
    f_nadir = calculate_frequencynadir(M, H, D, T, R, F, K, δp)
    tem[t, 1] = f_nadir
    tem[t, 2] = sum(A .* reshape([Hg; Fg; Rg], 3, 1)) + b[1, 1]
    tem[t, 3] = 0.50 - (sum(reshape(coff₁, 3, 1) .* reshape([Hg; F; R], 3, 1) .* Sampling_Statue) / sum(cunits.p_max .* Sampling_Statue) + coff₂) * 0.50
end
Plots.plot(tem[:, 1])
Plots.plot(tem[:, 2])
Plots.plot(tem[:, 3])

t = 1
Sampling_Statue = x[:, t]
M, H, D, T, R, F, K, δp, Hg, Dg = calculate_aggregatedfrequencyparameters(winds, cunits, NW, cNG, cluster_cunitsset, Sampling_Statue)
f_nadir = calculate_frequencynadir(M, H, D, T, R, F, K, δp)
tem[t, 3] = sum(reshape(coff₁, 3, 1) .* Sampling_Statue) / sum(cunits.p_max .* Sampling_Statue) + b[1, 1]
sum(reshape(coff₁, 3, 1) .* Sampling_Statue)
sum(cunits.p_max .* Sampling_Statue)
# @constraint(scuc, [t = 1:NT], sum(coff₁ .* x[:, t]) <= sum(coff₂ * cunits.p_max .* x[:, t]))