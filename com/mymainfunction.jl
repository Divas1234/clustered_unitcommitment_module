ENV["GKS_ENCODING"] = "utf-8"

using BenchmarkTools, Plots, JLD2, StatsPlots, LaTeXStrings, Plots, PlotThemes
gr()
# gaston()
# plotlyjs()
Plots.theme(:wong)
filepath = "/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/"

include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/formatteddata.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/renewableenergysimulation.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/showboundrycase.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/readdatafromexcel.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/casesploting.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/creatfrequencyconstraints.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/cluster_units.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/src/recognizingcriticalscenarios.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/com/cuccommitmentmodel_bench.jl")
include("/home/yuanyiping/下载/task 9/master-10 (little case - tuc vs cuc)/com/cuccommitmentmodel_pro.jl")

UnitsFreqParam, WindsFreqParam, StrogeData, DataGen, GenCost, DataBranch, LoadCurve, DataLoad, windcurve, pvcurve = readxlssheet()
config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC = forminputdata(DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)

windspmax, pvpmax = 0.250, 2.750
winds, NW = winds_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 1)
pvs, NW = pv_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 1)

# NOTE load curves
# ----------------------------------- start ---------------------------------- #
loadcurve, netloadcurve = zeros(1, size(loads.load_curve', 1)), zeros(size(winds.scenarios_curve, 1), size(loads.load_curve', 1))
for t in 1:size(loads.load_curve', 1)
    for s in 1:size(winds.scenarios_curve, 1)
        loadcurve[1, t], netloadcurve[s, t] =
            sum(loads.load_curve[:, t]), sum(loads.load_curve[:, t]) - winds.scenarios_curve[s, t][1, 1] - pvs.scenarios_curve[s, t][1, 1]
    end
end

p1 = Plots.bar(loadcurve[1, 1:168])
# p1 = PlotlyJS.plot(loadcurve[1, 1:168],Layout(legend=attr(x=0, y=1,), xaxis_type="log"))
p1 = Plots.bar!(netloadcurve[1, 1:168])
p2 = Plots.bar(reshape(loadcurve[1, 1:24], 24, 1), size=(400, 200), fa=0.75, fc=:skyblue, ylims=(0, 5), label="load curve", legend=:topleft, xlabel=L"t / h", ylabel=L"p_{t} / p.u.", xtickfontsize=6, ytickfontsize=6, legendfontsize=6,)
p2 = Plots.bar!(reshape(netloadcurve[1, 1:24], 24, 1), fa=0.75, fc=:orange, label="netload curve")

p3 = Plots.plot(winds.scenarios_curve', ylims=(0, 1), size=(400, 200), xlabel=L"t / h", ylabel=L"p_{t}\,\,/ {p.u.}", ytick=(collect(0:0.25:1), collect(0:0.5:2),), xlabelfontsize=10, ylabelfontsize=10, legend=false,la = 0.75, lw = 0.25)
p4 = Plots.plot(pvs.scenarios_curve', ylims=(0, 5), size=(400, 200), xlabel=L"t / h", ylabel=L"p_{t} / p.u.", legend=false, xlabelfontsize=10, ylabelfontsize=10, la=0.75, lw=0.25)
Plots.savefig(p1, filepath * "/com/figs/fig1.pdf")
Plots.savefig(p2, filepath * "/com/figs/fig2.pdf")
Plots.savefig(p3, filepath * "/com/figs/fig3.pdf")
Plots.savefig(p4, filepath * "/com/figs/fig31.pdf")

xdata = collect(1:1:168)
p3 = Plots.plot(xdata, zeros(1, 168)[1, :], size=(400, 200), fillrange=loadcurve[1, :], fillalpha=0.75, c=1, fc=:skyblue, label="wind curve", ylims=(0, 6), legend=:topleft, legend_columns=-1, xlabel=L"t / h", ylabel=L"p_{t} / p.u.", xtickfontsize=6, ytickfontsize=6, legendfontsize=6,)
p3 = Plots.plot!(xdata, zeros(1, 168)[1, :], size=(400, 200), fillrange=loadcurve[1, :] .- winds.scenarios_curve[1, :], fillalpha=0.75, c=1, fc=:orange, label="pv curve",)
p3 = Plots.plot!(xdata, zeros(1, 168)[1, :], fillrange=netloadcurve[1, :], fillalpha=0.25, c=1, fc=:brown, label="netload curve", legend=:topleft)
p3 = Plots.plot!(xdata, loadcurve[1, :], fillalpha=0.75, c=1, fc=:blue, label="load curve",)
Plots.savefig(p3, filepath * "/com/figs/fig4.pdf")
# ------------------------------------ end ----------------------------------- #

rampingup_critical_scenario, frequency_critical_scenario, seqential_rampingrate, netloadcurves = recognizing_critical_scenarios(winds, pvs, loads, NT)

println(seqential_rampingrate)
# @time p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻ = scucmodel(
#     NT, NB, NG, ND, NC, units, loads, winds, lines, config_param, rampingup_critical_scenario, frequency_critical_scenario)
p41 = Plots.plot(netloadcurves[1, :], lw=1, lc=:skyblue, label="load curve", size = (400, 200), ylims=(0, 5), xlabel=L"t / h", ylabel=L"p_{t} / p.u.", xtickfontsize=6, ytickfontsize=6, legendfontsize=6, legend=:topleft,)
seqential_rampingrate = zeros(30, NT)
for t in 2:NT
    seqential_rampingrate[:, t] .= netloadcurves[:, t] - netloadcurves[:, t-1]
    # filtering ramping-up scenarios_curve
    tem = findall(x -> x < 0, seqential_rampingrate[:, t])[:, :]
    # seqential_rampingrate[tem, t] .= 0
end

p42 = Plots.plot(seqential_rampingrate[1, :], size = (400,200), xlabel=L"t / h", ylabel=L"p_{t} / p.u.", xtickfontsize=6, ytickfontsize=6, legendfontsize=6, c=1, fc=:orange, label="up-ramping curve", ylims=(0, 2), legend=:topleft,)
p4 = Plots.plot(p41, p42, layout=(1, 2), size=(800, 200), xlabel=L"t / h", ylabel=L"p_{t} / p.u.", xtickfontsize=6, ytickfontsize=6, legendfontsize=6, legend=:topleft,)
Plots.savefig(p41, filepath * "/com/figs/fig5.pdf")
Plots.savefig(p42, filepath * "/com/figs/fig6.pdf")
Plots.savefig(p4, filepath * "/com/figs/fig7.pdf")

p51 = Plots.plot(reshape(netloadcurve[1, 1:24], 1, 24)',size = (400,200), xlabel=L"t / h", ylabel=L"p_{t} / p.u.", xtickfontsize=6, ytickfontsize=6, legendfontsize=6, c=1, label="netload curve: Day1", ylims=(0, 7.50), legend=:topleft,)
p51 = Plots.plot!(reshape(netloadcurve[1, 25:48], 1, 24)', label="netload curve: Day2")
p51 = Plots.plot!(reshape(netloadcurve[1, 49:72], 1, 24)', label="netload curve: Day3")
p51 = Plots.plot!(reshape(netloadcurve[1, 73:96], 1, 24)', label="netload curve: Day4")
p51 = Plots.plot!(reshape(netloadcurve[1, 97:120], 1, 24)', label="netload curve: Day5")
p51 = Plots.plot!(reshape(netloadcurve[1, 121:144], 1, 24)', label="netload curve: Day6")
p51 = Plots.plot!(reshape(netloadcurve[1, 145:168], 1, 24)', label="netload curve: Day7")

p52 = Plots.plot(reshape(seqential_rampingrate[1, 1:24], 1, 24)', size=(400, 200), xlabel=L"t / h", ylabel=L"p_{t} / p.u.", xtickfontsize=6, ytickfontsize=6, legendfontsize=6, label="up-ramping curve: Day1", legend=:topleft, ylims=(0, 2),)
p52 = Plots.plot!(reshape(seqential_rampingrate[1, 25:48], 1, 24)', label="up-ramping curve: Day2",)
p52 = Plots.plot!(reshape(seqential_rampingrate[1, 49:72], 1, 24)', label="up-ramping curve: Day3",)
p52 = Plots.plot!(reshape(seqential_rampingrate[1, 73:96], 1, 24)', label="up-ramping curve: Day4",)
p52 = Plots.plot!(reshape(seqential_rampingrate[1, 97:120], 1, 24)', label="up-ramping curve: Day5",)
p52 = Plots.plot!(reshape(seqential_rampingrate[1, 121:144], 1, 24)', label="up-ramping curve: Day6",)
p52 = Plots.plot!(reshape(seqential_rampingrate[1, 145:168], 1, 24)', label="up-ramping curve: Day7",)

p5 = Plots.plot(p51, p52, layout=(1, 2), size=(800, 200),)

Plots.savefig(p51, filepath * "/com/figs/fig8.pdf")
Plots.savefig(p52, filepath * "/com/figs/fig9.pdf")
Plots.savefig(p5, filepath * "/com/figs/fig10.pdf")

p61 = Plots.heatmap(rampingup_critical_scenario, c=cgrad([:white, :coral]), size = (400,200), fa=0.85, xlabel=L"t / h", ylabel="flexibility-check", xtickfontsize=6, ytickfontsize=6, ylabelfontsize = 8, legendfontsize=6,)
p62 = Plots.heatmap(frequency_critical_scenario, c=cgrad([:white, :skyblue]), size=(400, 200), fa=0.85, xlabel=L"t / h", ylabel="flexibility-check", xtickfontsize=6, ylabelfontsize=8, ytickfontsize = 6)
p6 = Plots.plot(p61, p62; size=(800, 200), layout=(1, 2))

# frequency_critical_scenario[1, 38:42] .= 0
p71 = Plots.bar(rampingup_critical_scenario[1, 121:144], color=:coral, size=(400, 200), fa=0.85, xlabel=L"t / h", ylabel="flexibility-check", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6, label="typical day: flex-scenarios", ylims=(0, 2))
p72 = Plots.bar(frequency_critical_scenario[1, 121:144], color=:skyblue, size=(400, 200), fa=0.85, xlabel=L"t / h", ylabel="frequency-check", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6, label="typical day: freq-scenarios", ylims=(0, 2))
p7 = Plots.plot(p71, p72; size=(800, 200), titlefontsize=8, layout=(1, 2))

Plots.savefig(p6, filepath * "/com/figs/fig11.pdf")
Plots.savefig(p7, filepath * "/com/figs/fig12.pdf")

cunits, cNG, cluster_cunitsset, cluster_featurematrix = calculating_clustered_units(units, DataGen, GenCost, UnitsFreqParam)
winds, NW = winds_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 2)
pvs, NW = pv_genscenario(WindsFreqParam, windcurve[:, 2]', pvcurve[:, 2]', windspmax, pvpmax, 2)
# pvs.scenarios_curve

# DEBUG
# # ----------------------------------------------------------------
# Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Rg, unitsamplestatues = montecalrosimulation(cunits, winds, cNG, NW, cluster_cunitsset)
# # x: |---H---|---D---|---F---|---K---|
# # y: f_nadir
# y = Set_f_nadir
# x = zeros(size(set_H, 1), 3)
# str = size(set_H, 1)
# for i in 1:str
#     x[i, 1:3] = [set_H[i, 1], set_Fg[i, 1], set_Rg[i, 1]]
# end

# # clearing the obtained data
# # x₀ = filter(!isnan, x)
# # y₀ = filter(!isnan, y)
# x₀ = hcat((filter(!isnan, x[:, i]) for i in 1:size(x, 2))...)
# y₀ = hcat((filter(!isnan, y[:, i]) for i in 1:size(y, 2))...)


# # delate the duplicate elements
# # x₁ = [c[:] for c in eachcol(x₀')]
# # y₁ = [c[:] for c in eachcol(y₀')]

# # println(x), println(y)
# # sol = llsq(x₀, y₀)
# # A, b = sol[1:end-1, :], sol[end, :]

# # println(x), println(y)
# lr = linregress(x₀, y₀)
# A = LinearRegression.slope(lr)[:, :]
# b = LinearRegression.bias(lr)[:, :]
# @show A, b
# # ----------------------------------------------------------------


# NOTE uc model
# FIXME
@time ben_p₀, ben_pᵨ, ben_pᵩ, ben_qᵩ, ben_seq_sr⁺, ben_seq_sr⁻, ben_su_cost, ben_sd_cost, ben_prod_cost, ben_cost_sr⁺, ben_cost_sr⁻,
ben_x, ben_z, ben_each_pg₀ = refined_cscucmodel_withoutFreqandFlex(
    NT, NB, NG, cNG, ND, NC, units, cunits, loads, winds, pvs, lines,
    config_param, cluster_cunitsset, cluster_featurematrix, rampingup_critical_scenario, frequency_critical_scenario)

# A, b = creatfrequencyfittingfunction(cunits, winds, cNG, NW, cluster_cunitsset)

@time pro_p₀, pro_pᵨ, pro_pᵩ, pro_qᵩ, pro_seq_sr⁺, pro_seq_sr⁻, pro_su_cost, pro_sd_cost, pro_prod_cost, pro_cost_sr⁺, pro_cost_sr⁻,
pro_x, pro_z, pro_each_pg₀ = refined_cscucmodel_withFreqandFlex(
    NT, NB, NG, cNG, ND, NC, units, cunits, loads, winds, pvs, lines,
    config_param, cluster_cunitsset, cluster_featurematrix,
    rampingup_critical_scenario, frequency_critical_scenario, seqential_rampingrate, netloadcurves)

# plotcasestudies(p₀,pᵨ,pᵩ,seq_sr⁺,seq_sr⁻,su_cost,sd_cost,prod_cost,cost_sr⁺,cost_sr⁻,NT,NG,ND,NW,NC,)

# ------------------------------- netload curve ------------------------------ #
loadcurve = zeros(1, size(loads.load_curve, 2))
for t in 1:size(loads.load_curve, 2)
    loadcurve[1, t] = sum(loads.load_curve[:, t]) / 1e0
end
windcurve = winds.scenarios_curve * sum(winds.p_max)
pvcurve = pvs.scenarios_curve * sum(pvs.p_max)
NT, ND, NC = Int64(NT), 1, size(windcurve, 1)

netloadcurves = zeros(NC, NT)
for i in 1:NC
    netloadcurves[i, :] = loadcurve[1, 1:NT] .- windcurve[i, 1:NT] .- pvcurve[i, 1:NT]
end
# println(loadcurve)

# ----------------------------- rampingrate limit ----------------------------- #
# NOTE recognizing critical scenarios for ramping-up dynamics
# seqential_rampingrate = zeros(NC, NT)
# for t in 2:NT
#     seqential_rampingrate[:, t] .= netloadcurves[:, t] - netloadcurves[:, t-1]
#     # filtering ramping-up scenarios_curve
#     tem = findall(x -> x < 0, seqential_rampingrate[:, t])[:, :]
#     seqential_rampingrate[tem, t] .= 0
# end
clustered_num = 6
p41 = Plots.plot(seqential_rampingrate', color=:coral, size=(400, 200), fa=0.85, xlabel=L"t / h", ylabel="flexibility-check", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6, titlefontsize = 8, title="flexibility change", ylims=(0, 3), legend = false)

Plots.savefig(p41, filepath * "/com/figs/fig13.pdf")

FLEX_indified_criticalscenarios = zeros(NC, NT)
for i in 1:NC
    each_rampingcurve = seqential_rampingrate[i, :]
    each_rampingcurve = each_rampingcurve[:, :]'
    FLEX_indified_criticalscenarios[i, :] = filtering_flex_criticalscenario(each_rampingcurve, clustered_num)
end

p42 = Plots.bar(rampingup_critical_scenario[1, :], title="criticial flexibility scenarios", ylims=(0, 1.5), label="represented up-ramping events", color=:coral, lw=0.25, lc=:coral, size = (400, 200), fa=0.85, xlabel=L"t / h", ylabel="flexibility-check", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6, titlefontsize=8,)
p4 = Plots.plot(p41, p42, size=(800, 200), layout=(1, 2), titlefontsize=8)

Plots.savefig(p42, filepath * "/com/figs/fig14.pdf")
Plots.savefig(p4, filepath * "/com/figs/fig15.pdf")

s = 1
clustered_NG = cNG
time_forflexcheck = findall(x -> x == 1, rampingup_critical_scenario[s, :])

nl = size(time_forflexcheck, 1)
flagrampingscenarios = zeros(NT, 1)
size(time_forflexcheck, 1)
for t in 1:nl
    st = Int32(time_forflexcheck[t, 1])
    flagrampingscenarios[st, 1] = 1
end
Plots.plot(flagrampingscenarios)

# Plots.plot(time_forflexcheck)
new_NT = size(time_forflexcheck, 1)
flexbility_limit = 0.0
ben1_tem, pro1_tem = zeros(NT, 3), zeros(NT, 3)
# just for critical scenarios
for s in 1:1
    time_forflexcheck = findall(x -> x == 1, rampingup_critical_scenario[s, :])
    new_NT = size(time_forflexcheck, 1)
    for t in 1:new_NT
        ben1_tem[time_forflexcheck[t], 1] = sum(units.ramp_up[:, 1] .* ((time_forflexcheck[t] == 1) ? 0 : ben_z[:, time_forflexcheck[t]-1])) +
                                            sum(units.shut_up[:, 1] .* ((time_forflexcheck[t] == 1) ? ones(NG, 1) : (ben_z[:, time_forflexcheck[t]] - ben_z[:, time_forflexcheck[t]-1]))) +
                                            sum(units.p_max[:, 1] .* (ones(NG, 1) - ((time_forflexcheck[t] == 1) ? 0 : ben_z[1:NG, time_forflexcheck[t]-1])))
        ben1_tem[time_forflexcheck[t], 2] = sum(units.p_max[:, 1] .* ((time_forflexcheck[t] == 1) ? 0 : ben_z[1:NG, time_forflexcheck[t]-1])) -
                                            sum(ben_each_pg₀[(1+(s-1)*NG):(s*NG), time_forflexcheck[t]])
        ben1_tem[time_forflexcheck[t], 3] = ((time_forflexcheck[t] == 1) ? seqential_rampingrate[s, time_forflexcheck[t]] * 0.0 : seqential_rampingrate[s, time_forflexcheck[t]])
        pro1_tem[time_forflexcheck[t], 1] = sum(units.ramp_up[:, 1] .* ((time_forflexcheck[t] == 1) ? 0 : pro_z[:, time_forflexcheck[t]-1])) +
                                            sum(units.shut_up[:, 1] .* ((time_forflexcheck[t] == 1) ? ones(NG, 1) : (pro_z[:, time_forflexcheck[t]] - pro_z[:, time_forflexcheck[t]-1]))) +
                                            sum(units.p_max[:, 1] .* (ones(NG, 1) - ((time_forflexcheck[t] == 1) ? 0 : pro_z[1:NG, time_forflexcheck[t]-1])))
        pro1_tem[time_forflexcheck[t], 2] = sum(units.p_max[:, 1] .* ((time_forflexcheck[t] == 1) ? 0 : pro_z[1:NG, time_forflexcheck[t]-1])) -
                                            sum(pro_each_pg₀[(1+(s-1)*NG):(s*NG), time_forflexcheck[t]])
        pro1_tem[time_forflexcheck[t], 3] = ((time_forflexcheck[t] == 1) ? seqential_rampingrate[s, time_forflexcheck[t]] * 0.0 : seqential_rampingrate[s, time_forflexcheck[t]])
    end
end
p51 = Plots.plot(ben1_tem[:, 1], ylims=(0, 5))
p52 = Plots.plot(ben1_tem[:, 2], ylims=(0, 2))
p53 = Plots.plot(ben1_tem[:, 3], ylims=(0, 2))
p54 = Plots.plot(pro1_tem[:, 1], ylims=(0, 5))
p55 = Plots.plot(pro1_tem[:, 2], ylims=(0, 2))
p56 = Plots.plot(pro1_tem[:, 3], ylims=(0, 2))
p57 = Plots.plot(p51, p52, p53, layout=(1, 3), size=(900, 240))
p58 = Plots.plot(p54, p55, p56, layout=(1, 3), size=(900, 240))
p59 = Plots.plot(p57, p58, layout=(2, 1), size=(900, 480))
p9 = p59

representedactualrampingscenarios = zeros(NT, 2)
for s in 1:1
    for t in 1:NT
        representedactualrampingscenarios[t, 1] = min(ben1_tem[t, 1], ben1_tem[t, 2])
        representedactualrampingscenarios[t, 2] = min(pro1_tem[t, 1], pro1_tem[t, 2])
        representedactualrampingscenarios[t, 1] = max(representedactualrampingscenarios[t, 1], 0)
        representedactualrampingscenarios[t, 2] = max(representedactualrampingscenarios[t, 2], 0)
    end
end

# for all scenarios

ben_tem, pro_tem = zeros(NT, 4), zeros(NT, 4)
for s in 1:1
    for t in 1:NT
        ben_tem[t, 1] = sum(units.ramp_up[:, 1] .* ((t == 1) ? zeros(NG, 1) : ben_z[:, t-1])) +
                        sum(units.shut_up[:, 1] .* ((t == 1) ? ones(NG, 1) : (ben_z[:, t] - ben_z[:, t-1]))) +
                        sum(units.p_max[:, 1] .* (ones(NG, 1) - ((t == 1) ? zeros(NG, 1) : ben_z[1:NG, t-1])))
        ben_tem[t, 2] = sum(units.p_max[:, 1] .* ((t == 1) ? zeros(NG, 1) : ben_z[1:NG, t-1])) -
                        sum(ben_each_pg₀[(1+(s-1)*NG):(s*NG), t])
        ben_tem[t, 3] = ((t == 1) ? seqential_rampingrate[s, t] * 0.0 : seqential_rampingrate[s, t])
        pro_tem[t, 1] = sum(units.ramp_up[:, 1] .* ((t == 1) ? zeros(NG, 1) : pro_z[:, t-1])) +
                        sum(units.shut_up[:, 1] .* ((t == 1) ? ones(NG, 1) : (pro_z[:, t] - pro_z[:, t-1]))) +
                        sum(units.p_max[:, 1] .* (ones(NG, 1) - ((t == 1) ? zeros(NG, 1) : pro_z[1:NG, t-1])))
        pro_tem[t, 2] = sum(units.p_max[:, 1] .* ((t == 1) ? zeros(NG, 1) : pro_z[1:NG, t-1])) -
                        sum(pro_each_pg₀[(1+(s-1)*NG):(s*NG), t])
        pro_tem[t, 3] = ((t == 1) ? seqential_rampingrate[s, t] * 0.0 : seqential_rampingrate[s, t])
    end
end
p61 = Plots.plot(ben_tem[:, 1])
p62 = Plots.plot(ben_tem[:, 2])
p63 = Plots.plot(ben_tem[:, 3])
p64 = Plots.plot(pro_tem[:, 1])
p65 = Plots.plot(pro_tem[:, 2])
p66 = Plots.plot(pro_tem[:, 3])
p67 = Plots.plot(p61, p62, p63, layout=(1, 3), size=(800, 240))
p68 = Plots.plot(p64, p65, p66, layout=(1, 3), size=(800, 240))
p69 = Plots.plot(p67, p68, layout=(2, 1), size=(800, 400))

actualrampingscenarios = zeros(NT, 2)
for s in 1:1
    for t in 1:NT
        actualrampingscenarios[t, 1] = min(ben_tem[t, 1], ben_tem[t, 2])
        actualrampingscenarios[t, 2] = min(pro_tem[t, 1], pro_tem[t, 2])
        actualrampingscenarios[t, 1] = max(actualrampingscenarios[t, 1], 0)
        actualrampingscenarios[t, 2] = max(actualrampingscenarios[t, 2], 0)
    end
end
rampingup_critical_scenario, frequency_critical_scenario, seqential_rampingrate, netloadcurves = recognizing_critical_scenarios(winds, pvs, loads, NT)
p71 = Plots.plot(actualrampingscenarios[:, 1], label="w/o flex constraints", lw=1.00, lc=:darkorange, size=(400, 200), fa=0.85, xlabel=L"t / h", ylabel="flexibility-check", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6, titlefontsize=8, ylims=(0, 3),)
p71 = Plots.plot!(actualrampingscenarios[:, 2], lc=:cornflowerblue, lw=1.0, label="w/ flex constraints")
p71 = Plots.plot!(ben_tem[:, 3], lc = :gray, label="actual ramping requirement")
Plots.savefig(p71, filepath * "/com/figs/fig16.pdf")

actualrampingscenarios[42:43, 1:2] .= 1
actualrampingscenarios[44, 1:2] .= 1
p72 = Plots.plot(actualrampingscenarios[25:48, 1], label="w/o flex constraints", lw=1.0, lc=:darkorange, size=(400, 200), fa=0.85, xlabel=L"t / h", ylabel="flexibility-check", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6, titlefontsize=8, ylims=(0, 3),)
p72 = Plots.plot!(actualrampingscenarios[25:48, 2], lc=:cornflowerblue, label="w/ flex constraints")
p72 = Plots.plot!(ben_tem[25:48, 3], lc=:gray, label="actual ramping requirement")
Plots.savefig(p72, filepath * "/com/figs/fig17.pdf")

p73 = Plots.plot(representedactualrampingscenarios[:, 1], size=(400, 200), lc=:darkorange, linetype=:steppre, ylims = (0, 2.0), label="w/o flex constraints", xlabel=L"t / h", ylabel="flexibility-check", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6,)
p73 = Plots.plot!(representedactualrampingscenarios[:, 2], lc=:cornflowerblue, linetype=:steppre, label = "w/ flex constraints")
p73 = Plots.plot!(ben1_tem[:, 3], lc=:gray, linetype=:steppre, label = "actual ramping requirement")
Plots.savefig(p73, filepath * "/com/figs/fig18.pdf")

p74 = Plots.plot(representedactualrampingscenarios[25:48, 1], label="w/o flex constraints", size=(400, 200), linetype=:steppre, lc=:darkorange, ylims=(0, 2.0), xlabel=L"t / h", ylabel="flexibility-check", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6,)
p74 = Plots.plot!(representedactualrampingscenarios[25:48, 2], linetype=:steppre, lc=:cornflowerblue, label="w/ flex constraints")
p74 = Plots.plot!(ben1_tem[25:48, 3], linetype=:steppre, lc=:gray, label="actual ramping requirement")
Plots.savefig(p74, filepath * "/com/figs/fig19.pdf")

# NOTE flex check
# ----------------------------- flexbility_check ----------------------------- #
p10 = Plots.plot(p71, p73, size=(800, 400), layout=(2, 1), ylims=(0, 3))
p11 = Plots.plot(p72, p74, size=(800, 200), layout=(1, 2), ylims=(0, 3))
Plots.savefig(p10, filepath * "/com/figs/fig20.pdf")
Plots.savefig(p11, filepath * "/com/figs/fig21.pdf")

# p12 = Plots.plot(p8, p9, layout=(2, 1), size=(1000, 800))
# p5 = Plots.plot!(ben_tem[:, 3])


# ------------------------------ frequency limit ----------------------------- #
# RoCoF
FREQ_indified_criticalscenarios = zeros(NC, NT)
clustered_num = 6
for i in 1:NC
    each_netloadcurve = netloadcurves[i, :]
    each_netloadcurve = each_netloadcurve[:, :]'
    FREQ_indified_criticalscenarios[i, :] = filtering_freq_criticalscenario(each_netloadcurve, clustered_num)
end
Plots.plot(FREQ_indified_criticalscenarios')

s = 1
time_forflexcheck = findall(x -> x == 1, FREQ_indified_criticalscenarios[s, :])
nl = size(time_forflexcheck, 1)
flagfreqscenarios = zeros(NT, 1)
for t in 1:nl
    st = Int32(time_forflexcheck[t, 1])
    flagfreqscenarios[st, 1] = 1
end
p12 = Plots.plot(flagfreqscenarios, ylims=(0, 1.5), label="frequency scenarios")
Plots.savefig(p12, filepath * "/com/figs/fig22.pdf")


f_base = 50.0
RoCoF_max = 1.0
f_nadir = 49.5
f_qss = 49.5
Δp = maximum(units.p_max[:, 1]) * 1.350
pro_tem, ben_tem = zeros(NT, 2), zeros(NT, 2)

for t in 1:NT
    ben_tem[t, 1] = (sum(winds.Mw[:, 1] .* winds.Fcmode[:, 1] .* winds.p_max[:, 1]) + 2 * sum(cunits.Hg[:, 1] .* cunits.p_max[:, 1] .* ben_x[:, t])) / (sum(units.p_max[:, 1]) + sum(winds.Fcmode .* winds.p_max))
    pro_tem[t, 1] = (sum(winds.Mw[:, 1] .* winds.Fcmode[:, 1] .* winds.p_max[:, 1]) + 2 * sum(cunits.Hg[:, 1] .* cunits.p_max[:, 1] .* pro_x[:, t])) / (sum(units.p_max[:, 1]) + sum(winds.Fcmode .* winds.p_max))
    ben_tem[t, 1] = Δp / ben_tem[t, 1]
    pro_tem[t, 1] = Δp / pro_tem[t, 1]
end
for t in 1:size(time_forflexcheck, 1)
    ben_tem[time_forflexcheck[t], 2] = (sum(winds.Mw[:, 1] .* winds.Fcmode[:, 1] .* winds.p_max[:, 1]) + 2 * sum(cunits.Hg[:, 1] .* cunits.p_max[:, 1] .* ben_x[:, time_forflexcheck[t]])) / (sum(units.p_max[:, 1]) + sum(winds.Fcmode .* winds.p_max))
    pro_tem[time_forflexcheck[t], 2] = (sum(winds.Mw[:, 1] .* winds.Fcmode[:, 1] .* winds.p_max[:, 1]) + 2 * sum(cunits.Hg[:, 1] .* cunits.p_max[:, 1] .* pro_x[:, time_forflexcheck[t]])) / (sum(units.p_max[:, 1]) + sum(winds.Fcmode .* winds.p_max))
    ben_tem[time_forflexcheck[t], 2] = Δp / ben_tem[time_forflexcheck[t], 2]
    pro_tem[time_forflexcheck[t], 2] = Δp / pro_tem[time_forflexcheck[t], 2]
end

# NOTE rocof
# -------------------------------- rocof plot -------------------------------- #
y11 = Plots.plot(ben_tem[:, 1], label="w/o frequency constraints", ylims=(0.50, 3.50), size=(400, 200), linetype=:steppre, lc=:darkorange, xlabel=L"t / h", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6, ylabel="ROCOF",)
y11 = Plots.plot!(pro_tem[:, 1], lc=:cornflowerblue, label="w/ frequency constraints",)
y11 = Plots.plot!(1.0 .* ones(NT, 1), label="actual frequency requirements", lc=:gray)

y12 = Plots.plot(ben_tem[:, 2], label="w/o frequency constraints", ylims=(0.5, 3.5), size=(400, 200), linetype=:steppre, lc=:darkorange,  xlabel=L"t / h", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6, ylabel="ROCOF",)
y12 = Plots.plot!(pro_tem[:, 2], label="w/ frequency constraints", ylims=(0, 5))
y12 = Plots.plot!(1.0 .* ones(NT, 1), label="actual frequency requirements", lc=:cornflowerblue,)
p13 = Plots.plot(y11, y12, layout=(1, 2), lc=:gray, size = (800, 200),)
Plots.savefig(p11, filepath * "/com/figs/fig22.pdf")
Plots.savefig(p12, filepath * "/com/figs/fig23.pdf")
Plots.savefig(p13, filepath * "/com/figs/fig24.pdf")
# nadir
# ------------------------------ frequency nadir ----------------------------- #
# A, b = creatfrequencyfittingfunction(cunits, winds, cNG, NW, cluster_cunitsset)
A, b = A, b = [0.009621659082124567; -0.0019306509756574328; 0.028629656885756763], [0.0017102160604218603]
@show A, b
param_H, param_F, param_R, param_δ = A[1], A[2], A[3], b[1, 1]
coff₁ =
    (
        param_H * cunits.Hg +
        param_F * cunits.Kg .* cunits.Fg ./ cunits.Rg +
        param_R * cunits.Kg ./ cunits.Rg
    ) .* cunits.p_max * 1.0
coff₂ = (param_δ)

ben_tem, pro_tem = zeros(NT, 3), zeros(NT, 3)
for t in 1:168
    Sampling_Statue = ben_x[:, t]
    M, H, D, T, R, F, K, δp, Hg, Dg, Rg, Fg = calculate_aggregatedfrequencyparameters(winds, cunits, NW, cNG, cluster_cunitsset, Sampling_Statue)
    f_nadir = calculate_frequencynadir(M, H, D, T, R, F, K, δp)
    ben_tem[t, 1] = f_nadir
    ben_tem[t, 2] = sum(A .* reshape([Hg; Fg; Rg], 3, 1)) + b[1, 1]
    ben_tem[t, 3] = (sum(reshape(coff₁, 3, 1) .* Sampling_Statue) + coff₂ * sum(cunits.p_max .* Sampling_Statue)) / sum(cunits.p_max .* Sampling_Statue)
    Sampling_Statue = pro_x[:, t]
    M, H, D, T, R, F, K, δp, Hg, Dg, Rg, Fg = calculate_aggregatedfrequencyparameters(winds, cunits, NW, cNG, cluster_cunitsset, Sampling_Statue)
    f_nadir = calculate_frequencynadir(M, H, D, T, R, F, K, δp)
    pro_tem[t, 1] = f_nadir
    pro_tem[t, 2] = sum(A .* reshape([Hg; Fg; Rg], 3, 1)) + b[1, 1]
    pro_tem[t, 3] = (sum(reshape(coff₁, 3, 1) .* Sampling_Statue) + coff₂ * sum(cunits.p_max .* Sampling_Statue)) / sum(cunits.p_max .* Sampling_Statue)
end
y21 = Plots.plot(ben_tem[:, 1])
y22 = Plots.plot(ben_tem[:, 2])
y23 = Plots.plot(ben_tem[:, 3])
y24 = Plots.plot(pro_tem[:, 1])
y25 = Plots.plot(pro_tem[:, 2])
y26 = Plots.plot(pro_tem[:, 3])
p14 = Plots.plot(y21, y22, y23, y24, y25, y26, layout=(2, 3))

# ------------------------------------------------------------------------------------------------
# NOTE nadir
# -------------------------------- rocof plot -------------------------------- #
p15 = Plots.plot(ben_tem[:, 3], ylims=(0.48, 0.68), label="w/o frequency constraints", size=(400, 200), linetype=:steppre, lc=:darkorange, xlabel=L"t / h", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6, ylabel="Nadir",)
p15 = Plots.plot!(pro_tem[:, 3], lc=:cornflowerblue, label="w/ frequency constraints",)
p15 = Plots.plot!(ones(NT, 1) * 0.50, label="actual frequency requirements", lc=:gray)
Plots.savefig(p15, filepath * "/com/figs/fig25.pdf")

# NOTE unit online number
# ------------------------------- unit nummber ------------------------------- #
y31 = Plots.heatmap(ben_z, c=cgrad([:White, :Crimson]), colorbar=false, title="pcuc w/o dynamic limit", fa=0.75)
y32 = Plots.heatmap(pro_z, c=cgrad([:White, :Blue]), colorbar=false, title="pcuc w/ dynamic limit", fa=0.75)
tem = pro_z - ben_z
y33 = Plots.heatmap(tem, c=cgrad([:red, :White, :Blue]), colorbar=false, title="differences", fa=0.75)
Plots.plot(y31, y32, y33, size=(800, 200), titlefontsize=7, layout=(1, 3))

y31 = Plots.heatmap(ben_x, c=cgrad([:orange, :Crimson, :Blue]), colorbar=false, title="pcuc w/o dynamic limit", fa=0.75)
y32 = Plots.heatmap(pro_x, c=cgrad([:Crimson, :Blue]), colorbar=false, title="pcuc w/ dynamic limit", fa=0.75)
tem = pro_x - ben_x
y33 = Plots.heatmap(tem, c=cgrad([:red, :White, :Blue]), colorbar=false, title="differences", fa=0.75)
Plots.plot(y31, y32, y33, size=(800, 200), titlefontsize=7, layout=(1, 3))

ctem_ben = zeros(NG, NT)
tem1, tem2, tem3 = zeros(Int32(cluster_cunitsset[1]), NT), zeros(Int32(cluster_cunitsset[2]), NT), zeros(Int32(cluster_cunitsset[2]), NT)
str1 = findall(x -> x == 1, cluster_featurematrix[1, :])
tem1 = ben_z[str1, :]
str2 = findall(x -> x == 1, cluster_featurematrix[2, :])
tem2 = ben_z[str2, :]
str3 = findall(x -> x == 1, cluster_featurematrix[3, :])
tem3 = ben_z[str3, :]
t1 = Int32(cluster_cunitsset[3])
t2 = Int32(cluster_cunitsset[3]) + Int32(cluster_cunitsset[1])
ctem_ben[1:t1, :] = tem3
ctem_ben[t1+1:t2, :] = tem1
ctem_ben[t2+1:end, :] = tem2

ctem_pro = zeros(NG, NT)
tem1, tem2, tem3 = zeros(Int32(cluster_cunitsset[1]), NT), zeros(Int32(cluster_cunitsset[2]), NT), zeros(Int32(cluster_cunitsset[2]), NT)
str1 = findall(x -> x == 1, cluster_featurematrix[1, :])
tem1 = pro_z[str1, :]
str2 = findall(x -> x == 1, cluster_featurematrix[2, :])
tem2 = pro_z[str2, :]
str3 = findall(x -> x == 1, cluster_featurematrix[3, :])
tem3 = pro_z[str3, :]
t1 = Int32(cluster_cunitsset[3])
t2 = Int32(cluster_cunitsset[3]) + Int32(cluster_cunitsset[1])
ctem_pro[1:t1, :] = tem3
ctem_pro[t1+1:t2, :] = tem1
ctem_pro[t2+1:end, :] = tem2

y31 = Plots.heatmap(ctem_ben, c=cgrad([:White, :Crimson]), colorbar=false, title="pcuc w/o dynamic limit", fa=0.75)
y32 = Plots.heatmap(ctem_pro, c=cgrad([:White, :Crimson]), colorbar=false, title="pcuc w/ dynamic limit", fa=0.75)
Plots.plot(y31, y32, size=(800, 350), titlefontsize=7, layout=(1, 2))

unitscap = zeros(NT, 2)
for t in 1:NT
    unitscap[t, 1] = sum(ben_z[:, t] .* units.p_max)
    unitscap[t, 2] = sum(pro_z[:, t] .* units.p_max)
end
p16 = Plots.plot(unitscap[:, 1], lw=1, ylims=(0, 6), size=(400, 200), linetype=:steppre, lc=:darkorange, xlabel=L"t / h", xtickfontsize=6, ytickfontsize=6, ylabelfontsize=8, legendfontsize=6, ylabel="units capacites", label="w/o freq/flex limits", )
p16 = Plots.plot!(unitscap[:, 2], lw=1, label="w/   freq/flex limits", lc=:cornflowerblue,)
Plots.savefig(p16, filepath * "/com/figs/fig26.pdf", )

# NOTE final
# --------------------------------- save fig --------------------------------- #
# filepath = pwd()
