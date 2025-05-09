using Plots
include("/home/yuanyiping/下载/task 9/master-5/tem/unitsol.jl")
include("/home/yuanyiping/下载/task 9/master-5/src/formatteddata.jl")
include("/home/yuanyiping/下载/task 9/master-5/src/renewableenergysimulation.jl")
include("/home/yuanyiping/下载/task 9/master-5/src/readdatafromexcel.jl")
include("/home/yuanyiping/下载/task 9/master-5/src/cluster_units.jl")
include("/home/yuanyiping/下载/task 9/master-5/src/recognizingcriticalscenarios.jl")
UnitsFreqParam, WindsFreqParam, StrogeData, DataGen, GenCost, DataBranch, LoadCurve, DataLoad = readxlssheet()

function checkflexinfo()

    config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC = forminputdata(DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)
    pro_unitsschedulingpower_withflex, pro_unitsschedulingpower_withoutflex, pro_unitsonline_withflex, pro_unitsonline_withoutflex = creatingflexbilitycheckresults()
    cunits, cNG, cluster_cunitsset, cluster_featurematrix = calculating_clustered_units(units, DataGen, GenCost, UnitsFreqParam)
    winds, NW = genscenario(WindsFreqParam, 1)

    NT = 24 * 1
    rampingup_critical_scenario, frequency_critical_scenario = recognizing_critical_scenarios(winds, loads, NT)
    p1 = Plots.heatmap(rampingup_critical_scenario, c=cgrad([:white, :red]),
        # title = "critical scenarios for flexility-check",
        ylabel="scenarios", xlabel="time")
    p2 = Plots.heatmap(frequency_critical_scenario, c=cgrad([:white, :blue]),
        #  title="critical scenarios for frequency-dynamics",
        ylabel="scenarios", xlabel="time")
    p3 = Plots.plot(p1, p2;
        size=(800, 300),
        titlefontsize=8,
        layout=(1, 2)
    )

    # rampingup_critical_scenario


    pro1, pro2, pro0, statisticnumber = calculate_FREQandFLEX()
    pro_unitsschedulingpower_withflex, pro_unitsschedulingpower_withoutflex, pro_unitsonline_withflex, pro_unitsonline_withoutflex = creatingflexbilitycheckresults()

    # cunits, cNG, cluster_cunitsset, cluster_featurematrix = calculating_clustered_units(units, DataGen, GenCost, UnitsFreqParam)
    onoffinit, Lupmin, Ldownmin = zeros(NG, 1), zeros(NG, 1), zeros(NG, 1)
    for i in 1:cNG
        onoffinit[i] = ((units.x_0[i, 1] > 0.5) ? 1 : 0)
        Lupmin[i] = min(NT, units.min_shutup_time[i] * onoffinit[i])
        Ldownmin[i] = min(NT, (units.min_shutdown_time[i, 1]) * (1 - onoffinit[i]))
    end

    flexsol_1, flexsol_0 = zeros(NT, 3), zeros(NT, 3)
    NT = 24
    for t in 1:NT
        a = sum(units.ramp_up[:, 1] .* ((t == 1) ? onoffinit[:, 1] : pro2[:, t-1]))
        b = sum(units.shut_up[:, 1] .* ((t == 1) ? ones(NG, 1) : (pro2[:, t-1] - pro2[:, t-1])))
        c = sum(units.p_max[:, 1] .* (ones(NG, 1) - ((t == 1) ? onoffinit[:, 1] : pro2[:, t-1])))
        flexsol_1[t, 1] = a + b + c
        a = sum(units.ramp_up[:, 1] .* ((t == 1) ? onoffinit[:, 1] : pro0[:, t-1]))
        b = sum(units.shut_up[:, 1] .* ((t == 1) ? ones(NG, 1) : (pro0[:, t-1] - pro0[:, t-1])))
        c = sum(units.p_max[:, 1] .* (ones(NG, 1) - ((t == 1) ? onoffinit[:, 1] : pro0[:, t-1])))
        flexsol_0[t, 1] = a + b + c
    end

    for t in 1:NT
        d = sum(units.p_max[:, 1] .* ((t == 1) ? onoffinit[:, 1] : pro2[:, t])) - sum(pro_unitsschedulingpower_withflex[:, t])
        flexsol_1[t, 2] = d
        d = sum(units.p_max[:, 1] .* ((t == 1) ? onoffinit[:, 1] : pro0[:, t])) - sum(pro_unitsschedulingpower_withoutflex[:, t])
        flexsol_0[t, 2] = d
    end
    flexsol_1[1, 2] = flexsol_1[2, 2] * 0.80
    flexsol_1[1, 1] = flexsol_1[2, 1] * 0.80
    flexsol_0[1, 2] = flexsol_0[2, 2] * 0.80
    flexsol_0[1, 1] = flexsol_0[2, 1] * 0.80
    for t in 1:NT
        flexsol_1[t, 3] = min(flexsol_1[t, 2], flexsol_1[t, 1])
        flexsol_0[t, 3] = min(flexsol_0[t, 2], flexsol_0[t, 1])
    end

    # Plots.plot(flexsol_1[1:24, 1])
    # Plots.plot!(flexsol_0[1:24, 1])
    # Plots.plot!(flexsol_1[1:24, 2])
    # Plots.plot!(flexsol_0[1:24, 2])
    Plots.plot(flexsol_1[1:24, 3])
    Plots.plot!(flexsol_0[1:24, 3])
    # # Plots.plot!(flexsol[:, 3])

    loadcurve = zeros(1, min(size(loads.load_curve, 2), NT))
    for t in 1:min(size(loads.load_curve, 2), NT)
        loadcurve[1, t] = sum(loads.load_curve[:, t])
    end

    tem = zeros(NT, 1)
    # loadcurve = loadcurve - reshape(winds.scenarios_curve[1, :] * sum(winds.p_max), 1, 24)
    loadcurve
    # reshape(winds.scenarios_curve[1, :] * sum(winds.p_max), 1, 24)
    for t in 1:NT
        tem[t, 1] = sum(loadcurve[:, t] .- ((t == 1) ? loadcurve[:, t] * 1.0 : loadcurve[:, t-1])) * 1.00
    end
    tem
    Plots.plot(tem)
    # Plots.plot!(rampingup_critical_scenario[1, :])
    # rampingup_critical_scenario[1, :]

    # windcurve = winds.scenarios_curve * sum(winds.p_max)
    # NT, ND, NC = Int64(NT), 1, size(windcurve, 1)

    # netloadcurves = zeros(NC, NT)
    # for i in 1:NC
    #     netloadcurves[i, :] = loadcurve[1, 1:NT] .- windcurve[i, 1:NT]
    # end
    # Plots.plot(loadcurve')
    # Plots.plot!(netloadcurves[1, :])


    for t in 1:NT
        flexsol_1[t, 3] = max(flexsol_1[t, 3], 0)
        flexsol_0[t, 3] = max(flexsol_0[t, 3], 0)
        tem[t, 1] = max(tem[t, 1], 0)
    end


    f1 = Plots.plot(flexsol_1[1:24, 3], lw=1, marker=:circle, markersize=3, lc=:blue, la=0.5, size=(300, 200), label="pcuc w/ dynamic limit", title = "flexbility check")
    f1 = Plots.plot!(tem, lw=1, marker=:diamond, markersize=3, lc=:black, la=0.5, label="dynamic limit")
    f1 = Plots.plot!(flexsol_0[1:24, 3], lw=1, marker=:diamond, markersize=3, lc=:orange, la=0.5, label="pcuc w/o dynamic limit")
    # f2 = Plots.plot!(tem, lw=1, marker=:diamond, markersize=6, lc=:green, la=0.5)

    return p3, f1

end