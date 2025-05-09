using Clustering, Statistics

function recognizing_critical_scenarios(winds, pvs, loads, NT)
    loadcurve = zeros(1, size(loads.load_curve, 2))
    for t in 1:size(loads.load_curve, 2)
        loadcurve[1, t] = sum(loads.load_curve[:, t]) / 1e0
    end
    windcurve = winds.scenarios_curve * 1.0
    pvcurve = pvs.scenarios_curve * 1.0
    NT, ND, NC = Int64(NT), 1, size(windcurve, 1)

    netloadcurves = zeros(NC, NT)
    for i in 1:NC
        netloadcurves[i, :] = loadcurve[1, 1:NT] .- windcurve[i, 1:NT] .- pvcurve[i, 1:NT]
    end
    # println(loadcurve)

    FREQ_indified_criticalscenarios = zeros(NC, NT)
    clustered_num = 6
    for i in 1:NC
        each_netloadcurve = netloadcurves[i, :]
        each_netloadcurve = each_netloadcurve[:, :]'
        FREQ_indified_criticalscenarios[i, :] = filtering_freq_criticalscenario(each_netloadcurve, clustered_num)
    end

    # return frequency_critical_scenario

    # NOTE recognizing critical scenarios for ramping-up dynamics
    seqential_rampingrate = zeros(NC, NT)
    for t in 2:NT
        seqential_rampingrate[:, t] .= netloadcurves[:, t] - netloadcurves[:, t-1]
        # filtering ramping-up scenarios_curve
        tem = findall(x -> x < 0, seqential_rampingrate[:, t])[:, :]
        seqential_rampingrate[tem, t] .= 0
    end

    # Plots.plot(seqential_rampingrate)
    FLEX_indified_criticalscenarios = zeros(NC, NT)
    for i in 1:NC
        each_rampingcurve = seqential_rampingrate[i, :]
        each_rampingcurve = each_rampingcurve[:, :]'
        FLEX_indified_criticalscenarios[i, :] = filtering_flex_criticalscenario(each_rampingcurve, clustered_num)
    end

    return FLEX_indified_criticalscenarios, FREQ_indified_criticalscenarios, seqential_rampingrate, netloadcurves
end


function filtering_freq_criticalscenario(each_netloadcurve, clustered_num)
    each_netloadcurve = each_netloadcurve[:, :]

    # println(each_netloadcurve)
    R = kmeans(each_netloadcurve, clustered_num, maxiter=2000, display=:none)
    @assert nclusters(R) == clustered_num # verify the number of clusters

    mapping_location = assignments(R) # get the assignments of points to clusters
    each_subclustered_num = counts(R) # get the cluster sizes
    M = adjoint((R.centers)) # get the cluster centers
    # println(M)
    # println(each_subclustered_num)
    # println(mapping_location)

    sorted_M = sort(M[:, 1])
    lc1 = findall(x -> x == sorted_M[1], M)[1][1]
    lc2 = findall(x -> x == sorted_M[2], M)[1][1]
    # println([lc1, lc2])
    critical_scenario_1 = findall(x -> x == lc1, mapping_location) # get the critical scenarios for frequency dynamic checkvaildity
    critical_scenario_2 = findall(x -> x == lc2, mapping_location)
    tem = zeros(1, size(each_netloadcurve, 2))
    # tem[1, Int64.(critical_scenario_1)] .= 1
    # tem[1, Int64.(critical_scenario_2)] .= 1

    critical_scenario_3 = findall(x -> x <= 0.750, each_netloadcurve[1, :])
    tem[1, Int64.(critical_scenario_3)] .= 1

    return tem

end


function filtering_flex_criticalscenario(each_netloadcurve, clustered_num)
    each_netloadcurve = each_netloadcurve[:, :]

    R = kmeans(each_netloadcurve, clustered_num, maxiter=2000, display=:none)
    @assert nclusters(R) == clustered_num # verify the number of clusters

    mapping_location = assignments(R) # get the assignments of points to clusters
    each_subclustered_num = counts(R) # get the cluster sizes
    M = adjoint((R.centers)) # get the cluster centers

    sorted_M = sort(M[:, 1], rev=true)
    lc1 = findall(x -> x == sorted_M[1], M)[1][1]
    lc2 = findall(x -> x == sorted_M[2], M)[1][1]
    critical_scenario_1 = findall(x -> x == lc1, mapping_location) # get the critical scenarios for frequency dynamic checkvaildity
    critical_scenario_2 = findall(x -> x == lc2, mapping_location)
    tem = zeros(1, size(each_netloadcurve, 2))
    # tem[1, Int64.(critical_scenario_1)] .= 1
    # tem[1, Int64.(critical_scenario_2)] .= 1

    critical_scenario_3 = findall(x -> x >= 0.50, each_netloadcurve[1, :])
    tem[1, Int64.(critical_scenario_3)] .= 1

    return tem

end
