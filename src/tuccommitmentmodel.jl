using JuMP, Gurobi, Test, DelimitedFiles
# using CPLEX
include("linearization.jl")
include("powerflowcalculation.jl")

function traditionalscucmodel(
    NT::Int64,
    NB::Int64,
    NG::Int64,
    ND::Int64,
    NC::Int64,
    units::unit,
    loads::load,
    winds::wind,
    pvs::solar,
    lines::transmissionline,
    config_param::config,
    rampingup_critical_scenario,
    frequency_critical_scenario
)
    println("Step-3: Creating dispatching model")

    if config_param.is_NetWorkCon == 1
        Adjacmatrix_BtoG, Adjacmatrix_B2D, Gsdf =
            linearpowerflow(units, lines, loads, NG, NB, ND, NL)
        Adjacmatrix_BtoW = zeros(NB, length(winds.index))
        for i in 1:length(winds.index)
            Adjacmatrix_BtoW[winds.index[i, 1], i] = 1
        end
    end

    NS = winds.scenarios_nums
    NW = length(winds.index)

    # creat scucsimulation_model

    # scuc = JuMP.Model(Gurobi.Optimizer)

    # scuc = JuMP.Model(Gurobi.Optimizer)
    scuc = JuMP.direct_model(Gurobi.Optimizer())
    # --------------------------------- setting 3 -------------------------------- #
    # set_optimizer_attribute(scuc, "Method", 2)
    # set_optimizer_attribute(scuc, "IterationLimit", 10)
    set_optimizer_attribute(scuc, "NoRelHeurTime", 50)
    set_optimizer_attribute(scuc, "NoRelHeurWork", 50)
    # set_optimizer_attribute(scuc, "Threads", 72)
    # set_optimizer_attribute(scuc, "IterationLimit", 100)
    # set_optimizer_attribute(scuc, "MIPGap", 0.01)
    # set_optimizer_attribute(scuc, "MIPFocus", 1)
    # set_optimizer_attribute(scuc, "FlowCoverCuts", 2)
    # set_optimizer_attribute(scuc, "MIRCuts", 2)
    # set_optimizer_attribute(scuc, "Presolve", -1)
    # set_optimizer_attribute(scuc, "Method", 4) # binfa
    # set_optimizer_attribute(scuc, "ImproveStartGap", 0.95)


    # env = Gurobi.Env(; start = false)
    # GRBsetintparam(env.ptr_env, GRB_INT_PAR_RECORD, 1)
    # GRBstartenv(env.ptr_env)
    # model = Model(() -> Gurobi.Optimizer(env))


    # scuc = JuMP.Model(CPLEX.Optimizer)
    # set_optimizer_attribute(scuc, "Method", 3)
    # set_optimizer_attribute(scuc, "IterationLimit", 1000)
    # set_optimizer_attribute(scuc, "NoRelHeurTime", 50)
    # set_optimizer_attribute(scuc, "NoRelHeurWork", 10)
    # set_optimizer_attribute(scuc, "MIPGap", 0.002)

    # NS = 1 # for test
    # --------------------------------- 2657.84 s -------------------------------- #
    # binary variables
    @variable(scuc, x[1:NG, 1:NT], Bin)
    @variable(scuc, u[1:NG, 1:NT], Bin)
    @variable(scuc, v[1:NG, 1:NT], Bin)

    # continuous variables
    @variable(scuc, pg₀[1:(NG*NS), 1:NT] >= 0)
    @variable(scuc, pgₖ[1:(NG*NS), 1:NT, 1:3] >= 0)
    @variable(scuc, su₀[1:NG, 1:NT] >= 0)
    @variable(scuc, sd₀[1:NG, 1:NT] >= 0)
    @variable(scuc, sr⁺[1:(NG*NS), 1:NT] >= 0)
    @variable(scuc, sr⁻[1:(NG*NS), 1:NT] >= 0)
    @variable(scuc, Δpd[1:(ND*NS), 1:NT] >= 0)
    @variable(scuc, Δpw[1:(NW*NS), 1:NT] >= 0)
    @variable(scuc, Δpr[1:(NW*NS), 1:NT] >= 0)
    # # pss variables
    # @variable(scuc, κ⁺[1:(NC * NS), 1:NT], Bin) # charge status
    # @variable(scuc, κ⁻[1:(NC * NS), 1:NT], Bin) # discharge status
    # @variable(scuc, pc⁺[1:(NC * NS), 1:NT] >= 0)# charge power
    # @variable(scuc, pc⁻[1:(NC * NS), 1:NT] >= 0)# discharge power
    # @variable(scuc, qc[1:(NC * NS), 1:NT] >= 0) # cumsum power
    # # @variable(scuc, pss_sumchargeenergy[1:NC * NS, 1] >= 0)

    # @variable(scuc, ζ[1:NG, 1:NT] >= 0)
    # @variable(scuc, z[1:(NG^2), 1:NT], Bin)

    refcost, eachslope = linearizationfuelcurve(units, NG)

    c₀ = config_param.is_CoalPrice
    pₛ = scenarios_prob
    plentycoffi_1 = config_param.is_LoadsCuttingCoefficient * 1e4
    plentycoffi_2 = config_param.is_WindsCuttingCoefficient * 1e4
    ρ⁺ = c₀ * 2
    ρ⁻ = c₀ * 2

    println("start...")
    println(
        "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
    )

    # model-1: MIQP with quartic fuel equation
    # @objective(scuc, Min, sum(sum(su₀[i, t] + sd₀[i, t] for i in 1:NG) for t in 1:NT)+
    #     pₛ*c₀/100*(sum(sum(sum(sum(units.coffi_a[i,1] * (pg₀[i + (s - 1) * NG,t]^2) + units.coffi_b[i,1] * pg₀[i + (s - 1) * NG,t] + units.coffi_c[i,1] for t in 1:NT)) for s in 1:NS) for i in 1:NG))+
    #     pₛ*plentycoffi_1*sum(sum(sum(Δpd[1+(s-1)*ND : s*ND, t]) for t in 1:NT) for s in 1:NS)+
    #     pₛ*plentycoffi_2*sum(sum(sum(Δpw[1+(s-1)*NW : s*NW, t]) for t in 1:NT) for s in 1:NS))

    # model-2:MILP with piece linearization equation of nonliear equation
    @objective(
        scuc,
        Min,
        100 * sum(sum(su₀[i, t] + sd₀[i, t] for i in 1:NG) for t in 1:NT) +
        pₛ *
        c₀ *
        (
            sum(
                sum(
                    sum(sum(pgₖ[i+(s-1)*NG, t, :] .* eachslope[:, i] for t in 1:NT)) for
                    s in 1:NS
                ) for i in 1:NG
            ) +
            sum(sum(sum(x[:, t] .* refcost[:, 1] for t in 1:NT)) for s in 1:NS) +
            sum(
                sum(
                    sum(
                        ρ⁺ * sr⁺[i+(s-1)*NG, t] + ρ⁻ * sr⁻[i+(s-1)*NG, t] for i in 1:NG
                    ) for t in 1:NT
                ) for s in 1:NS
            )
        ) +
        pₛ *
        plentycoffi_1 *
        sum(sum(sum(Δpd[(1+(s-1)*ND):(s*ND), t]) for t in 1:NT) for s in 1:NS) +
        pₛ *
        plentycoffi_2 *
        sum(sum(sum(Δpw[(1+(s-1)*NW):(s*NW), t]) for t in 1:NT) for s in 1:NS) +
        pₛ * plentycoffi_2 *
        sum(sum(sum(Δpr[(1+(s-1)*NW):(s*NW), t]) for t in 1:NT) for s in 1:NS)
    ) / 1e6

    #
    # for test
    # @objective(scuc, Min, 0)
    println("objective_function")
    println("\t MILP_type objective_function \t\t\t\t\t\t done")

    println("subject to.")

    # minimum shutup and shutdown ductration limits
    onoffinit, Lupmin, Ldownmin = zeros(NG, 1), zeros(NG, 1), zeros(NG, 1)
    for i in 1:NG
        onoffinit[i] = ((units.x_0[i, 1] > 0.5) ? 1 : 0)
        Lupmin[i] = min(NT, units.min_shutup_time[i] * onoffinit[i])
        Ldownmin[i] = min(NT, (units.min_shutdown_time[i, 1]) * (1 - onoffinit[i]))
    end

    # @constraint(scuc, [i = 1:NG, t = 1:Int64((Lupmin[i] + Ldownmin[i]))], x[i, t] == onoffinit[i])

    for i in 1:NG
        for t in Int64(max(1, Lupmin[i])):NT
            LB = Int64(max(t - units.min_shutup_time[i, 1] + 1, 1))
            @constraint(scuc, sum(u[i, r] for r in LB:t) <= x[i, t])
        end
        for t in Int64(max(1, Ldownmin[i])):NT
            LB = Int64(max(t - units.min_shutup_time[i, 1] + 1, 1))
            @constraint(scuc, sum(v[i, r] for r in LB:t) <= (1 - x[i, t]))
        end
    end
    println("\t constraints: 1) minimum shutup/shutdown time limits\t\t\t done")

    # binary variable logic
    @constraint(
        scuc,
        [i = 1:NG, t = 1:NT],
        u[i, t] - v[i, t] == x[i, t] - ((t == 1) ? onoffinit[i] : x[i, t-1])
    )
    @constraint(scuc, [i = 1:NG, t = 1:NT], u[i, t] + v[i, t] <= 1)
    println("\t constraints: 2) binary variable logic\t\t\t\t\t done")

    # shutup/shutdown cost
    shutupcost = units.coffi_cold_shutup_1
    shutdowncost = units.coffi_cold_shutdown_1
    @constraint(scuc, [t = 1], su₀[:, t] .>= shutupcost .* (x[:, t] - onoffinit[:, 1]))
    @constraint(scuc, [t = 1], sd₀[:, t] .>= shutdowncost .* (onoffinit[:, 1] - x[:, t]))
    @constraint(scuc, [t = 2:NT], su₀[:, t] .>= shutupcost .* u[:, t])
    @constraint(scuc, [t = 2:NT], sd₀[:, t] .>= shutdowncost .* v[:, t])
    println("\t constraints: 3) shutup/shutdown cost\t\t\t\t\t done")

    # loadcurtailments and spoliedwinds limits
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT],
        Δpw[(1+(s-1)*NW):(s*NW), t] .<= winds.scenarios_curve[s, t] * winds.p_max[:, 1]
    )
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT],
        Δpr[(1+(s-1)*NW):(s*NW), t] .<= pvs.scenarios_curve[s, t] * pvs.p_max[:, 1]
    )
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT],
        Δpd[(1+(s-1)*ND):(s*ND), t] .<= loads.load_curve[:, t]
    )
    # @constraint(scuc, [s=1:NS, t = 1:NT], Δpw[1+(s-1)*NW:s*NW, t] .== zeros(NW,1))
    # @constraint(scuc, [s=1:NS, t = 1:NT], Δpd[1+(s-1)*ND:s*ND, t] .== zeros(ND,1))
    println("\t constraints: 4) loadcurtailments and spoliedwinds\t\t\t done")

    # generatos power limits
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT],
        pg₀[(1+(s-1)*NG):(s*NG), t] + sr⁺[(1+(s-1)*NG):(s*NG), t] .<=
        units.p_max[:, 1] .* x[:, t]
    )
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT],
        pg₀[(1+(s-1)*NG):(s*NG), t] - sr⁻[(1+(s-1)*NG):(s*NG), t] .>=
        units.p_min[:, 1] .* x[:, t]
    )
    println("\t constraints: 5) generatos power limits\t\t\t\t\t done")

    # forcast_error = 0.05
    # forcast_reserve = winds.scenarios_curve * sum(winds.p_max[:, 1]) * forcast_error
    forcast_error = 0.025
    forcast_reserve = winds.scenarios_curve * sum(winds.p_max[:, 1]) * forcast_error + pvs.scenarios_curve * sum(pvs.p_max[:, 1]) * forcast_error
    # config_param.is_Alpha, config_param.is_Belta = 0, 0.025
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT, i = 1:NG],
        sum(sr⁺[(1+(s-1)*NG):(s*NG), t]) >=
        0.5 * (
            config_param.is_Alpha * forcast_reserve[s, t] +
            config_param.is_Belta * sum(loads.load_curve[:, t])
        )
    )
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT],
        sum(sr⁻[(1+(s-1)*NG):(s*NG), t]) >=
        0.5 * (
            config_param.is_Alpha * forcast_reserve[s, t] +
            config_param.is_Belta * sum(loads.load_curve[:, t])
        )
    )

    println("\t constraints: 6) system reserves limits\t\t\t\t\t done")

    # power balance constraints
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT],
        sum(pg₀[(1+(s-1)*NG):(s*NG), t]) +
        sum(winds.scenarios_curve[s, t] * winds.p_max[:, 1] - Δpw[(1+(s-1)*NW):(s*NW), t]) + sum(pvs.scenarios_curve[s, t] .* pvs.p_max[:, 1] - Δpr[(1+(s-1)*NW):(s*NW), t]) -
        sum(loads.load_curve[:, t] - Δpd[(1+(s-1)*ND):(s*ND), t]) .== 0
    )
    println("\t constraints: 7) power balance constraints\t\t\t\t done")

    # ramp-up and ramp-down constraints
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT],
        pg₀[(1+(s-1)*NG):(s*NG), t] -
        ((t == 1) ? units.p_0[:, 1] : pg₀[(1+(s-1)*NG):(s*NG), t-1]) .<=
        units.ramp_up[:, 1] .* ((t == 1) ? onoffinit[:, 1] : x[:, t-1]) +
        units.shut_up[:, 1] .* ((t == 1) ? ones(NG, 1) : u[:, t-1]) +
        units.p_max[:, 1] .* (ones(NG, 1) - ((t == 1) ? onoffinit[:, 1] : x[:, t-1]))
    )
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT],
        ((t == 1) ? units.p_0[:, 1] : pg₀[(1+(s-1)*NG):(s*NG), t-1]) -
        pg₀[(1+(s-1)*NG):(s*NG), t] .<=
        units.ramp_down[:, 1] .* x[:, t] +
        units.shut_down[:, 1] .* v[:, t] +
        units.p_max[:, 1] .* (x[:, t])
    )
    println("\t constraints: 8) ramp-up/ramp-down constraints\t\t\t\t done")

    # PWL constraints
    eachseqment = (units.p_max - units.p_min) / 3
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT, i = 1:NG],
        pg₀[i+(s-1)*NG, t] .== units.p_min[i, 1] * x[i, t] + sum(pgₖ[i+(s-1)*NG, t, :])
    )
    @constraint(
        scuc,
        [s = 1:NS, t = 1:NT, i = 1:NG, k = 1:3],
        pgₖ[i+(s-1)*NG, t, k] <= eachseqment[i, 1] * x[i, t]
    )
    println("\t constraints: 9) piece linearization constraints\t\t\t done")

# ----------------------------------------------------------------


    # transmissionline power limits for basline states
    # if config_param.is_NetWorkCon == 1
    #     for l in 1:NL
    #         subGsdf_units = Gsdf[l, units.locatebus]
    #         subGsdf_winds = Gsdf[l, winds.index]
    #         subGsdf_loads = Gsdf[l, loads.locatebus]
    #         subGsdf_psses = Gsdf[1, stroges.locatebus]
    #         @constraint(
    #             scuc,
    #             [s = 1:NS, t = 1:NT],
    #             sum(subGsdf_units[i] * pg₀[i + (s - 1) * NG, t] for i in 1:NG) + sum(
    #                 subGsdf_winds[w] *
    #                 (winds.scenarios_curve[s, t] * winds.p_max[w, 1] - Δpw[(s - 1) * NW + w, t]) for
    #                 w in 1:NW
    #             ) - sum(
    #                 subGsdf_loads[d] * (loads.load_curve[d, t] - Δpd[(s - 1) * ND + d, t]) for
    #                 d in 1:ND
    #             ) <= lines.p_max[l, 1]
    #         )
    #         @constraint(
    #             scuc,
    #             [s = 1:NS, t = 1:NT],
    #             sum(subGsdf_units[i] * pg₀[i + (s - 1) * NG, t] for i in 1:NG) + sum(
    #                 subGsdf_winds[w] *
    #                 (winds.scenarios_curve[s, t] * winds.p_max[w, 1] - Δpw[(s - 1) * NW + w, t]) for
    #                 w in 1:NW
    #             ) - sum(
    #                 subGsdf_loads[d] * (loads.load_curve[d, t] - Δpd[(s - 1) * ND + d, t]) for
    #                 d in 1:ND
    #             ) >= lines.p_min[l, 1]
    #         )
    #     end
    #     println("\t constraints: 10) transmissionline limits for basline\t\t\t done")
    # end


    # NOTE: this is not a complete implementation of flexibly-scheduled constraints
    flexbility_limit = 0.0

    # for s in 1:NS
    #     time_forflexcheck = findall(x -> x == 1, rampingup_critical_scenario[s, :])
    #     new_NT = size(time_forflexcheck, 1)
    #     @constraint(
    #         scuc,
    #         [t = 2:new_NT],
    #         sum(units.ramp_up[:, 1] .* ((time_forflexcheck[t] == 1) ? onoffinit[:, 1] : x[:, time_forflexcheck[t]-1])) +
    #         sum(units.shut_up[:, 1] .* ((time_forflexcheck[t] == 1) ? ones(NG, 1) : (x[:, time_forflexcheck[t]-1] - x[:, time_forflexcheck[t]-1]))) +
    #         sum(units.p_max[:, 1] .* (ones(NG, 1) - ((time_forflexcheck[t] == 1) ? onoffinit[:, 1] : x[:, time_forflexcheck[t]-1]))) >=
    #         sum(loads.load_curve[:, time_forflexcheck[t]] .- ((time_forflexcheck[t] == 1) ? loads.load_curve[:, time_forflexcheck[t]] * 0.75 : loads.load_curve[:, time_forflexcheck[t]-1])) * 0.50
    #     )
    #     @constraint(
    #         scuc,
    #         [t = 2:new_NT],
    #         sum(units.p_max[:, 1] .* ((time_forflexcheck[t] == 1) ? onoffinit[:, 1] : x[:, time_forflexcheck[t]-1])) -
    #         sum(pg₀[(1+(s-1)*NG):(s*NG), time_forflexcheck[t]]) >=
    #         sum(loads.load_curve[:, time_forflexcheck[t]] .- ((time_forflexcheck[t] == 1) ? loads.load_curve[:, time_forflexcheck[t]] * 0.75 : loads.load_curve[:, time_forflexcheck[t]-1])) * 0.50
    #     )
    # end
    # println("\t constraints: 12) flexibly-scheduled constraints\t\t\t done")

    # # frequency constrol process
    # f_base = 50.0
    # RoCoF_max = 2.0
    # f_nadir = 49.5
    # f_qss = 49.5
    # Δp = maximum(units.p_max[:, 1]) * 1.0

    # # RoCoF constraint
    # @constraint(
    #     scuc,
    #     [t = 1:NT],
    #     (sum(winds.Mw[:, 1] .* winds.Fcmode[:, 1] .* winds.p_max[:, 1]) + 2 * sum(units.Hg[:, 1] .* units.p_max[:, 1] .* x[:, t]))
    #     /
    #     (sum(units.p_max[:, 1]) + sum(winds.Fcmode .* winds.p_max))
    #     >=
    #     Δp * f_base / RoCoF_max / 50.0 * 0.80
    # )

    # f_nadir constraint
    # |---H---|---D---|---F---|---K---|---δp---|
    # MODEL = 1
    # if MODEL == 1
    #     # A, ϵ = creatfrequencyfittingfunction(units, winds, NG, NW)
    #     param_H, param_D, param_F, param_K, param_δ = -0.08, 128.14, 15.69, -12.635, 239.625
    #     ϵ = -472.0816
    #     coff₁ =
    #         (
    #             param_H * units.Hg +
    #             param_D * units.Dg +
    #             param_F * units.Kg .* units.Fg ./ units.Rg +
    #             param_K * units.Kg ./ units.Rg
    #         ) .* units.p_max
    #     coff₂ = (f_base - f_nadir) * Δp * f_base - (param_δ * Δp + ϵ[1, 1])
    #     # println(coff₁)
    #     # println(coff₂)
    # else
    #     coff₁ = [
    #         2.780932120714308e17,
    #         3.482699236849364e17,
    #         3.359168850810112e17,
    #         3.997594260256427e17,
    #         1.6866162819599613e17,
    #     ]
    #     coff₂ = 3.787335736962308e17
    # end
    # @constraint(scuc, [t = 1:NT], sum(coff₁ .* x[:, t]) <= sum(coff₂ * units.p_max .* x[:, t]))

    println("\t constraints: 12) frequency responce constraints limits\t\t\t done")

    println("\n")
    println("Model has been loaded")
    println("Step-4: calculation...")

    JuMP.optimize!(scuc)
    # JuMP.write_to_file(scuc, "issue477.mps")


    println("callback gurobisolver\t\t\t\t\t\t\t done")
    @test JuMP.termination_status(scuc) == MOI.OPTIMAL
    @test JuMP.primal_status(scuc) == MOI.FEASIBLE_POINT
    println("#TEST: termination_status\t\t\t\t\t\t pass")

    println(
        "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
    )
    # println("pg₀>>\n", JuMP.value.(pg₀[1:NG,1:NT]))
    # println("x>>  \n", JuMP.value.(x))
    # println("Δpd>>\n", JuMP.value.(Δpd[1:ND,1:NT]))
    # println("Δpw>>\n", JuMP.value.(Δpw[1:NW,1:NT]))
    su_cost = sum(JuMP.value.(su₀))
    sd_cost = sum(JuMP.value.(sd₀))
    pᵪ = JuMP.value.(pgₖ)
    p₀ = JuMP.value.(pg₀)
    x₀ = JuMP.value.(x)
    r⁺ = JuMP.value.(sr⁺)
    r⁻ = JuMP.value.(sr⁻)
    pᵨ = JuMP.value.(Δpd)
    pᵩ = JuMP.value.(Δpw)
    qᵩ = JuMP.value.(Δpr)
    # pss_charge_state⁺ = JuMP.value.(κ⁺)
    # pss_charge_state⁻ = JuMP.value.(κ⁻)
    # pss_charge_p⁺ = JuMP.value.(pc⁺)
    # pss_charge_p⁻ = JuMP.value.(pc⁻)
    # pss_charge_q = JuMP.value.(qc)

    prod_cost =
        pₛ *
        c₀ *
        (
            sum(
                sum(
                    sum(sum(pᵪ[i+(s-1)*NG, t, :] .* eachslope[:, i] for t in 1:NT)) for
                    s in 1:NS
                ) for i in 1:NG
            ) + sum(sum(sum(x₀[:, t] .* refcost[:, 1] for t in 1:NT)) for s in 1:NS)
        )
    cr⁺ =
        pₛ *
        c₀ *
        sum(sum(sum(ρ⁺ * r⁺[i+(s-1)*NG, t] for i in 1:NG) for t in 1:NT) for s in 1:NS)
    cr⁻ =
        pₛ *
        c₀ *
        sum(sum(sum(ρ⁺ * r⁻[i+(s-1)*NG, t] for i in 1:NG) for t in 1:NT) for s in 1:NS)
    seq_sr⁺ = pₛ * c₀ * sum(ρ⁺ * r⁺[i, :] for i in 1:NG)
    seq_sr⁻ = pₛ * c₀ * sum(ρ⁺ * r⁻[i, :] for i in 1:NG)
    𝜟pd = pₛ * sum(sum(sum(pᵨ[(1+(s-1)*ND):(s*ND), t]) for t in 1:NT) for s in 1:NS)
    𝜟pw = pₛ * sum(sum(sum(pᵩ[(1+(s-1)*NW):(s*NW), t]) for t in 1:NT) for s in 1:NS)
    str = zeros(1, 7)
    str[1, 1] = su_cost * 1
    str[1, 2] = sd_cost * 1
    str[1, 3] = prod_cost
    str[1, 4] = cr⁺
    str[1, 5] = cr⁻
    str[1, 6] = 𝜟pd
    str[1, 7] = 𝜟pw

    filepath = pwd()
    open(
        "D:/GithubClonefiles/clustered_unitcommitment_module/res/ben_calculation_result.txt",
        "w",
    ) do io
        writedlm(io, [" "])
        writedlm(io, ["su_cost" "sd_cost" "prod_cost" "cr⁺" "cr⁻" "𝜟pd" "𝜟pw"], '\t')
        writedlm(io, str, '\t')
        writedlm(io, [" "])
        writedlm(io, ["list 1: units stutup/down states"])
        writedlm(io, JuMP.value.(x), '\t')
        writedlm(io, [" "])
        writedlm(io, ["list 2: units dispatching power in scenario NO.1"])
        writedlm(io, JuMP.value.(pg₀[1:NG, 1:NT]), '\t')
        writedlm(io, [" "])
        writedlm(io, ["list 3: spolied wind power"])
        writedlm(io, JuMP.value.(Δpw[1:NW, 1:NT]), '\t')
        writedlm(io, [" "])
        writedlm(io, ["list 4: forced load curtailments"])
        writedlm(io, JuMP.value.(Δpd[1:ND, 1:NT]), '\t')
        writedlm(io, ["list 5: spolied pvs power"])
        writedlm(io, JuMP.value.(Δpr[1:NW, 1:NT]), '\t')
        # writedlm(io, [" "])
        # writedlm(io, ["list 5: pss charge state"])
        # writedlm(io, pss_charge_state⁺[1:NC, 1:NT], '\t')
        # writedlm(io, [" "])
        # writedlm(io, ["list 6: pss discharge state"])
        # writedlm(io, pss_charge_state⁻[1:NC, 1:NT], '\t')
        # writedlm(io, [" "])
        # writedlm(io, ["list 7: pss charge power"])
        # writedlm(io, pss_charge_p⁺[1:NC, 1:NT], '\t')
        # writedlm(io, [" "])
        # writedlm(io, ["list 8: pss discharge power"])
        # writedlm(io, pss_charge_p⁻[1:NC, 1:NT], '\t')
        # writedlm(io, [" "])
        # writedlm(io, ["list 9: pss strored energy"])
        # writedlm(io, pss_charge_q[1:NC, 1:NT], '\t')
        writedlm(io, [" "])
        writedlm(io, ["list 10: sr⁺"])
        writedlm(io, r⁺[1:NG, 1:NT], '\t')
        writedlm(io, [" "])
        writedlm(io, ["list 11: sr⁻"])
        return writedlm(io, r⁻[1:NG, 1:NT], '\t')
    end

    println("the calculation_result has been saved into | calculation_result.txt |\t done")
    println(
        "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
    )

    # return p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, pss_charge_p⁺, pss_charge_p⁻, su_cost, sd_cost, prod_cost, cr⁺, cr⁻
    return p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cr⁺, cr⁻
end
