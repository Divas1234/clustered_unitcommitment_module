# windsimulation
using Distributions

struct wind
    index::Vector{Int64}
    locatebus::Vector{Int64}
    p_max::Vector{Float64}
    scenarios_prob::Float64
    scenarios_nums::Int64
    scenarios_curve::Array{Float64}
    # frequency constrol process of renewable energy
    Fcmode::Vector{Float64}
    # model-1(drop constrol):Kw/(Rw * (1 + sTw))
    Kw::Vector{Float64}
    Rw::Vector{Float64}
    # model-2(Virtual inertia mechine): (sMw + Dw)/(1 + sTw)
    Mw::Vector{Float64}
    Dw::Vector{Float64}
    Tw::Vector{Float64}
    # wind(Fcmode,Kw,Rw,Mw,Dw) = new(Fcmode,Kw,Rw,Mw,Dw)
    function checkvaildity(Fcmode)
        if Fcmode == 1
            if Kw != 0 || Rw != 0
                println("dismathch for Fcmode-1 and related params")
            end
        end
        return if Fcmode == 2
            if Kw != 0 || Rw != 0
                println("dismathch for Fcmode-2 and related params")
            end
        end
    end
    function wind(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve, Fcmode, Kw, Rw, Mw, Dw, Tw)
        return new(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve, Fcmode, Kw, Rw, Mw, Dw, Tw)
    end
end

struct solar
    index::Vector{Int64}
    locatebus::Vector{Int64}
    p_max::Vector{Float64}
    scenarios_prob::Float64
    scenarios_nums::Int64
    scenarios_curve::Array{Float64}
    # frequency constrol process of renewable energy
    Fcmode::Vector{Float64}
    # model-1(drop constrol):Kw/(Rw * (1 + sTw))
    Kw::Vector{Float64}
    Rw::Vector{Float64}
    # model-2(Virtual inertia mechine): (sMw + Dw)/(1 + sTw)
    Mw::Vector{Float64}
    Dw::Vector{Float64}
    Tw::Vector{Float64}
    # wind(Fcmode,Kw,Rw,Mw,Dw) = new(Fcmode,Kw,Rw,Mw,Dw)
    function checkvaildity(Fcmode)
        if Fcmode == 1
            if Kw != 0 || Rw != 0
                println("dismathch for Fcmode-1 and related params")
            end
        end
        return if Fcmode == 2
            if Kw != 0 || Rw != 0
                println("dismathch for Fcmode-2 and related params")
            end
        end
    end
    function solar(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve, Fcmode, Kw, Rw, Mw, Dw, Tw)
        return new(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve, Fcmode, Kw, Rw, Mw, Dw, Tw)
    end
end

# predefine value
index = [1; 2]
locatebus = [3; 4]
NW = length(index)

# NOTE = "NOTE: the number of scenarios should be equal to the number of scenarios in the wind and solar data"
scenarios_nums = 5

scenarios_prob = 1 / scenarios_nums
# scenarios_curve = zeros(scenarios_nums, NT)
# NT = 24
NT = 168
# assum the capacity of each wind is same and reforced as 0.5 p.u.
cap = [1.0] * 1.0
p_max = cap .* ones(NW, 1)
p_max = p_max[:, 1]
scenarios_curve = zeros(scenarios_nums, NT)


# scenarios_curvebase = reshape(scenarios_curvebase, 1, NT)

function winds_genscenario(WindsFreqParam, winds_scenarios_curvebase, pv_scenarios_curvebase, windspmax, pvpmax, flag)
    if flag == 1
        rand(123)
        scenarios_nums = 5
        winds_scenarioscurve = stochastizedrenewablecurve(winds_scenarios_curvebase, scenarios_nums, NT) .* windspmax
        pv_scenarioscurve = stochastizedrenewablecurve(pv_scenarios_curvebase, scenarios_nums, NT) .* pvpmax
        scenarios_curve = winds_scenarioscurve + pv_scenarioscurve * 0.0
        scenarios_prob = 1 / scenarios_nums
    else
        scenarios_nums = 1
        winds_scenarioscurve = winds_scenarios_curvebase .* windspmax
        pv_scenarioscurve = pv_scenarios_curvebase .* pvpmax
        scenarios_curve = winds_scenarioscurve + pv_scenarioscurve .* 0.0
        scenarios_prob = 1 / scenarios_nums
    end

    FCmode = WindsFreqParam[:, 1]
    KW = WindsFreqParam[:, 2]
    RW = WindsFreqParam[:, 3]
    MW = WindsFreqParam[:, 4]
    DW = WindsFreqParam[:, 5]
    TW = WindsFreqParam[:, 6]
    scenarios_nums = size(scenarios_curve, 1)
    winds = wind(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve, FCmode, KW, RW, MW, DW, TW)

    return winds, NW
end


function pv_genscenario(pvFreqParam, winds_scenarios_curvebase, pv_scenarios_curvebase, windspmax, pvpmax, flag)
    if flag == 1
        rand(123)
        scenarios_nums = 5
        winds_scenarioscurve = stochastizedrenewablecurve(winds_scenarios_curvebase, scenarios_nums, NT) .* windspmax
        pv_scenarioscurve = stochastizedrenewablecurve(pv_scenarios_curvebase, scenarios_nums, NT) .* pvpmax
        scenarios_curve = winds_scenarioscurve .* 0.0 + pv_scenarioscurve
        scenarios_prob = 1 / scenarios_nums
    else
        scenarios_nums = 1
        winds_scenarioscurve = winds_scenarios_curvebase .* windspmax
        pv_scenarioscurve = pv_scenarios_curvebase .* pvpmax
        scenarios_curve = winds_scenarioscurve .* 0.0 + pv_scenarioscurve
        scenarios_prob = 1 / scenarios_nums
    end

    FCmode = pvFreqParam[:, 1]
    KW = pvFreqParam[:, 2]
    RW = pvFreqParam[:, 3]
    MW = pvFreqParam[:, 4]
    DW = pvFreqParam[:, 5]
    TW = pvFreqParam[:, 6]
    scenarios_nums = size(scenarios_curve, 1)
    solarsation = solar(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve, FCmode, KW, RW, MW, DW, TW)

    return solarsation, NW
end

function stochastizedrenewablecurve(scenarios_curvebase, scenarios_nums, NT)
    sample_sets = rand(Weibull(), scenarios_nums * NT) * 0.0025
    scenarios_error = reshape(sample_sets, scenarios_nums, NT)
    for i = 1:scenarios_nums
        for j in 1:NT
            sample_temp = rand()
            if sample_temp > 0.5
                scenarios_curve[i, j] = scenarios_curvebase[1, j] + scenarios_error[i, j]
            else
                scenarios_curve[i, j] = scenarios_curvebase[1, j] - scenarios_error[i, j]
                if scenarios_curve[i, j] <= 0
                    scenarios_curve[i, j] = 0
                end
            end
        end
    end
    return scenarios_curve
end



# # windsimulation
# using Distributions

# struct wind
#     index::Vector{Int64}
#     locatebus::Vector{Int64}
#     p_max::Vector{Float64}
#     scenarios_prob::Float64
#     scenarios_nums::Int64
#     scenarios_curve::Array{Float64}
#     # frequency constrol process of renewable energy
#     Fcmode::Vector{Float64}
#     # model-1(drop constrol):Kw/(Rw * (1 + sTw))
#     Kw::Vector{Float64}
#     Rw::Vector{Float64}
#     # model-2(Virtual inertia mechine): (sMw + Dw)/(1 + sTw)
#     Mw::Vector{Float64}
#     Dw::Vector{Float64}
#     Tw::Vector{Float64}
#     # wind(Fcmode,Kw,Rw,Mw,Dw) = new(Fcmode,Kw,Rw,Mw,Dw)
#     function checkvaildity(Fcmode)
#         if Fcmode == 1
#             if Kw != 0 || Rw != 0
#                 println("dismathch for Fcmode-1 and related params")
#             end
#         end
#         return if Fcmode == 2
#             if Kw != 0 || Rw != 0
#                 println("dismathch for Fcmode-2 and related params")
#             end
#         end
#     end
#     function wind(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve, Fcmode, Kw, Rw, Mw, Dw, Tw)
#         return new(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve, Fcmode, Kw, Rw, Mw, Dw, Tw)
#     end
# end

# # predefine value
# index = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10]
# locatebus = [3; 4; 5; 3; 4; 5; 3; 4; 5; 3]
# NW = length(index)
# scenarios_nums = 5
# scenarios_prob = 1 / scenarios_nums
# # scenarios_curve = zeros(scenarios_nums, NT)
# # NT = 24
# NT = 168
# # assum the capacity of each wind is same and reforced as 0.5 p.u.
# cap = [0.10] * 1.0
# p_max = cap .* ones(NW, 1)
# p_max = p_max[:, 1]
# scenarios_curve = zeros(scenarios_nums, NT)


# # scenarios_curvebase = reshape(scenarios_curvebase, 1, NT)

# function genscenario(WindsFreqParam, winds_scenarios_curvebase, pv_scenarios_curvebase, windspmax, pvpmax, flag)
#     if flag == 1
#         rand(123)
#         scenarios_nums = 5
#         winds_scenarioscurve = stochastizedrenewablecurve(winds_scenarios_curvebase, scenarios_nums, NT) .* windspmax
#         pv_scenarioscurve = stochastizedrenewablecurve(pv_scenarios_curvebase, scenarios_nums, NT) .* pvpmax
#         scenarios_curve = winds_scenarioscurve + pv_scenarioscurve
#         scenarios_prob = 1 / scenarios_nums
#     else
#         scenarios_nums = 1
#         winds_scenarioscurve = winds_scenarios_curvebase .* windspmax
#         pv_scenarioscurve = pv_scenarios_curvebase .* pvpmax
#         scenarios_curve = winds_scenarioscurve + pv_scenarioscurve
#         scenarios_prob = 1 / scenarios_nums
#     end

#     FCmode = WindsFreqParam[:, 1]
#     KW = WindsFreqParam[:, 2]
#     RW = WindsFreqParam[:, 3]
#     MW = WindsFreqParam[:, 4]
#     DW = WindsFreqParam[:, 5]
#     TW = WindsFreqParam[:, 6]
#     scenarios_nums = size(scenarios_curve, 1)
#     winds = wind(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve, FCmode, KW, RW, MW, DW, TW)

#     return winds, NW
# end

# function stochastizedrenewablecurve(scenarios_curvebase, scenarios_nums, NT)
#     sample_sets = rand(Weibull(), scenarios_nums * NT) * 0.025
#     scenarios_error = reshape(sample_sets, scenarios_nums, NT)
#     for i = 1:scenarios_nums
#         for j in 1:NT
#             sample_temp = rand()
#             if sample_temp > 0.5
#                 scenarios_curve[i, j] = scenarios_curvebase[1, j] + scenarios_error[i, j]
#             else
#                 scenarios_curve[i, j] = scenarios_curvebase[1, j] - scenarios_error[i, j]
#             end
#         end
#     end
#     return scenarios_curve
# end