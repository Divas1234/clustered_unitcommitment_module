function calculate_frequencynadir(M, H, D, T, R, F, K, δp)
    ωₙ = (D * R + 1) / (2 * H * R * T)
    ι = (D * R * T + 2 * H * R + F * T) / (2 * (D * R + 1)) * ωₙ
    ωᵣ = ωₙ * sqrt(1 - ι^2)
    α = sqrt((1 - 2 * T * ι * ωₙ + T * T * ωᵣ * ωᵣ) / (1 - ι^2))
    θ = atan((ωᵣ * T) / (1 - ι * ωₙ * T)) - atan(sqrt(1 - ι^2) / ι * (-1))
    t_nadir = 1 / ωᵣ * atan(ωᵣ * T / (ι * ωᵣ * T - 1))
    Δf = R * δp / (D * R + 1) * (1 + sqrt(1 - ι^2) * α * exp(-1.0 * ι * ωₙ * t_nadir))
    f_nadir = 50 - abs(Δf)
    return f_nadir
end

# through cluster
function calculate_aggregatedfrequencyparameters(winds, cunits, NW, cNG, cluster_cunitsset, Sampling_Statue)
    # normalized winds parameters through COI
    vsmFC_number = sum(winds.Fcmode[:, 1])
    doopFC_number = length(winds.Fcmode[:, 1]) - vsmFC_number
    adjustablewindsVSCpower = winds.Fcmode .* winds.p_max * maximum(winds.scenarios_curve[1, :])
    # current_Kw    = sum(winds.Kw ./ winds.Rw .* (ones(NW, 1) - winds.Fcmode) .* winds.p_max) / sum(((ones(NW, 1) - winds.Fcmode) .* winds.p_max)) # Kw
    inverse_winds_Rw = zeros(NW, 1)
    for i in 1:NW
        if abs(winds.Fcmode[i, 1] - 1) >= 0.50
            inverse_winds_Rw[i, 1] = 1 / winds.Rw[i, 1]
        end
    end
    # @show inverse_winds_Rw

    # @show inverse_winds_Rw
    current_Rw = sum(winds.Kw .* inverse_winds_Rw .* (ones(NW, 1) - winds.Fcmode) .* winds.p_max) / sum(((ones(NW, 1) - winds.Fcmode) .* winds.p_max))
    current_Rw = 1 / current_Rw
    current_Dw = sum(winds.Dw .* adjustablewindsVSCpower) / sum(adjustablewindsVSCpower) # Dw
    current_Mw = sum(winds.Mw .* adjustablewindsVSCpower) / sum(adjustablewindsVSCpower) # Mw
    current_Hw = current_Mw / 2
    current_Kw = 1.0

    # units parameters
    adjustabletheramlpower = cunits.p_max .* Sampling_Statue
    current_Rg = sum(cunits.Kg ./ cunits.Rg .* adjustabletheramlpower) / sum(adjustabletheramlpower) # Kg
    current_Rg = 1 / current_Rg
    current_Tg = sum(cunits.Tg .* adjustabletheramlpower) / sum(adjustabletheramlpower) # T
    current_Fg_div_Rg = sum(cunits.Kg .* cunits.Fg ./ cunits.Rg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
    current_Fg = 1 / current_Fg_div_Rg
    current_Kg = 1
    current_Dg = sum(cunits.Dg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
    current_Hg = sum(cunits.Hg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
    current_Mg = current_Hg * 2

    # potential prob. # δp
    p_step = maximum(units.p_max) * 1.5

    # sumD and sumH
    current_sumD = (current_Dg + current_Dw) + (1 / current_Rw * 0.10)
    current_sumH = current_Hg + current_Hw
    current_sumM = current_sumH * 2

    # simplify parameters
    Mg, Hg, Dg, Tg, Rg, Fg, Kg = current_sumM, current_sumH, current_sumD, current_Tg, current_Rg, current_Fg, current_Kg
    M, H, D, T, R, F, K, δp = Mg, Hg, Dg, Tg, Rg, Fg, Kg, p_step
    return M, H, D, T, R, F, K, δp
end


# function frequencydynamic_ASFR(temflag)
#     Mg, Hg, Dg, Tg, Rg, Fg, Kg, δp, endtime = formparameter(temflag)
#     D, M, H, Tᵣ, Fₕ, R, K = Dg, Mg, Hg, Tg, Fg, Rg, Kg
#     R = Rg / Kg
#     δp = δp * 1.0
#     wₙ = sqrt((D + K / R) / (2 * H * Tᵣ))
#     ζ = (2 * H + (Fₕ * (K / R) + D) * Tᵣ) / (2 * sqrt(2 * H * Tᵣ) * sqrt(D + K / R))
#     # ζ = (D * R * Tᵣ + 2 * H * R + Fₕ * Tᵣ) / 2 / (D * R + 1) * wₙ
#     wᵣ = wₙ * sqrt(1 - ζ^2)
#     ψ = asin(sqrt(1 - ζ^2))

#     xdata = collect(0:δt:endtime)
#     ydata = zeros(size(xdata, 1), 1)
#     f_base = 50

#     newNT = size(xdata, 1)
#     for i in 1:newNT
#         t = xdata[i, 1]
#         δf = δp / (2 * H * Tᵣ * (wₙ^2))
#         δf =
#             δf +
#             δp / (2 * Hg * wᵣ) * exp(-ζ * wₙ * t) * (sin(wᵣ * t) - 1 / (wₙ * Tᵣ) * sin(wᵣ * t + ψ))
#         ydata[i, 1] = f_base - δf
#     end
#     ydata[1:2, 1] = ydata[1:2, 1] .* 1.00
#     ydata[:, 1] = ydata[:, 1] .+ 0.0
#     different_ASFR = zeros(size(xdata, 1), 1)
#     increment_Padd = zeros(size(xdata, 1), 1)
#     for i in 1:Int64(size(xdata, 1) - 1)
#         different_ASFR[i, 1] = (ydata[i+1, 1] - ydata[i, 1]) / δt
#         str = (Fg + (1 - Fg) / Tg) * exp(-1.0 / Tg * i * δt)
#         increment_Padd[i, 1] = (Kg / Rg) * str * (f_base - ydata[i, 1])
#     end
#     return ydata, different_ASFR, increment_Padd
# end