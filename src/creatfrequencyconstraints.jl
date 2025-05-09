using MultivariateStats, LinearRegression
# include("D:/OneDriveFles/OneDrive/.uestc/AcademicResearchWorks/Prevoious Code/task 9/master-10 (little case - tuc vs cuc)/src/calculatefrequencynadirvalue.jl")

function creatfrequencyfittingfunction(cunits, winds, cNG, NW, cluster_cunitsset)

    Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Rg, unitsamplestatues = montecalrosimulation(cunits, winds, cNG, NW, cluster_cunitsset)
    # x: |---H---|---D---|---F---|---K---|
    # y: f_nadir
    y = Set_f_nadir
    x = zeros(size(set_H, 1), 3)
    str = size(set_H, 1)
    for i in 1:str
        x[i, 1:3] = [set_H[i, 1], set_Fg[i, 1], set_Rg[i, 1]]
    end

    # clearing the obtained data
    # x₀ = filter(!isnan, x)
    # y₀ = filter(!isnan, y)
    x₀ = hcat((filter(!isnan, x[:, i]) for i in 1:size(x, 2))...)
    y₀ = hcat((filter(!isnan, y[:, i]) for i in 1:size(y, 2))...)

    # delate the duplicate elements
    # x₁ = Rmduplicateelements(x₀, 1)
    # y₁ = Rmduplicateelements(y₀, 1)

    # println(x), println(y)
    # sol = llsq(x₀, y₀)
    # A, b = sol[1:end-1, :], sol[end, :]

    # println(x), println(y)
    lr = linregress(x₀, y₀)
    A = LinearRegression.slope(lr)[:,:]
    b = LinearRegression.bias(lr)[:,:]

    return A, b

end

function montecalrosimulation(cunits, winds, cNG, NW, cluster_cunitsset)
    sampleNumber = 10000
    unitsamplestatues = zeros(cNG, sampleNumber)
    for i in 1:sampleNumber
        for g in 1:cNG
            unitsamplestatues[g, i] = rand(1:cluster_cunitsset[g])
            unitsamplestatues[3, i] = 0
        end
    end

    # just for a test
    # i = 5
    # Sampling_Statue = unitsamplestatues[:, i]
    # M, H, D, T, R, F, K, δp = calculate_aggregatedfrequencyparameters(winds, cunits, NW, cNG, Sampling_Statue)
    # f_nadir = calculate_frequencynadir(M, H, D, T, R, F, K, δp)

    Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Kg, set_Rg = zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1)

    for i in 1:sampleNumber
        Sampling_Statue = unitsamplestatues[:, i]
        M, H, D, T, R, F, K, δp, Hg, Dg, Rg, Fg = calculate_aggregatedfrequencyparameters(winds, cunits, NW, cNG, cluster_cunitsset, Sampling_Statue)
        f_nadir = calculate_frequencynadir(M, H, D, T, R, F, K, δp)
        Set_f_nadir[i, 1], set_H[i, 1], set_δp[i, 1], set_Dg[i, 1], set_Fg[i, 1], set_Kg[i, 1], set_Rg[i, 1] = f_nadir, Hg, δp, Dg, Fg, K, Rg
    end

    return Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Rg, unitsamplestatues

end

function calculate_frequencynadir(M, H, D, T, R, F, K, δp)
    ωₙ = (D * R + 1) / (2 * H * R * T)
    ι = (D * R * T + 2 * H * R + F * T) / (2 * (D * R + 1)) * ωₙ
    ωᵣ = ωₙ * sqrt(1 - ι^2)
    α = sqrt((1 - 2 * T * ι * ωₙ + T * T * ωᵣ * ωᵣ) / (1 - ι^2))
    θ = atan((ωᵣ * T) / (1 - ι * ωₙ * T)) - atan(sqrt(1 - ι^2) / ι * (-1))
    t_nadir = 1 / ωᵣ * atan(ωᵣ * T / (ι * ωᵣ * T - 1))
    Δf = R * δp / (D * R + 1) * (1 + sqrt(1 - ι^2) * α * exp(-1.0 * ι * ωₙ * t_nadir))
    # f_nadir = 50 - abs(Δf)
    # tem = abs(Δf)
    return Δf
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
    return M, H, D, T, R, F, K, δp, current_Hg, current_Dg, 1 / current_Rg, current_Fg_div_Rg
end
