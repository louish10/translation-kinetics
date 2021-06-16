using Plots, LaTeXStrings
include("canonical_two_stage_model_with_division_functions.jl")
include("model_with_ribosomes_functions.jl")

function main()
    global n = 500
    global alpha = 3.42
    global gamma = 0.06
    global T = 27.5

    for i in [1.0, 10.0, 27.5, 50.0]
        rhoT = range(0.5*T, 100.0*T, length=n)
        betaT = range(0.5*T, 100.0*T, length=n)
        T = i
        path = mkpath("data/model-comparison/T-$T")

        heatmap(rhoT, betaT, difference_in_means, c = :thermal, title = "Relative Difference in Mean Protein Number")
        xlabel!(L"\rho T")
        ylabel!(L"\beta T")
        savefig(string(path, "/protein_number_diff.svg"))

        heatmap(rhoT, betaT, difference_in_variances, c = :thermal, title = "Relative Difference in Protein Number Variance")
        xlabel!(L"\rho T")
        ylabel!(L"\beta T")
        savefig(string(path, "/protein_variance_diff.svg"))

        heatmap(rhoT, betaT, fano_factor_model_1, c = :thermal, title = "Fano factor for Canonical Two Stage Model")
        xlabel!(L"\rho T")
        ylabel!(L"\beta T")
        savefig(string(path, "/fano_factor_model_1.svg"))

        heatmap(rhoT, betaT, fano_factor_model_2, c = :thermal, title = "Fano Factor for Model Including Ribosome Kinetics")
        xlabel!(L"\rho T")
        ylabel!(L"\beta T")
        savefig(string(path, "/fano_factor_model_2.svg"))

        heatmap(rhoT, betaT, fano_factor_difference, c = :thermal, title = "Relative Difference in Fano Factor")
        xlabel!(L"\rho T")
        ylabel!(L"\beta T")
        savefig(string(path, "/fano_factor_model_diff.svg"))

        heatmap(rhoT, betaT, protein_model_1, c = :thermal, title = "Protein Number for Canonical Two Stage Model")
        xlabel!(L"\rho T")
        ylabel!(L"\beta T")
        savefig(string(path, "/protein_number_model_1.svg"))

        heatmap(rhoT, betaT, protein_model_2, c = :thermal, title = "Protein Number for Model Including Ribosome Kinetics")
        xlabel!(L"\rho T")
        ylabel!(L"\beta T")
        savefig(string(path, "/protein_number_model_2.svg"))

        heatmap(rhoT, betaT, protein_var_model_1, c = :thermal, title = "Protein Variance for Canonical Two Stage Model")
        xlabel!(L"\rho T")
        ylabel!(L"\beta T")
        savefig(string(path, "/protein_variance_model_1.svg"))

        heatmap(rhoT, betaT, protein_var_model_2, c = :thermal, title = "Protein Variance for Model Including Ribosome Kinetics")
        xlabel!(L"\rho T")
        ylabel!(L"\beta T")
        savefig(string(path, "/protein_variance_model_2.svg"))
    end
end

function protein_model_1(rhoT, betaT)
    CanonicalTwoStageModel.p_time_av(alpha, betaT/T, gamma, T)
end

function protein_model_2(rhoT, betaT)
    ModelWithRibosomes.p_time_av(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T)
end

function protein_var_model_1(rhoT, betaT)
    CanonicalTwoStageModel.p_var_time_av(alpha, betaT/T, gamma, T)
end

function protein_var_model_2(rhoT, betaT)
    ModelWithRibosomes.p_var_time_av(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T)
end

function difference_in_means(rhoT, betaT)
    mod1 =CanonicalTwoStageModel.p_time_av(alpha, betaT/T, gamma, T)
    mod2 =ModelWithRibosomes.p_time_av(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T)
    return abs((mod1-mod2)/mod2) 
end

function difference_in_variances(rhoT, betaT)
    mod1 = CanonicalTwoStageModel.p_var_time_av(alpha, betaT/T, gamma, T)
    mod2 = ModelWithRibosomes.p_var_time_av(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T)
    return abs((mod1-mod2)/mod2) 
end

function fano_factor_model_1(rhoT, betaT)
    return CanonicalTwoStageModel.p_var_time_av(alpha, betaT/T, gamma, T)/CanonicalTwoStageModel.p_time_av(alpha, betaT/T, gamma, T)
end

function fano_factor_model_2(rhoT, betaT)
    return ModelWithRibosomes.p_var_time_av(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T)/ModelWithRibosomes.p_time_av(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T)
end

function fano_factor_difference(rhoT, betaT)
    mod1 = fano_factor_model_1(T,betaT/T)
    mod2 = fano_factor_model_2(T, model_2_beta(betaT/T, rhoT))
    return abs((mod1-mod2)/mod2) 
end

function model_2_beta(beta, rhoT)
    return 2/3*beta/(rhoT)
end

main()
