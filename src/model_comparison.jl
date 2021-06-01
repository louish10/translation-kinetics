using Plots, LaTeXStrings
include("canonical_two_stage_model_with_division_functions.jl")
include("model_with_ribosomes_functions.jl")

function main()
    global n = 500
    global alpha = 7.0
    global gamma = 0.5
    global T = 1.0
    rhoT = range(0.1, 50.0, length=n)
    betaT = range(0.1, 50.0, length=n)

    for T in [0.1, 0.5, 1.0, 10.0]
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
        savefig(string(path, "/protein_variance_model_1.svg"))

        heatmap(rhoT, betaT, protein_model_2, c = :thermal, title = "Protein Number for Model Including Ribosome Kinetics")
        xlabel!(L"\rho T")
        ylabel!(L"\beta T")
        savefig(string(path, "/protein_variance_model_2.svg"))
    end
end

function protein_model_1(rhoT, betaT)
    CanonicalTwoStageModel.approximate_protein_mean(alpha, betaT/T, gamma, T, 0)
end

function protein_model_2(rhoT, betaT)
    ModelWithRibosomes.p(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T, 0)
end

function protein_var_model_1(rhoT, betaT)
    CanonicalTwoStageModel.approximate_protein_variance(alpha, betaT, gamma, T, 0)
end

function protein_var_model_2(rhoT, betaT)
    ModelWithRibosomes.p_var(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T, 0)
end

function difference_in_means(rhoT, betaT)
    mod1 =CanonicalTwoStageModel.approximate_protein_mean(alpha, betaT/T, gamma, T, 0)
    mod2 =ModelWithRibosomes.p(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T, 0)
    return abs((mod1-mod2)/mod2) 
end

function difference_in_variances(rhoT, betaT)
    mod1 = CanonicalTwoStageModel.approximate_protein_variance(alpha, betaT/T, gamma, T, 0)
    mod2 = ModelWithRibosomes.p_var(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T, 0)
    return abs((mod1-mod2)/mod2) 
end

function fano_factor_model_1(rhoT, betaT)
    return CanonicalTwoStageModel.approximate_protein_variance(alpha, betaT/T, gamma, T, 0)/CanonicalTwoStageModel.approximate_protein_mean(alpha, betaT/T, gamma, T, 0)
end

function fano_factor_model_2(rhoT, betaT)
    return ModelWithRibosomes.p_var(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T, 0)/ModelWithRibosomes.p(alpha, model_2_beta(betaT/T, rhoT), gamma, rhoT/T, T, 0)
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
