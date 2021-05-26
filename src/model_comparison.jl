using Plots
include("canonical_two_stage_model_with_division_functions.jl")
include("model_with_ribosomes_functions.jl")

function main()
    global n = 50
    global alpha = 7.0
    global gamma = 0.5
    global rho = 5.0
    T = range(1.5, 50, length=n)
    beta = range(0.01, 0.5, length=n)

    for rho in [0.1, 0.5, 1.0, 10.0]
        path = mkpath("data/model-comparison/rho-$rho")

        heatmap(T, beta, difference_in_means, c = :thermal, title = "Difference in Mean Protein Number")
        xlabel!("T")
        ylabel!("beta")
        savefig(string(path, "/protein_number_diff.svg"))

        heatmap(T, beta, percentage_difference_in_means, c = :thermal, title = "Percentage Difference in Mean Protein Number")
        xlabel!("T")
        ylabel!("beta")
        savefig(string(path, "/protein_number_percentage_diff.svg"))

        heatmap(T, beta, difference_in_variances, c = :thermal, title = "Difference in Protein Number Variance")
        xlabel!("T")
        ylabel!("beta")
        savefig(string(path, "/protein_variance_diff.svg"))

        heatmap(T, beta, fano_factor_model_1, c = :thermal, title = "Fano factor for Canonical Two Stage Model")
        xlabel!("T")
        ylabel!("beta")
        savefig(string(path, "/fano_factor_model_1.svg"))

        heatmap(T, beta, fano_factor_model_2, c = :thermal, title = "Fano Factor for Model Including Ribosome Kinetics")
        xlabel!("T")
        ylabel!("beta")
        savefig(string(path, "/fano_factor_model_2.svg"))

        heatmap(T, beta, fano_factor_difference, c = :thermal, title = "Difference in Fano Factor")
        xlabel!("T")
        ylabel!("beta")
        savefig(string(path, "/fano_factor_model_diff.svg"))

        heatmap(T, beta, protein_model_1, c = :thermal, title = "Protein Number for Canonical Two Stage Model")
        xlabel!("T")
        ylabel!("beta")
        savefig(string(path, "/protein_variance_model_1.svg"))

        heatmap(T, beta, protein_model_2, c = :thermal, title = "Protein Number for Model Including Ribosome Kinetics")
        xlabel!("T")
        ylabel!("beta")
        savefig(string(path, "/protein_variance_model_2.svg"))
    end
end

function protein_model_1(T, beta)
    CanonicalTwoStageModel.approximate_protein_mean(alpha, beta, gamma, T, 0)
end

function protein_model_2(T, beta)
    ModelWithRibosomes.p(alpha, model_2_beta(beta, T), gamma, rho, T, 0)
end

function protein_var_model_1(T, beta)
    CanonicalTwoStageModel.approximate_protein_variance(alpha, beta, gamma, T, 0)
end

function protein_var_model_2(T, beta)
    ModelWithRibosomes.p_var(alpha, model_2_beta(beta, T), gamma, rho, T, 0)
end

function difference_in_means(T, beta)
    return abs(ModelWithRibosomes.p(alpha, model_2_beta(beta, T), gamma, rho, T, 0) - CanonicalTwoStageModel.approximate_protein_mean(alpha, beta, gamma, T, 0))
end

function percentage_difference_in_means(T, beta)
    return abs((ModelWithRibosomes.p(alpha, model_2_beta(beta, T), gamma, rho, T, 0) - CanonicalTwoStageModel.approximate_protein_mean(alpha, beta, gamma, T, 0))/ModelWithRibosomes.p(alpha, model_2_beta(beta, T), gamma, rho, T, 0) )
end

function difference_in_variances(T, beta)
    return abs(ModelWithRibosomes.p_var(alpha, model_2_beta(beta, T), gamma, rho, T, 0) - CanonicalTwoStageModel.approximate_protein_variance(alpha, beta, gamma, T, 0))
end

function fano_factor_model_1(T, beta)
    return CanonicalTwoStageModel.approximate_protein_variance(alpha, beta, gamma, T, 0)/CanonicalTwoStageModel.approximate_protein_mean(alpha, beta, gamma, T, 0)
end

function fano_factor_model_2(T, beta)
    return ModelWithRibosomes.p_var(alpha, model_2_beta(beta, T), gamma, rho, T, 0)/ModelWithRibosomes.p(alpha, model_2_beta(beta, T), gamma, rho, T, 0)
end

function fano_factor_difference(T, beta)
    return fano_factor_model_1(T,beta) - fano_factor_model_2(T, beta)
end

function model_2_beta(beta, T)
    return 2/3*beta/(rho*T)
end

main()
