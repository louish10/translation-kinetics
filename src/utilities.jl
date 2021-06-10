module Utilities
    using Statistics
    function print_summary(path, n_simulations, proteins, mrnas, alpha, beta, gamma, delta)
        protein_mean = mean(proteins)
        protein_var = var(proteins)
        protein_mean_theory = beta*alpha/(delta*gamma)
        protein_var_theory = beta*alpha/(delta*gamma)*(1+beta/(delta+gamma))

        mrna_mean = mean(mrnas)
        mrna_var = var(mrnas)
        mrna_mean_theory = alpha/gamma
        mrna_var_theory = alpha/gamma

        protein_mrna_cov = cov(proteins, mrnas)
        protein_mrna_cov_theory = beta*alpha/(gamma*(delta+gamma))

        filename = string(path, "/results.txt")
        touch(filename)
        open(filename, "w") do io
            write(io, "Number of simulations: $n_simulations\n")
            write(io, "-----------\n")

            write(io, "Proteins\n")
            write(io, "Mean\n")
            write(io, "Theory: $protein_mean_theory\n")
            write(io, "Simulation: $protein_mean\n")
            write(io, "\nVariance\n")
            write(io, "Theory: $protein_var_theory\n")
            write(io, "Simulation: $protein_var\n")

            write(io, "-----------\n")

            write(io, "mRNAs\n")
            write(io, "Mean\n")
            write(io, "Theory: $mrna_mean_theory\n")
            write(io, "Simulation: $mrna_mean\n")
            write(io, "\nVariance\n")
            write(io, "Theory: $mrna_var_theory\n")
            write(io, "Simulation: $mrna_var\n")

            write(io, "-----------\n")

            write(io, "Covariance\n")
            write(io, "Theory: $protein_mrna_cov_theory\n")
            write(io, "Simulation: $protein_mrna_cov\n")
        end
    end

    export print_summary
end
