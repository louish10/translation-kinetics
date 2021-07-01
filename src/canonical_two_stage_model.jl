using Plots, Statistics, SpecialFunctions
include("simple_canonical_model_functions.jl")

function main(args)
    alpha = parse(Float64, args[1])
    beta = parse(Float64, args[2])
    gamma = parse(Float64, args[3])
    delta = parse(Float64, args[4])
    n_simulations = parse(Int, args[5])

    proteins = zeros(n_simulations)
    mrnas = zeros(n_simulations)
    path = mkpath("data/canonical-two-stage-model/alph-$alpha-bet-$beta-gam-$gamma-del-$delta")

    for i in 1:n_simulations
        m, p = simulate_network(alpha, beta, gamma, delta)
        proteins[i] = p
        mrnas[i] = m

        if i % 1000 == 0
            write(stdout, "\r$i")
        end
    end
    write(stdout, "\n")
    nbins = maximum([1, Int(maximum(proteins))])
    histogram(proteins, nbins=Int(nbins), normed=true)
    xlabel!("p")
    ylabel!("P(p)")
    x=collect(0:.1:nbins)
    protein_mean = beta*alpha/(delta*gamma)
    protein_var = beta*alpha/(delta*gamma)*(1+beta/(delta+gamma))

    plot!(x, gaussian(x, protein_mean, sqrt(protein_var)), label="gaussian fit", lw=3)
    plot!(x, negative_binomial(x, protein_mean, sqrt(protein_var)), label="negative binomial fit", lw=3)
    savefig(string(path, "/translation_proteins.svg"))

    mrna_mean = alpha/gamma
    mrna_std = sqrt(alpha/gamma)
    nbins = Int(maximum(mrnas))
    histogram(mrnas, nbins=Int(nbins), normed=true)
    vline!([alpha/gamma], label="theoretical mean")
    vline!([mean(mrnas)], label="statistical mean")
    xlabel!("m")
    ylabel!("P(m)")
    plot!(x, gaussian(x, mrna_mean, mrna_std), label="gaussian fit", lw=3)
    plot!(x, poisson(x, protein_mean), label="poisson fit", lw=3)
    savefig(string(path, "/translation_mrnas.svg"))

    print_summary(path, n_simulations, proteins, mrnas, alpha, beta, gamma, delta)
end

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


function simulate_network(alpha, beta, gamma, delta)
    # This is heuristic but seems right
    t_final = 1/(delta+gamma)*100
    t = 0
    m = 0
    p = 0
    
    species = [0, 0]
    propensities = [
        species -> alpha,
        species -> beta * species[1],
        species -> gamma * species[1],
        species -> delta * species[2],
    ]

    while true
        instantaneous_propensities = map(f -> f(species), propensities)
        a = sum(instantaneous_propensities)

        prop_1 = alpha
        prop_2 = beta * m
        prop_3 = gamma * m
        prop_4 = delta * p

        dt = 1/a*log(1/(1-rand()))
        t = t + dt
        if t > t_final
            break
        end

        R = rand()
        if R < instantaneous_propensities[1]/a
            species[1] += 1
        elseif R < sum(instantaneous_propensities[1:2])/ a
            species[2] += 1
        elseif R < sum(instantaneous_propensities[1:3])/ a
            species[1] -= 1
        else
            species[2] -= 1
        end
    end
    return species[1], species[2]
end

function gaussian(x, mean, std)
    1/(std*sqrt(2*pi))*exp.(-1/2*(x.-mean).^2 ./(std^2))
end

function negative_binomial(x, mean, std)
    p = 1-mean/std^2
    r = mean^2/(std^2-mean)
    return SpecialFunctions.gamma.(x .+ r) ./ (SpecialFunctions.gamma.(x .+ 1) .* SpecialFunctions.gamma.(r)) .* (1-p) .^ r .* p .^x
end

function poisson(x, lambda)
    lambda .^ x .* exp(-1*lambda) ./ SpecialFunctions.gamma.(x .+ 1)
end
main(ARGS)
