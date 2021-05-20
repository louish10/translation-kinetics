using Plots, Statistics, Distributions, SpecialFunctions

function main(args)
    global alpha = parse(Float64, args[1])
    global beta = parse(Float64, args[2])
    global gamma = parse(Float64, args[3])
    global delta = parse(Float64, args[4])
    global rho = parse(Float64, args[5])
    global omega = parse(Float64, args[5])
    global T = parse(Float64, args[6])

    global path = mkpath("data/model-with-ribosomes/alph-$alpha-bet-$beta-gam-$gamma-del-$delta-rho-$rho-T-$T")

    m_before, m_after, p_before, p_after, r_before, r_after = simulate_network()

    generate_summary(m_before, m_after, p_before, p_after, r_before, r_after)
    create_mrna_plots(m_before, m_after)
    create_protein_plots(p_before, p_after)
    create_ribosome_plots(r_before, r_after)
end

function generate_summary(m_before, m_after, p_before, p_after, r_before, r_after)
    # mRNA
    m_zero = mrna(0)
    m_T = mrna(T)

    mrna_after_stat_mean = mean(m_after)
    mrna_before_stat_mean = mean(m_before)

    mrna_after_stat_var = var(m_after)
    mrna_before_stat_var = var(m_before)

    mrna_after_theory_var = mrna_var(0)
    mrna_before_theory_var = mrna_var(T)

    # Protein
    p_before_stat_mean = mean(p_before)
    p_after_stat_mean = mean(p_after)

    p_zero = p(0)
    p_T = p(T)

    p_after_stat_var = var(p_after)
    p_after_theory_var = p_var(0)

    p_before_theory_var = p_var(T)
    p_before_stat_var = var(p_before)


    # Ribosome
    r_before_stat_mean = mean(r_before)
    r_after_stat_mean = mean(r_after)

    r_zero = r(0)
    r_T = r(T)

    r_after_stat_var = var(r_after)
    r_after_theory_var = r_var(0)

    r_before_theory_var = r_var(T)
    r_before_stat_var = var(r_before)

    # Covariances
    m_p_stat_covar_before = cov(m_before, p_before)
    m_p_stat_covar_after = cov(m_after, p_after)

    m_p_theory_covar_before = m_p_covar(T)
    m_p_theory_covar_after = m_p_covar(0)

    m_r_stat_covar_before = cov(m_before, r_before)
    m_r_stat_covar_after = cov(m_after, r_after)

    m_r_theory_covar_before = m_r_covar(T)
    m_r_theory_covar_after = m_r_covar(0)

    p_r_stat_covar_before = cov(p_before, r_before)
    p_r_stat_covar_after = cov(p_after, r_after)

    p_r_theory_covar_before = p_r_covar(T)
    p_r_theory_covar_after = p_r_covar(0)


    filename = string(path, "/results.txt")
    touch(filename)
    
    open(filename, "w") do io
        write(io, "After division\n")
        write(io, "-----------\n")

        write(io, "mRNA\n\n")
        write(io, "Mean\n")
        write(io, "Theory: $m_zero\n")
        write(io, "Simulation: $mrna_after_stat_mean\n")
        write(io, "\nVariance\n")
        write(io, "Theory: $mrna_after_theory_var\n")
        write(io, "Simulation: $mrna_after_stat_var\n")

        write(io, "-----------\n")

        write(io, "Protein\n\n")
        write(io, "Mean\n")
        write(io, "Theory: $p_zero\n")
        write(io, "Simulation: $p_after_stat_mean\n")
        write(io, "\nVariance\n")
        write(io, "Theory: $p_after_theory_var\n")
        write(io, "Simulation: $p_after_stat_var\n")

        write(io, "-----------\n")

        write(io, "Ribosome\n\n")
        write(io, "Mean\n")
        write(io, "Theory: $r_zero\n")
        write(io, "Simulation: $r_after_stat_mean\n")
        write(io, "\nVariance\n")
        write(io, "Theory: $r_after_theory_var\n")
        write(io, "Simulation: $r_after_stat_var\n")

        write(io, "-----------\n")
        write(io, "Covariances\n\n")
        write(io, "mRNA, protein\n")
        write(io, "Theory: $m_p_theory_covar_after\n")
        write(io, "Simulation: $m_p_stat_covar_after\n")
        write(io, "mRNA, ribosome\n")
        write(io, "Theory: $m_r_theory_covar_after\n")
        write(io, "Simulation: $m_r_stat_covar_after\n")
        write(io, "protein, ribosome\n")
        write(io, "Theory: $p_r_theory_covar_after\n")
        write(io, "Simulation: $p_r_stat_covar_after\n")

        write(io, "-----------\n\n")
        write(io, "Before division\n")
        write(io, "-----------\n")

        write(io, "mRNA\n\n")
        write(io, "Mean\n")
        write(io, "Theory: $m_T\n")
        write(io, "Simulation: $mrna_before_stat_mean\n")
        write(io, "\nVariance\n")
        write(io, "Theory: $mrna_before_theory_var\n")
        write(io, "Simulation: $mrna_before_stat_var\n")

        write(io, "-----------\n")

        write(io, "Protein\n\n")
        write(io, "Mean\n")
        write(io, "Theory: $p_T\n")
        write(io, "Simulation: $p_before_stat_mean\n")
        write(io, "\nVariance\n")
        write(io, "Theory: $p_before_theory_var\n")
        write(io, "Simulation: $p_before_stat_var\n")

        write(io, "-----------\n")

        write(io, "Ribosome\n\n")
        write(io, "Mean\n")
        write(io, "Theory: $r_T\n")
        write(io, "Simulation: $r_before_stat_mean\n")
        write(io, "\nVariance\n")
        write(io, "Theory: $r_before_theory_var\n")
        write(io, "Simulation: $r_before_stat_var\n")

        write(io, "-----------\n")
        write(io, "Covariances\n\n")
        write(io, "mRNA, protein\n")
        write(io, "Theory: $m_p_theory_covar_before\n")
        write(io, "Simulation: $m_p_stat_covar_before\n")
        write(io, "mRNA, ribosome\n")
        write(io, "Theory: $m_r_theory_covar_before\n")
        write(io, "Simulation: $m_r_stat_covar_before\n")
        write(io, "protein, ribosome\n")
        write(io, "Theory: $p_r_theory_covar_before\n")
        write(io, "Simulation: $p_r_stat_covar_before\n")

    end
end

function create_mrna_plots(m_before, m_after)
    mean_before = mrna(T)
    std_before = sqrt(mrna_var(T))

    nbins = maximum([1, Int(maximum(m_before))])
    x=collect(0:.1:nbins)
    histogram(m_before, nbins=Int(nbins), normed=true, linecolor=:match)
    plot!(x, gaussian(x, mean_before, std_before), label="gaussian fit", lw=3)
    plot!(x, poisson(x, mean_before), label="poisson fit", lw=3)
    xlabel!("m")
    ylabel!("P(m)")
    savefig(string(path, "/mrnas_before.svg"))

    mean_after = mrna(0)
    std_after = sqrt(mrna_var(0))

    nbins = maximum([1, Int(maximum(m_after))])
    x=collect(0:.1:nbins)
    histogram(m_after, nbins=Int(nbins), normed=true, linecolor=:match)
    plot!(x, gaussian(x, mean_after, std_after), label="gaussian fit", lw=3)
    plot!(x, poisson(x, mean_after), label="poisson fit", lw=3)
    xlabel!("m")
    ylabel!("P(m)")
    savefig(string(path, "/mrnas_after.svg"))
end

function create_ribosome_plots(r_before, r_after)
    mean_before = r(T)
    std_before = sqrt(r_var(T))

    nbins = maximum([1, Int(maximum(r_before))])
    x=collect(0:.1:nbins)
    histogram(r_before, nbins=Int(nbins), normed=true, linecolor=:match)
    plot!(x, gaussian(x, mean_before, std_before), label="gaussian fit", lw=3)
    plot!(x, poisson(x, mean_before), label="poisson fit", lw=3)
    xlabel!("r")
    ylabel!("P(r)")
    savefig(string(path, "/ribosomes_before.svg"))

    mean_after = r(0)
    std_after = sqrt(r_var(0))

    nbins = maximum([1, Int(maximum(r_after))])
    x=collect(0:.1:nbins)
    histogram(r_after, nbins=Int(nbins), normed=true, linecolor=:match)
    plot!(x, gaussian(x, mean_after, std_after), label="gaussian fit", lw=3)
    plot!(x, poisson(x, mean_after), label="poisson fit", lw=3)
    xlabel!("r")
    ylabel!("P(r)")
    savefig(string(path, "/ribosomes_after.svg"))
end

function create_protein_plots(p_before, p_after)
    mean_before = p(T)
    std_before = sqrt(p_var(T))

    nbins = maximum([1, Int(maximum(p_before))])
    x=collect(0:.1:nbins)
    histogram(p_before, nbins=Int(nbins), normed=true, linecolor=:match)
    plot!(x, gaussian(x, mean_before, std_before), label="gaussian fit", lw=3)
    plot!(x, negative_binomial(x, mean_before, std_before), label="negative binomial fit", lw=3)
    xlabel!("p")
    ylabel!("P(p)")
    savefig(string(path, "/proteins_before.svg"))

    mean_after = p(0)
    std_after = sqrt(p_var(0))

    nbins = maximum([1, Int(maximum(p_after))])
    x=collect(0:.1:nbins)
    histogram(p_after, nbins=Int(nbins), normed=true, linecolor=:match)
    plot!(x, gaussian(x, mean_after, std_after), label="gaussian fit", lw=3)
    plot!(x, negative_binomial(x, mean_after, std_after), label="negative binomial fit", lw=3)
    xlabel!("p")
    ylabel!("P(p)")
    savefig(string(path, "/proteins_after.svg"))
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

function n_choose_k(n, k)

end

function simulate_network()
    write(stdout, "Starting simulation")
    # This is heuristic but seems right
    t_final = T*100000
    time_since_division = 0
    t = 0
    p_after_division = []
    p_before_division = []
    m_after_division = []
    m_before_division = []
    r_after_division = []
    r_before_division = []
    
    species = [0, 0, 0]
    propensities = [
        species -> alpha,
        species -> beta * species[1] * species[3],
        species -> gamma * species[1],
        species -> delta * species[2],
        species -> rho
    ]

    write(stdout, "\n")
    while true

        if t%1000 == 0
            write(stdout, "\r$t of $t_final")
        end

        instantaneous_propensities = map(f -> f(species), propensities)
        a = sum(instantaneous_propensities)

        dt = 1/a*log(1/(1-rand()))
        t = t + dt
        time_since_division += dt
        if t > t_final
            break
        end

        if time_since_division > T
            t = t - t%T
            time_since_division = 0

            if t > t_final/1000
                append!(m_before_division, species[1])
                append!(p_before_division, species[2])
                append!(r_before_division, species[3])
            end

            species = partition_species(species)

            if t > t_final/1000
                append!(m_after_division, species[1])
                append!(p_after_division, species[2])
                append!(r_after_division, species[3])
            end
        else
            R = rand()
            if R < instantaneous_propensities[1]/a
                species[1] += 1
            elseif R < sum(instantaneous_propensities[1:2])/ a
                species[2] += 1
            elseif R < sum(instantaneous_propensities[1:3])/ a
                species[1] -= 1
            elseif R < sum(instantaneous_propensities[1:4])/a
                species[2] -= 1
            else
                species[3] += 1
            end
        end
    end

    return m_before_division, m_after_division, p_before_division, p_after_division, r_before_division, r_after_division
end

function mrna(t)
    return alpha/gamma*(1-exp(-gamma * t)/(2-exp(-gamma*T)))
end

function p(t)
    return alpha*beta*(2+4*T*gamma - (t^2 + 2*t*T + 3*T^2)*gamma^2 + 2*exp((T-t)*gamma)*(1+(t+T)*gamma)+2*exp(T*gamma)*(-2+gamma*(t^2*gamma+T*(-2+2*t*gamma+3*T*gamma))))*rho / (2 * (-1 + 2*exp(T*gamma))*gamma^3)
end

function r(t)
    return rho*(t+T)
end

function mrna_var(t)
    return alpha*(exp(gamma*(T-t)) - 2*exp(gamma*T) + 1)/(gamma*(1-2*exp(gamma*T)))
end

function r_var(t)
    return rho*(t+T)
end

function p_var(t)
    return (alpha*beta*rho*exp(-gamma*t))/(18*gamma^5*(1-2*exp(gamma*T))^2)*(exp(gamma*t)*(alpha*beta*(2*gamma*(gamma^2*(9*t^2*T+3*t^3+9*t*T^2+7*T^3)-4*gamma*T*(3*t+4*T)-6*t+8*T)+25)-2*beta*rho*(gamma*(3*gamma*t^2*(3-2*gamma*t)+6*gamma*T^2*(7-3*gamma*t)-18*T*(gamma*t-2)*(gamma*t+1)+4*gamma^2*T^3)+18)+9*gamma^2*(gamma*(gamma*t^2+T*(2*gamma*t+3*gamma*T-4))-2))-4*exp(gamma*(t+T))*(beta*rho*(gamma*(6*gamma*t^2*(2*gamma*t-3)+9*gamma*T^2*(4*gamma*t-7)+36*T*(gamma*t-2)*(gamma*t+1)+10*gamma^2*T^3)-54)+2*alpha*beta*(3*gamma^3*t^3+9*gamma^3*t^2*T+3*gamma*t*(gamma*T*(3*gamma*T-4)-3)+gamma*T*(gamma*T*(7*gamma*T-13)+3)+13)+9*gamma^2*(gamma*(gamma*t^2+T*(2*gamma*t+3*gamma*T-3))-2))+4*exp(gamma*(t+2*T))*(alpha*beta*(2*gamma*(gamma^2*(9*t^2*T+3*t^3+9*t*T^2+7*T^3)-2*gamma*T*(6*t+5*T)-12*t+T)+31)+2*beta*rho*(gamma*(3*gamma*t^2*(2*gamma*t-3)+3*gamma*T^2*(6*gamma*t-7)+18*T*(gamma*t-2)*(gamma*t+1)+14*gamma^2*T^3)-36)+9*gamma^2*(gamma*(gamma*t^2+T*(2*gamma*t+3*gamma*T-2))-2))+12*exp(2*gamma*T)*(18*beta*rho+3*gamma*(gamma+beta*gamma*rho*(3*t^2+6*t*T+2*T^2)+gamma^2*(t+T)*(beta*rho*t*(t+2*T)+1)+6*beta*rho*(t+T))+alpha*beta*(3*gamma^2*(t+T)^2-4*gamma*T-10))-6*exp(gamma*T) *(18*beta*rho+3*gamma *(gamma+beta*gamma*rho*(3*t^2+6*t*T+2*T^2)+gamma^2*(t+T)*(beta*rho* t*(t+2*T)+1)+6*beta*rho*(t+T))+alpha*beta*(3*gamma^2*(t+T)^2-4*gamma*T-8))+9*alpha*beta*exp(-(gamma*(t-2*T)))*(2*gamma*(t+T)+3))
end

function m_r_covar(t)
    return 0
end

function p_r_covar(t)
   return alpha*beta*rho/(6*gamma^3*(2*exp(gamma*T)-1))*(-3*gamma^2*(t+T)^2+2*exp(gamma*T)*(3*gamma^2*(t+T)^2 - 4*gamma*T - 4) + 6*exp(gamma*(T-t))*(gamma*(t+T)+1)+4*gamma*T+2) 
end

function m_p_covar(t)
    return exp(-t*gamma)*alpha*beta*rho*(-exp(gamma*T)*t*(t+2*T)*gamma^2+2*exp(gamma*t)*(gamma*(t+T)-1)*(2*exp(gamma*T)-1))/(2*(2*exp(T*gamma)-1)*gamma^3)
end

function partition_species(species)
    dist1 = Binomial(species[1])
    dist2 = Binomial(species[2])
    dist3 = Binomial(species[3])
    return [rand(dist1), rand(dist2), rand(dist3)]
end

main(ARGS)
