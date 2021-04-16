using Plots, Statistics, Distributions
include("utilities.jl")

function main(args)
    global alpha = parse(Float64, args[1])
    global beta = parse(Float64, args[2])
    global gamma = parse(Float64, args[3])
    global delta = parse(Float64, args[4])
    global T = parse(Float64, args[5])

    global path = mkpath("data/canonical-two-stage-model-with-division/alph-$alpha-bet-$beta-gam-$gamma-del-$delta-T-$T")

    m_before, m_after, p_before, p_after = simulate_network()

    write(stdout, "\n")

    m_zero = mrna(0)
    m_T = mrna(T)
    A = alpha/gamma - m_zero

    nbins = maximum([1, Int(maximum(m_after))])
    histogram(m_after, nbins=Int(nbins), normed=true)
    xlabel!("m")
    ylabel!("P(m)")
    m_zero = alpha/gamma*(1-exp(-gamma*T))/(2-exp(-gamma*T))
    vline!([m_zero], label="theoretical mean")
    vline!([mean(m_after)], label="statistical mean")
    savefig(string(path, "/mrna_after_division.svg"))

    mrna_stat_mean = mean(m_after)
    stat_var = var(m_after)
    theory_var = alpha/gamma - A

    write(stdout, "\nStatistical mean of mrna after division: $mrna_stat_mean")
    write(stdout, "\nTheoretical mean of mrna after division: $m_zero")
    write(stdout, "\nStatistical variance of mrna after division: $stat_var")
    write(stdout, "\nTheoretical variance of mrna after division: $theory_var")

    statistical_covariance = cov(m_after, p_after)
    cov_int_const = 1/(4-exp(-(gamma+delta)*T))*(-3*alpha/(gamma*(gamma+delta)) + (4-exp(-gamma*T))*(A/delta))
    theoretical_covariance = beta*(alpha/(gamma*(gamma+delta)) - A/delta+cov_int_const)
    write(stdout, "\nStatistical covariance after division: $statistical_covariance")
    write(stdout, "\nTheoretical covariance after division: $theoretical_covariance")

    nbins = maximum([1, Int(maximum(m_before))])
    histogram(m_before, nbins=Int(nbins), normed=true)
    xlabel!("m")
    ylabel!("P(m)")
    m_T = alpha/gamma*(1-exp(-gamma*T)) + m_zero*exp(-gamma*T)
    vline!([m_T], label="theoretical mean")
    vline!([mean(m_before)], label="statistical mean")
    savefig(string(path, "/mrna_before_division.svg"))

    nbins = maximum([1, Int(maximum(p_after))])
    histogram(p_after, nbins=Int(nbins), normed=true)
    xlabel!("p")
    ylabel!("P(p)")
    A = alpha/gamma - m_zero
    p_zero = 1/(2-exp(-delta*T))*(beta*alpha/(delta*gamma)-A*beta/(delta-gamma)*exp(-gamma*T)-(beta*alpha/(gamma*delta)+A*beta/(gamma-delta))*exp(-delta*T))
    vline!([p_zero], label="theoretical mean")
    vline!([mean(p_after)], label="statistical mean")
    savefig(string(path, "/protein_after_division.svg"))

    nbins = maximum([1, Int(maximum(p_before))])
    histogram(p_before, nbins=Int(nbins), normed=true)
    xlabel!("p")
    ylabel!("P(p)")
    p_T = beta*alpha/(delta*gamma) - A*beta/(delta-gamma)*exp(-gamma*T) + (p_zero - beta*alpha/(gamma*delta)-A*beta/(gamma-delta))*exp(-delta*T)
    vline!([p_T], label="theoretical mean")
    vline!([mean(p_before)], label="statistical mean")
    savefig(string(path, "/protein_before_division.svg"))

    protein_stat_var = var(p_after)
    protein_int_const = 1/(4-exp(-2*delta*T))*(F(cov_int_const, T) - 4*F(cov_int_const, 0) + p_T)
    protein_theory_var = F(cov_int_const, 0) + protein_int_const

    print_summary(p_before, p_after, m_before, m_after)
end

function print_summary(p_before, p_after, m_before, m_after)
    mrna_after_stat_mean = mean(m_after)
    m_zero = alpha/gamma*(1-exp(-gamma*T))/(2-exp(-gamma*T))
    A = alpha/gamma - m_zero
    mrna_after_stat_var = var(m_after)
    mrna_after_theory_var = alpha/gamma - A 

    covar_after_stat = cov(m_after, p_after)
    c_prime = 1/(4-exp(-(gamma+delta)*T))*(-3*alpha/(gamma*(gamma+delta)) + (4-exp(-gamma*T))*(A/delta))
    covar_after_theory = beta*(alpha/(gamma*(gamma+delta)) - A/delta+c_prime)
    covar_before_theory = beta*(alpha/(gamma*(gamma+delta)) - A/delta*exp(-gamma*T)+ c_prime*exp(-(gamma+delta)*T))
    covar_before_stat = cov(m_before, p_before)

    p_after_stat_mean = mean(p_after)
    p_zero = 1/(2-exp(-delta*T))*(beta*alpha/(delta*gamma)-A*beta/(delta-gamma)*exp(-gamma*T)-(beta*alpha/(gamma*delta)+A*beta/(gamma-delta))*exp(-delta*T))
    p_T = beta*alpha/(delta*gamma) - A*beta/(delta-gamma)*exp(-gamma*T) + (p_zero - beta*alpha/(gamma*delta)-A*beta/(gamma-delta))*exp(-delta*T)
    p_after_stat_var = var(p_after)
    c_double_prime = 1/(4-exp(-2*delta*T))*(F(c_prime, T) - 4*F(c_prime, 0) + p_T)
    p_after_theory_var = F(c_prime, 0) + c_double_prime

    mrna_before_stat_mean = mean(m_before)
    m_T = alpha/gamma*(1-exp(-gamma*T)) + m_zero*exp(-gamma*T)

    mrna_before_theory_var = alpha/gamma - A*exp(-gamma*T)
    mrna_before_stat_var = var(m_before)

    filename = string(path, "/results.txt")
    touch(filename)
    open(filename, "w") do io
        write(io, "After division")
        write(io, "-----------\n")

        write(io, "mRNA\n")
        write(io, "Mean\n")
        write(io, "Theory: $m_zero\n")
        write(io, "Simulation: $mrna_after_stat_mean\n")
        write(io, "\nVariance\n")
        write(io, "Theory: $mrna_after_theory_var\n")
        write(io, "Simulation: $mrna_after_stat_var\n")

        write(io, "-----------\n")

        write(io, "Protein\n")
        write(io, "Mean\n")
        write(io, "Theory: $p_zero\n")
        write(io, "Simulation: $p_after_stat_mean\n")
        write(io, "\nVariance\n")
        write(io, "Theory: $p_after_stat_var\n")
        write(io, "Simulation: $p_after_stat_var\n")

        write(io, "-----------\n")

        write(io, "Covariance\n")
        write(io, "Theory: $covar_after_stat\n")
        write(io, "Simulation: $covar_after_theory\n")

        write(io, "-----------\n")
        write(io, "Before division\n")
        write(io, "-----------\n")

        write(io, "mRNA\n")
        write(io, "Mean\n")
        write(io, "Theory: $m_T\n")
        write(io, "Simulation: $mrna_before_stat_mean\n")
        write(io, "\nVariance\n")
        write(io, "Theory: $mrna_before_theory_var\n")
        write(io, "Simulation: $mrna_before_stat_var\n")

        write(io, "Covariance\n")
        write(io, "Theory: $covar_before_stat\n")
        write(io, "Simulation: $covar_before_theory\n")

        approx_mrna_var_zero = mrna_variance(0)
        approx_mrna_var_T = mrna_variance(T)
        approx_covar_T = approximate_covariance(T)
        approx_covar_zero = approximate_covariance(0)
        c_double_prime_approx = alpha*beta/gamma * (beta/gamma*T-5*beta/(6*gamma^2)+T+1/gamma)
        approx_protein_var_zero = approximate_protein_variance(0)
        approx_protein_mean_zero = approximate_protein_mean(0)
        approx_protein_mean_T = approximate_protein_mean(T)

        write(io, "-----------\n")
        write(io, "Approximations\n")
        write(io, "mRNA variance after division: $approx_mrna_var_zero\n")
        write(io, "mRNA variance before division: $approx_mrna_var_T\n")
        write(io, "covariance before division: $approx_covar_T\n")
        write(io, "covariance after division: $approx_covar_zero\n")
        write(io, "protein variance after division: $approx_protein_var_zero\n")
        write(io, "protein mean before division: $approx_protein_mean_T\n")
        write(io, "protein mean after division: $approx_protein_mean_zero\n")
    end
end

function F(c_prime, t)
    m_zero = alpha/gamma*(1-exp(-gamma*T))/(2-exp(-gamma*T))
    A = alpha/gamma - m_zero
    p_zero = 1/(2-exp(-delta*T))*(beta*alpha/(delta*gamma)-A*beta/(delta-gamma)*exp(-gamma*T)-(beta*alpha/(gamma*delta)+A*beta/(gamma-delta))*exp(-delta*T))
    B = p_zero - beta*alpha/(gamma*delta) - A*beta/(gamma-delta)
    
    return beta^2*(alpha/(gamma*delta*(gamma+delta))-2*A/(delta*(2*delta - gamma))*exp(-gamma*t) + 2*c_prime/(delta-gamma)*exp(-(gamma+delta)*t)) + beta*alpha/(gamma*delta) - A*beta/(2*delta-gamma)*exp(-gamma*t) + beta*delta*A/((gamma-delta)*(2*delta-gamma))*exp(-gamma*t) + B*exp(-delta*t)
end

function mrna(t)
    return alpha/beta*(1-exp(-gamma * t)/(2-exp(-gamma*T)))
end

function mrna_variance(t)
    return alpha/gamma*(1-exp(-gamma* t)/(2-exp(-gamma * T)))
end

function approximate_covariance(t)
    return beta*alpha/gamma*(1/gamma - t*exp(-gamma*t)/(2-exp(-gamma*T)) - exp(-gamma * t)/(4-exp(-gamma*T))*(3/gamma + T*exp(-gamma*T)/(2-exp(-gamma*T))))
end

function approximate_protein_mean(t)
    return alpha*beta/gamma*(t + T + 1/gamma*(1+exp(gamma*(T-t))-exp(gamma*T))/(2-exp(-gamma*T)))
end

function approximate_protein_variance(t)
    c_double_prime = (approximate_F(T)- 4 * approximate_F(0) + approximate_protein_mean(T))/3
    return approximate_F(t) + c_double_prime
end

function approximate_F(t)
    2*beta^2*alpha/gamma*(t/gamma - 1/(2-exp(-gamma*T))*1/gamma*(1-t)*exp(-gamma*t) + 1/gamma*exp(-gamma*t)/(4-exp(-gamma*T))*(3/gamma + T*exp(-gamma*T)/(2-exp(-gamma*T)))) + beta*alpha/gamma * (t + 1/gamma*exp(-gamma * t)/(2-exp(-gamma*T)))
end

function simulate_network()
    write(stdout, "Starting simulation")
    # This is heuristic but seems right
    t_final = T*1000000
    time_since_division = 0
    t = 0
    m = 0
    A = alpha/gamma - m
    p = 0
    p_after_division = []
    p_before_division = []
    m_after_division = []
    m_before_division = []
    
    species = [m, p]
    propensities = [
        species -> alpha,
        species -> beta * species[1],
        species -> gamma * species[1],
        species -> delta * species[2],
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
            end

            species = partition_species(species)

            if t > t_final/1000
                append!(m_after_division, species[1])
                append!(p_after_division, species[2])
            end
        else
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
    end

    return m_before_division, m_after_division, p_before_division, p_after_division
end

function partition_species(species)
    dist1 = Binomial(species[1])
    dist2 = Binomial(species[2])
    return [rand(dist1), rand(dist2)]
end

main(ARGS)
