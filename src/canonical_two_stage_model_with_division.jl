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

    m_zero = alpha/gamma*(1-exp(-gamma*T))/(2-exp(-gamma*T))
    m_T = alpha/gamma*(1-exp(-gamma*T)) + m_zero*exp(-gamma*T)
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
    int_constant = 1/(4-exp(-2*gamma*T))*(4*A-2*A*exp(-gamma*T) - 2*alpha/gamma)
    theory_var = alpha/gamma - A + int_constant

    write(stdout, "\nStatistical mean of mrna after division: $mrna_stat_mean")
    write(stdout, "\nTheoretical mean of mrna after division: $m_zero")
    write(stdout, "\nStatistical variance of mrna after division: $stat_var")
    write(stdout, "\nTheoretical variance of mrna after division: $theory_var")

    statistical_covariance = cov(m_after, p_after)
    cov_int_const = 1/(4-exp(-(gamma+delta)*T))*(-3*alpha/(gamma*(gamma+delta)) + (4-exp(-gamma*T))*(A/delta-int_constant/(delta-gamma)))
    theoretical_covariance = beta*(alpha/(gamma*(gamma+delta)) - A/delta+int_constant/(delta-gamma) + cov_int_const)
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
    protein_int_const = 1/(4-exp(-2*delta*T))*(F(int_constant, cov_int_const, T) - 4*F(int_constant, cov_int_const, 0) + p_T)
    protein_theory_var = F(int_constant, cov_int_const, 0) + protein_int_const
    write(stdout, "\nStatistical variance of protein after division: $protein_stat_var")
    write(stdout, "\nTheoretical variance of protein after division: $protein_theory_var")

    print_summary(p_before, p_after, m_before, m_after)
end

function print_summary(p_before, p_after, m_before, m_after)
    mrna_after_stat_mean = mean(m_after)
    m_zero = alpha/gamma*(1-exp(-gamma*T))/(2-exp(-gamma*T))
    A = alpha/gamma - m_zero
    mrna_after_stat_var = var(m_after)
    c = 1/(4-exp(-2*gamma*T))*(4*A-2*A*exp(-gamma*T) - 2*alpha/gamma)
    mrna_after_theory_var = alpha/gamma - A + c

    covar_after_stat = cov(m_after, p_after)
    c_prime = 1/(4-exp(-(gamma+delta)*T))*(-3*alpha/(gamma*(gamma+delta)) + (4-exp(-gamma*T))*(A/delta-c/(delta-gamma)))
    covar_after_theory = beta*(alpha/(gamma*(gamma+delta)) - A/delta+c/(delta-gamma) + c_prime)

    p_after_stat_mean = mean(p_after)
    p_zero = 1/(2-exp(-delta*T))*(beta*alpha/(delta*gamma)-A*beta/(delta-gamma)*exp(-gamma*T)-(beta*alpha/(gamma*delta)+A*beta/(gamma-delta))*exp(-delta*T))
    p_T = beta*alpha/(delta*gamma) - A*beta/(delta-gamma)*exp(-gamma*T) + (p_zero - beta*alpha/(gamma*delta)-A*beta/(gamma-delta))*exp(-delta*T)
    p_after_stat_var = var(p_after)
    c_double_prime = 1/(4-exp(-2*delta*T))*(F(c, c_prime, T) - 4*F(c, c_prime, 0) + p_T)
    p_after_theory_var = F(c, c_prime, 0) + c_double_prime

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
    end
end

function F(c, c_prime, t)
    m_zero = alpha/gamma*(1-exp(-gamma*T))/(2-exp(-gamma*T))
    A = alpha/gamma - m_zero
    p_zero = 1/(2-exp(-delta*T))*(beta*alpha/(delta*gamma)-A*beta/(delta-gamma)*exp(-gamma*T)-(beta*alpha/(gamma*delta)+A*beta/(gamma-delta))*exp(-delta*T))
    B = p_zero - beta*alpha/(gamma*delta) - A*beta/(gamma-delta)
    
    return beta^2*(alpha/(gamma*delta*(gamma+delta))-2*A/(delta*(2*delta - gamma))*exp(-gamma*t) + c/((delta-gamma)^2)*exp(-2*gamma*t) + 2*c_prime/(delta-gamma)*exp(-(gamma+delta)*t)) + beta*alpha/(gamma*delta) - A*beta/(2*delta-gamma)*exp(-gamma*t) + beta*delta*A/((gamma-delta)*(2*delta-gamma))*exp(-gamma*t) + B*exp(-delta*t)
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
