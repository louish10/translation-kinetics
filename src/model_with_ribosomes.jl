using Plots, Statistics, Distributions
include("utilities.jl")

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

    nbins = maximum([1, Int(maximum(m_before))])
    histogram(m_before, nbins=Int(nbins), normed=true)
    vline!([mrna(T)], label="theoretical mean")
    vline!([mean(m_before)], label="statistical mean")
    xlabel!("m")
    ylabel!("P(m)")
    savefig(string(path, "/mrnas_before.svg"))

    nbins = maximum([1, Int(maximum(m_after))])
    histogram(m_after, nbins=Int(nbins), normed=true)
    vline!([mrna(0)], label="theoretical mean")
    vline!([mean(m_after)], label="statistical mean")
    xlabel!("m")
    ylabel!("P(m)")
    savefig(string(path, "/mrnas_after.svg"))

    nbins = maximum([1, Int(maximum(p_after))])
    histogram(p_after, nbins=Int(nbins), normed=true)
    vline!([p(0)], label="theoretical mean")
    vline!([mean(p_after)], label="statistical mean")
    xlabel!("p")
    ylabel!("P(p)")
    savefig(string(path, "/proteins_after.svg"))

    nbins = maximum([1, Int(maximum(p_before))])
    histogram(p_before, nbins=Int(nbins), normed=true)
    vline!([p(T)], label="theoretical mean")
    vline!([mean(p_before)], label="statistical mean")
    xlabel!("p")
    ylabel!("P(p)")
    savefig(string(path, "/proteins_before.svg"))

    nbins = maximum([1, Int(maximum(r_before))])
    histogram(r_before, nbins=Int(nbins), normed=true)
    vline!([rho*T*2], label="theoretical mean")
    vline!([mean(r_before)], label="statistical mean")
    xlabel!("r")
    ylabel!("P(r)")
    savefig(string(path, "/ribosomes_before.svg"))

    nbins = maximum([1, Int(maximum(r_after))])
    histogram(r_after, nbins=Int(nbins), normed=true)
    vline!([rho*T], label="theoretical mean")
    vline!([mean(r_after)], label="statistical mean")
    xlabel!("r")
    ylabel!("P(r)")
    savefig(string(path, "/ribosomes_after.svg"))

    statistical_variance = var(m_after)
    theoretical_variance = mrna_var(0)

    generate_summary(m_before, m_after, p_before, p_after, r_before, r_after)
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

function simulate_network()
    write(stdout, "Starting simulation")
    # This is heuristic but seems right
    t_final = T*1000000
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
    alpha*beta*rho/(2*gamma^3*(2*exp(gamma*T)-1)) * (2 + 4*T*gamma - (t^2 + 2*t*T + 3*T^2)*gamma^2 + 2*exp((T-t)*gamma)*(1+(t+T)*gamma) + 2*exp(T*gamma)*(-2 + gamma*(t^2*gamma + T*(-2 + 2*t*gamma + 3*T*gamma))))
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
    return exp(-2*gamma*t)*alpha*beta*rho/(6*(1-2*exp(gamma*T))^2*(4*exp(T*gamma)-1)*gamma^5)*(-3*exp(2*T*gamma)*alpha*beta*(2*(t+T)*gamma+3)+12*exp(3*T*gamma)*alpha*beta*(2*(t+T)*gamma+3)+4*exp(2*t*gamma+3*T*gamma)*(4*(gamma*(3*gamma*t^2+T*(6*(t+T)*gamma-5))-5)*gamma^2+alpha*beta*(2*gamma*(4*t*(t^2+3*T*t+3*T^2)*gamma^2-2*T*(8*t+3*T)*gamma-16*t+5*T)+31)+4*beta*(gamma*(-6*gamma*t^2+4*(t^2+3*T*t+3*T^2)*gamma^2*t-15*T-12*T*(t+T)*gamma)-9)*rho)-exp(2*(t+T)*gamma)*(4*(gamma*(15*gamma*t^2+T*(30*(t+T)*gamma-31))-23)*gamma^2+alpha*beta*(2*gamma*(20*t*(t^2+3*T*t+3*T^2)*gamma^2-2*T*(40*t+9*T)*gamma-64*t+21*T)+95)+4*beta*(gamma*(9*gamma^2*T^3+3*gamma*(20*t*gamma-11)*T^2+6*(10*t*gamma*(t*gamma-1)-7)*T+10*t^2*gamma*(2*t*gamma-3))-18)*rho)-4*exp((t+2T)*gamma)*(alpha*beta*(9*(t+T)^2*gamma^2-12*T*gamma-26)+18*beta*rho+3*gamma*((t+T)*gamma+1)*(3*gamma+beta*(6*t+14*T+3*(t^2+2*T*t-T^2)*gamma)*rho))+16*exp((t+3*T)*gamma)*(alpha*beta*(3*(t+T)^2*gamma^2-4*T*gamma-10)+9*beta*rho+3*gamma*((t+T)*gamma^2+gamma+beta*(t*(t+T)*(t+2*T)*gamma^2+(t+2*T)*(3*t+2*T)*gamma+3*t+5*T)*rho))+2*exp((2*t+T)*gamma)*(alpha*beta*(8*t*(t^2+3*T*t+3*T^2)*gamma^3-T*(32*t+3*T)*gamma^2+4*(T-5*t)*gamma+8)+gamma*(gamma*(gamma*(12*gamma*t^2+T*(24*(t+T)*gamma-29))-17)+beta*(9*gamma^2*T^3+3*gamma*(16*t*gamma-3)*T^2+12*(4*t*gamma*(t*gamma-1)-1)*T+8*t^2*gamma*(2*t*gamma-3))*rho))+2*exp((t+T)*gamma)*(alpha*beta*(3*(t+T)^2*gamma^2-4*T*gamma-8)+3*gamma*((t+T)*((t-T)*(t+3*T)*beta*rho+1)*gamma^2+(3*t+T)*(t+3*T)*beta*rho*gamma+gamma+4*T*beta*rho))-exp(2*t*gamma)*gamma*(2*beta*gamma^2*(alpha+2*rho)*t^3+3*gamma*(gamma^2+2*T*beta*(alpha+2*rho)*gamma-2*beta*rho)*t^2+2*(alpha*beta*(T*gamma*(3*T*gamma-4)-2)+3*T*gamma*(gamma^2+2*beta*(T*gamma-1)*rho))*t+2*gamma*(T*gamma*(3*T*gamma-4)-2)))

end

function m_r_covar(t)
    return 0
end

function p_r_covar(t)
   alpha*beta*rho/(6*gamma^3*(2*exp(gamma*T)-1))*(-3*gamma^2*(t+T)^2+2*exp(gamma*T)*(3*gamma^2*(t+T)^2 - 4*gamma*T - 4) + 6*exp(gamma*(T-t))*(gamma*(t+T)+1)+4*gamma*T+2) 
end

function m_p_covar(t)
    exp(-t*gamma)*alpha*beta*rho*(2*exp(t*gamma)*(gamma*(t+T)-1)-12*exp(gamma*(t+T))*(-1 +(t+T)*gamma)+16*exp((t+2*T)*gamma)*(-1+(t+T)*gamma)-4*exp(2*T*gamma)*(-3 + 2*T*gamma + t*(t+2T)*gamma^2)+exp(T*gamma)*(-6 +4*T*gamma + (t-T)*(t+3*T)*gamma^2))/(2*(1-6*exp(T*gamma)+8*exp(2*T*gamma))*gamma^3)
end

function partition_species(species)
    dist1 = Binomial(species[1])
    dist2 = Binomial(species[2])
    dist3 = Binomial(species[3])
    return [rand(dist1), rand(dist2), rand(dist3)]
end

main(ARGS)
