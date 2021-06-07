using Plots, Statistics, Distributions, SpecialFunctions
include("model_with_ribosomes_functions.jl")

function main(args)
    global alpha = parse(Float64, args[1])
    global beta = parse(Float64, args[2])
    global gamma = parse(Float64, args[3])
    global rho = parse(Float64, args[4])
    global T = parse(Float64, args[5])

    global path = mkpath("data/model-with-ribosomes-time-average/alph-$alpha-bet-$beta-gam-$gamma-rho-$rho-T-$T")

    p_time_averages = simulate_network()

    generate_summary(p_time_averages)
end

function generate_summary(p_time_averages)
    filename = string(path, "/results.txt")
    touch(filename)

    stat_mean = mean(p_time_averages)
    theory_mean = ModelWithRibosomes.p_time_av(alpha, beta, gamma, rho, T)
    
    open(filename, "w") do io
        write(io, "Protein\n\n")
        write(io, "Mean\n")
        write(io, "Theory: $theory_mean\n")
        write(io, "Simulation: $stat_mean\n")
    end
end

function simulate_network()
    write(stdout, "Starting simulation")
    # This is heuristic but seems right
    t_final = T*100000
    time_since_division = 0
    t = 0

    p_time_averages = []

    species = [0, 0, 0]
    propensities = [
        species -> alpha,
        species -> beta * species[1] * species[3],
        species -> gamma * species[1],
        species -> rho
    ]

    write(stdout, "\n")

    protein_numbers = []
    time_intervals = []
    last_protein_time = 0

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

            append!(time_intervals, t-last_protein_time)
            append!(protein_numbers, species[2])
            append!(p_time_averages, 1/T*sum(protein_numbers .* time_intervals))

            protein_numbers = []
            time_intervals = []

            species = partition_species(species)
        else
            R = rand()
            if R < instantaneous_propensities[1]/a
                species[1] += 1
            elseif R < sum(instantaneous_propensities[1:2])/ a
                append!(time_intervals, t-last_protein_time)
                append!(protein_numbers, species[2])
                last_protein_time = t

                species[2] += 1
            elseif R < sum(instantaneous_propensities[1:3])/ a
                species[1] -= 1
            else
                species[3] += 1
            end
        end
    end

    return p_time_averages
end

function partition_species(species)
    dist1 = Binomial(species[1])
    dist2 = Binomial(species[2])
    dist3 = Binomial(species[3])
    return [rand(dist1), rand(dist2), rand(dist3)]
end

main(ARGS)
