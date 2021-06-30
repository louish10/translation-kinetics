using Plots, Statistics, Distributions, SpecialFunctions

function main(args)
    global path = mkpath("data/simple-birth-death-process")

    global alpha = 5.0
    global gamma = 0.1
    global T = 10.0

    time_averaged_measurements = simulate_time_average()
    random_time_measurements = simulate_random_measurements()

    generate_summary(time_averaged_measurements, random_time_measurements)
end

function generate_summary(time_averaged_measurements, random_time_measurements)
    filename = string(path, "/results.txt")

    time_averaged_mean = mean(time_averaged_measurements)
    random_time_mean = mean(random_time_measurements)
    open(filename, "w") do io
        write(io, "Time average: $time_averaged_mean\n")
        write(io, "Random sample average: $random_time_mean\n")
    end
end

function simulate_time_average()
    write(stdout, "Starting simulation")
    t_final = T*100000
    time_since_division = 0
    t = 0

    time_averages = []

    species = [0]
    propensities = [
        species -> alpha,
        species -> gamma * species[1],
    ]

    write(stdout, "\n")

    copy_numbers = []
    time_intervals = []
    last_reaction_time = 0

    steady_state = false

    while true
        if t%1000 == 0
            write(stdout, "\r$t of $t_final")
        end

        instantaneous_propensities = map(f -> f(species), propensities)
        a = sum(instantaneous_propensities)

        dt = 1/a*log(1/(1-rand()))
        t = t + dt
        time_since_division += dt

        if t>t_final/10
            steady_state = true
        end

        if t > t_final
            break
        end

        if time_since_division > T
            t = t - t%T
            time_since_division = 0

            if steady_state
                append!(time_intervals, t-last_reaction_time)
                last_reaction_time = t

                append!(copy_numbers, species[1])
                append!(time_averages, 1/T*sum(copy_numbers .* time_intervals))
            end

            copy_numbers = []
            time_intervals = []

            species = partition_species(species)
        else
            R = rand()
            if R < instantaneous_propensities[1]/a
                append!(time_intervals, t-last_reaction_time)
                append!(copy_numbers, species[1])
                last_reaction_time = t

                species[1] += 1
            else
                append!(time_intervals, t-last_reaction_time)
                append!(copy_numbers, species[1])
                last_reaction_time = t

                species[1] -= 1
            end
        end
    end
    write(stdout, "\n")

    return time_averages
end

function simulate_random_measurements()
    write(stdout, "Starting simulation")
    t_final = T*100
    N = 10_000
    measurements = []

    write(stdout, "\n")

    for i in 1:N
        write(stdout, "\r$i of $N")
        time_since_division = 0
        t = 0
        species = [0]
        propensities = [
            species -> alpha,
            species -> gamma * species[1],
        ]
        while true
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

                species = partition_species(species)
            else
                R = rand()
                if R < instantaneous_propensities[1]/a
                    species[1] += 1
                else
                    species[1] -= 1
                end
            end
        end
        append!(measurements, species[1])
    end
    write(stdout, "\n")

    return measurements
end

function partition_species(species)
    dist1 = Binomial(species[1])
    return [rand(dist1)]
end

main(ARGS)
