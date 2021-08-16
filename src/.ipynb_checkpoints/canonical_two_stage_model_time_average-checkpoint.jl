using Plots, Statistics, Distributions, SpecialFunctions
include("canonical_two_stage_model_with_division_functions.jl")

function main(args)
    global alpha = parse(Float64, args[1])
    global beta = parse(Float64, args[2])
    global gamma = parse(Float64, args[3])
    global T = parse(Float64, args[4])

    global path = mkpath("data/canonical-two-stage-model-time-average/alph-$alpha-bet-$beta-gam-$gamma-T-$T")

    p_time_averages, p_squared_time_averages, mrna_time_averages, mrna_squared_time_averages = simulate_network()

    generate_summary(p_time_averages, p_squared_time_averages, mrna_time_averages, mrna_squared_time_averages )
end


function generate_summary(p_time_averages, p_squared_time_averages, mrna_time_averages, mrna_squared_time_averages)
    filename = string(path, "/results.txt")
    touch(filename)

    stat_mean = mean(p_time_averages)
    theory_mean = CanonicalTwoStageModel.p_time_av(alpha, beta, gamma, T)

    stat_var = mean(p_squared_time_averages) - mean(p_time_averages).^2
    theory_var = CanonicalTwoStageModel.p_var_time_av(alpha, beta, gamma, T)
    
    mrna_stat_mean = mean(mrna_time_averages)
    mrna_theory_mean = CanonicalTwoStageModel.m_time_av(alpha, gamma, T)

    mrna_stat_var = mean(mrna_squared_time_averages) - mean(mrna_time_averages).^2
    mrna_theory_var = CanonicalTwoStageModel.m_var_time_av(alpha, gamma, T)

    open(filename, "w") do io
        write(io, "Protein\n\n")
        write(io, "Mean\n")
        write(io, "Theory: $theory_mean\n")
        write(io, "Simulation: $stat_mean\n")
        write(io, "Variance\n")
        write(io, "Theory: $theory_var\n")
        write(io, "Simulation: $stat_var\n")

        write(io, "\n\nMrna\n\n")
        write(io, "Mean\n")
        write(io, "Theory: $mrna_theory_mean\n")
        write(io, "Simulation: $mrna_stat_mean\n")
        write(io, "Variance\n")
        write(io, "Theory: $mrna_theory_var\n")
        write(io, "Simulation: $mrna_stat_var\n")
    end
end

function simulate_network()
    write(stdout, "Starting simulation")
    # This is heuristic but seems right
    t_final = T*1000
    time_since_division = 0
    t = 0

    p_time_averages = []
    p_squared_time_averages = []
    mrna_time_averages = []
    mrna_squared_time_averages = []

    species = [0, 0, 0]
    propensities = [
        species -> alpha,
        species -> beta * species[1],
        species -> gamma * species[1],
    ]

    write(stdout, "\n")

    mrna_numbers = []
    mrna_time_intervals = []
    protein_numbers = []
    time_intervals = []
    last_protein_time = 0
    last_mrna_time = 0

    steady_state = false

    while true

        write(stdout, "\r$t of $t_final")

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
                append!(time_intervals, t-last_protein_time)
                last_protein_time = t

                append!(protein_numbers, species[2])
                append!(p_time_averages, 1/T*sum(protein_numbers .* time_intervals))
                append!(p_squared_time_averages, 1/T*sum(protein_numbers .^2 .*time_intervals))

                append!(mrna_time_intervals, t-last_mrna_time)
                last_mrna_time = t

                append!(mrna_numbers, species[1])
                append!(mrna_time_averages, 1/T*sum(mrna_numbers .* mrna_time_intervals))
                append!(mrna_squared_time_averages, 1/T*sum((mrna_numbers.^2) .* mrna_time_intervals))
            end

            protein_numbers = []
            mrna_numbers = []
            time_intervals = []
            mrna_time_intervals = []

            species = partition_species(species)
        else
            R = rand()
            if R < instantaneous_propensities[1]/a
                append!(mrna_time_intervals, t-last_mrna_time)
                append!(mrna_numbers, species[1])
                last_mrna_time = t

                species[1] += 1
            elseif R < sum(instantaneous_propensities[1:2])/ a
                append!(time_intervals, t-last_protein_time)
                append!(protein_numbers, species[2])
                last_protein_time = t

                species[2] += 1
            else
                append!(mrna_time_intervals, t-last_mrna_time)
                append!(mrna_numbers, species[1])
                last_mrna_time = t

                species[1] -= 1
            end
        end
    end
    write(stdout, "\n")

    return p_time_averages, p_squared_time_averages, mrna_time_averages, mrna_squared_time_averages
end

function partition_species(species)
    dist1 = Binomial(species[1])
    dist2 = Binomial(species[2])
    return [rand(dist1), rand(dist2)]
end

function negative_binomial(x, mean, std)
    p = 1-mean/std^2
    r = mean^2/(std^2-mean)
    return SpecialFunctions.gamma.(x .+ r) ./ (SpecialFunctions.gamma.(x .+ 1) .* SpecialFunctions.gamma.(r)) .* (1-p) .^ r .* p .^x
end

function gaussian(x, mean, std)
    1/(std*sqrt(2*pi))*exp.(-1/2*(x.-mean).^2 ./(std^2))
end


main(ARGS)
