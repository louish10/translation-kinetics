using DifferentialEquations, StatsBase, Plots, Distributions, LaTeXStrings, DataFrames, CSV, ProgressBars, Optim


function reaction(n)
    "
    ρ_m$(n), 0 --> M$(n)
    δ_m$(n), M$(n) --> 0
    σ_b$(n), M$(n) + R --> MR$(n)
    δ_u$(n), MR$(n) --> M$(n) + R
    ρ_p$(n), MR$(n) --> M$(n) + P$(n) + R
    δ_mr$(n), MR$(n) --> 0
    δ_p$(n), P$(n) --> 0
    "
end

function endstring(n)
    cumulative_string = "end ρ_r δ_r "
    for i in 1:n
        cumulative_string = string(cumulative_string, "ρ_m$(i) δ_m$(i) σ_b$(i) δ_u$(i) ρ_p$(i) δ_mr$(i) δ_p$(i) ")
    end
    cumulative_string
end

function generate_rn(n)
    rnstring = "rn = @reaction_network begin
        ρ_r, 0 --> R
        δ_r, R --> 0
    "
    for i in 1:n
        rnstring = string(rnstring, reaction(i))
    end
    rnstring = string(rnstring, endstring(n))
    rnexpression = Base.Meta.parse(rnstring)
    eval(rnexpression)
end

function relative_error(true_val, approx_val)
    map(abs, 1 .- approx_val./true_val)
end

function rn_can_f()
    return @reaction_network begin
        k1, 0 --> m
        k2, m --> 0
        k3, m --> m + p
        k6, p --> 0
    end k1 k2 k3 k6
end

function rn_rib_f()
    return @reaction_network begin
        k1, 0 --> r
        k2, r --> 0
        k3, 0 --> m
        k4, m --> 0
        k5, m + r --> m + r + p
        k6, p --> 0
    end k1 k2 k3 k4 k5 k6
end

function rn_rib_bound_f()
    return @reaction_network begin
        k1, 0 --> r
        k2, r --> 0
        k3, 0 --> m
        k4, m --> 0
        k5, m + r --> m + r + p
        k6, p --> 0
        k7, r --> rb
        k8, rb --> r
        k9, rb --> 0
    end k1 k2 k3 k4 k5 k6 k7 k8 k9
end

function rn_rib_bound_int_f()
    return @reaction_network begin
        k1, 0 --> r
        k2, r --> 0
        k3, 0 --> m
        k4, m --> 0
        k5, m + r --> rs
        k6, rs --> m + r
        k7, rs --> m + r + p
        k8, p --> 0
        k9, rs --> 0 
        k10, r --> rb
        k11, rb --> r
        k12, rb --> 0
    end k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12
end

function ps_can_f(ps, lnameans_bound) # ps frmo rib_bound
    return [
        ps[3],
        ps[4],
        ps[5]*lnameans_bound[1],
        ps[2]
    ]
end

function ps_rib_bound_f(ps, lnameans_bound_int) # ps from rib bound int
    
    k2 = ps[4] + (ps[5]*ps[2]/(ps[6]+ps[7]+ps[2]))*lnameans_bound_int[1]
    #k2 = ps[3]/(lnameans_bound_int[2] + lnameans_bound_int[3])
    kb = (ps[2]+ps[11]) * (ps[5] / (ps[6] + ps[7] + ps[2])*lnameans_bound_int[2] + ps[10]/(ps[11] + ps[2]))
    return [
        ps[1], # rib creation
        ps[2], # rib degredation
        ps[3], # mrna creation
        k2, # mrna degredation
        ps[5]*ps[7]/(ps[6]+ps[7]+ps[2]), # translation
        ps[2], # protein degredation = rib degredation
        kb, # rib binding
        ps[11], # rib unbinding
        ps[2], # bound ribosome degredation
    ]
end

function ps_rib_bound_int_f(ps, lnameans) # ps from whole cell
    k3_values = ps[5+7:7:length(ps)]
    k4_values = ps[6+7:7:length(ps)]
    k5_values= ps[7+7:7:length(ps)]
    k6 = ps[2]
    m_values = lnameans[2+3:3:length(lnameans)]
    sigma = sum((k3_values .- ((k4_values .+ k5_values) .* k3_values) ./ (k4_values .+ k5_values .+ k6)) .* m_values)
    
    rib_unbinding = 1
    rib_binding = sigma *(rib_unbinding + k6)/k6 
    return [
        ps[1], # rib creation
        ps[2], # rib degredation
        ps[3], # mrna creation
        ps[4], # mrna degredation
        ps[5],# rib binding
        ps[6],# rib unbinding
        ps[7], # translation
        ps[2], # protein degredation = rib degredation
        ps[2], # monosome degredation
        rib_binding, # rib binding
        rib_unbinding, # rib unbinding
        ps[2], # bound ribosome degredation
    ]
end

function negative_binomial_ys(xs, mean, var)
    p = mean/var
    r = var*p^2/(1-p)

    nb = NegativeBinomial(r,p)
    return Distributions.pdf(nb, xs)
end

function poisson_ys(xs, mean)
    pois = Poisson(mean)
    return Distributions.pdf(pois, xs)
end

function log_likelihood_can(m_data::Vector, p_data::Vector, k1::Float64, k2::Float64, k3::Float64, k6::Float64)
    rn_can = rn_can_f()
    prob = LNAProblem(rn_can, zeros(length(species(rn_can))), [k1, k2, k3, k6])
    sol = solve(prob)
    lnameans = mean(sol)
    lnacovs = StatsBase.cov(sol)
    
    norm = MvNormal(lnameans, lnacovs)
    
    log_likelihood = 0
    for i in 1:length(m_data)
        log_likelihood += loglikelihood(norm, [m_data[i], p_data[i]])
    end
    - log_likelihood
end

function log_likelihood_can(p_data::Vector, k1::Float64, k2::Float64, k3::Float64, k6::Float64)
    rn_can = rn_can_f()
    prob = LNAProblem(rn_can, zeros(length(species(rn_can))), [k1, k2, k3, k6])
    sol = solve(prob)
    lnameans = mean(sol)
    lnacovs = StatsBase.cov(sol)
    
    norm = Normal(lnameans[2], sqrt(lnacovs[2,2]))
    
    log_likelihood = loglikelihood(norm, p_data)
    - log_likelihood
end

function log_likelihood_rib_bound(r_data, m_data, p_data, rib_bound_data, k6, k7, k1, k2, k3, kb, kmb)  
    rn_bound = rn_rib_bound_f()
    prob = LNAProblem(rn_bound, zeros(length(species(rn_bound))), [k7, k6, k1, k2, k3, k6, kb, kmb, k6])
    sol = solve(prob)
    lnameans = mean(sol)
    lnacovs = StatsBase.cov(sol)
    
    norm = MvNormal(lnameans, lnacovs)
    
    log_likelihood = 0
    for i in 1:length(m_data)
        log_likelihood += loglikelihood(norm, [r_data[i], m_data[i], p_data[i], rib_bound_data[i]])
    end
    - log_likelihood
end

function log_likelihood_rib_bound(
        m_data::Vector,
        p_data::Vector,
        p1::Float64, p2::Float64, p3::Float64,
        p4::Float64,
        p5::Float64,
        p6::Float64,
        p7::Float64,
        p8::Float64,
        p9::Float64)
    rn_bound = rn_rib_bound_f()
    prob = LNAProblem(rn_bound, zeros(length(species(rn_bound))), [p1,p2,p3,p4,p5,p6,p7,p8,p9])
    sol = solve(prob)
    lnameans = mean(sol)
    lnacovs = StatsBase.cov(sol)
    
    lnameans = [lnameans[2], lnameans[3]]
    
    lnacovs = [lnacovs[2,2] lnacovs[2,3];
               lnacovs[3,2] lnacovs[3,3]]
    
    norm = MvNormal(lnameans, lnacovs)
    
    log_likelihood = 0
    for i in 1:length(m_data)
        log_likelihood += loglikelihood(norm, [m_data[i], p_data[i]])
    end
    - log_likelihood
end

function log_likelihood_rib_bound(
        p_data::Vector,
        p1::Float64, p2::Float64, p3::Float64,
        p4::Float64,
        p5::Float64,
        p6::Float64,
        p7::Float64,
        p8::Float64,
        p9::Float64)
    rn_bound = rn_rib_bound_f()
    prob = LNAProblem(rn_bound, zeros(length(species(rn_bound))), [p1,p2,p3,p4,p5,p6,p7,p8,p9])
    sol = solve(prob)
    lnameans = mean(sol)
    lnacovs = StatsBase.cov(sol)
    
    norm = Normal(lnameans[3], sqrt(lnacovs[3,3]))
    
    - loglikelihood(norm, p_data)
end

function log_likelihood_rib_bound_int(r_data, m_data, rib_int_data, p_data, rib_bound_data, k1, k2, k6, k7, k3, k4, k5, kb, kmb)  
    rn_bound_int = rn_rib_bound_int_f()
    prob = LNAProblem(rn_bound_int, zeros(length(species(rn_bound_int))), [k7, k6, k1, k2, k3, k4, k5, k6, k6, kb, kmb, k6])
    sol = solve(prob)
    lnameans = mean(sol)
    lnacovs = StatsBase.cov(sol)
    
    norm = MvNormal(lnameans, lnacovs)
    
    log_likelihood = 0
    for i in 1:length(m_data)
        log_likelihood += loglikelihood(norm, [r_data[i], m_data[i], rib_int_data[i], p_data[i], rib_bound_data[i]])
    end
    - log_likelihood
end

function log_likelihood_rib_bound_int(
        p_data::Vector, 
        k1::Float64,
        k2::Float64,
        k3::Float64,
        k4::Float64,
        k5::Float64,
        k6::Float64,
        k7::Float64,
        k8::Float64,
        k9::Float64,
        k10::Float64,
        k11::Float64,
        k12::Float64
)
    rn_bound_int = rn_rib_bound_int_f()
    prob = LNAProblem(rn_bound_int, zeros(length(species(rn_bound_int))), [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12])
    sol = solve(prob)
    lnameans = mean(sol)
    lnacovs = StatsBase.cov(sol)
    
    norm = Normal(lnameans[4], sqrt(lnacovs[4,4]))
    
    - loglikelihood(norm, p_data)
end

function log_likelihood_rib_bound_int(
        m_data::Vector,
        p_data::Vector,
        ps::Vector
    )  
    rn_bound_int = rn_rib_bound_int_f()
    prob = LNAProblem(rn_bound_int, zeros(length(species(rn_bound_int))), ps)
    sol = solve(prob)
    lnameans = mean(sol)
    lnacovs = StatsBase.cov(sol)
    
    lnameans = [lnameans[2], lnameans[4]]
    
    lnacovs = [lnacovs[2,2] lnacovs[2,4];
               lnacovs[4,2] lnacovs[4,4]]
    
    
    norm = MvNormal(lnameans, lnacovs)
    
    log_likelihood = 0
    for i in 1:length(m_data)
        log_likelihood += loglikelihood(norm, [m_data[i], p_data[i]])
    end
    - log_likelihood
end

function log_likelihood_whole_cell(data, ps)
    rn = generate_rn(10)
    prob = LNAProblem(rn, zeros(length(species(rn))), ps)
    sol = solve(prob)
    lnameans = mean(sol)
    lnacovs = StatsBase.cov(sol)
    
    norm = MvNormal(lnameans, lnacovs)
    
    log_likelihood = 0
    for i in 1:length(m_data)
        log_likelihood += loglikelihood(norm, data[i])
    end
    - log_likelihood
    
    return likelihood
end

function log_likelihood_whole_cell(p_data::Vector, ps::Vector)
    rn = generate_rn(10)
    prob = LNAProblem(rn, zeros(length(species(rn))), ps)
    sol = solve(prob)
    lnameans = mean(sol)
    lnacovs = StatsBase.cov(sol)
    
    norm = Normal(lnameans[4], sqrt(lnacovs[4,4]))
    
    
    -loglikelihood(norm, p_data)
end