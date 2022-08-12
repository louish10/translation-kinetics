using Catalyst
using Symbolics
using SteadyStateDiffEq
using DiffEqBase
using LinearAlgebra
using StatsBase
using MatrixEquations

struct LNAProblem{T,NT,uType,P}
    rn::T
    u0::uType
    p::P
    J::Matrix{NT}
    D::Matrix{NT}
end

DiffEqBase.remake(prob::LNAProblem; u0=prob.u0, p=prob.p) = LNAProblem(prob.rn, u0, p, prob.J, prob.D)

function LNAProblem(rn::ReactionSystem, u0, p=DiffEqBase.NullParameters(); combinatoric_ratelaw=true)
    if length(u0) < length(species(rn))
        error("u0 length error. Expected array of length $(length(u0)). Got $(length(species(rn))).")
    end
    states = Catalyst.get_states(rn)
    eqs = Catalyst.get_eqs(rn)
    fs = [ Catalyst.oderatelaw(rx; combinatoric_ratelaw=combinatoric_ratelaw) for rx in eqs ]
    S = Catalyst.netstoichmat(rn)
    
    J = Num[ sum(S[i,r] * Symbolics.derivative(fs[r], states[j]) for r in 1:length(eqs)) for i in 1:length(states), j in 1:length(states) ]
    D = Num[ sum(S[i,r] * S[j,r] * fs[r] for r in 1:length(eqs)) for i in 1:length(states), j in 1:length(states) ]
    
    LNAProblem(rn, u0, p, J, D)
end

get_drift(prob::LNAProblem) = prob.J
get_diff(prob::LNAProblem) = prob.D

struct LNASolution{T,NT,uType,P}
    rn::T
    u::uType
    p::P
    J::Matrix{NT}
    D::Matrix{NT}
end

DiffEqBase.solve(prob::LNAProblem, p=prob.p; kwargs...) = solve(remake(prob; p=p); kwargs...)

function DiffEqBase.solve(prob::LNAProblem; solver=SSRootfind())
    ssprob = SteadyStateProblem(prob.rn, prob.u0, prob.p)
    u = solve(ssprob, solver)
    
    states = Catalyst.get_states(prob.rn)
    ps = Catalyst.get_ps(prob.rn)
    
    subs_dict_vals = Dict(states[i] => u[i] for i in 1:length(states))
    subs_dict_ps = Dict(ps[i] => prob.p[i] for i in 1:length(ps))
    subs_dict = union(subs_dict_vals, subs_dict_ps)
    
    NT_new = eltype(u)
    
    subs = x -> convert(NT_new, Catalyst.value(substitute(x, subs_dict)))
    J = subs.(prob.J)
    D = subs.(prob.D)
    
    LNASolution(prob.rn, u, prob.p, J, D)
end

StatsBase.mean(lna::LNASolution) = lna.u
StatsBase.cov(lna::LNASolution) = lyapc(lna.J, lna.D)

function StatsBase.autocov(lna::LNASolution, lags)
    σ = cov(lna)
    
    map(lag -> LinearAlgebra.exp(lna.J * abs(lag)) * σ, lags)
end

function power(lna::LNASolution, ω)
    A = lna.J + UniformScaling(ω * im)
    ret = A \ lna.D
    ret = ret'
    ret = A \ ret
    real.(1 / (2 * pi) * ret)
end