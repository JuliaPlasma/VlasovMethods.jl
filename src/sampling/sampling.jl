

abstract type SamplingMethod end

struct NoSampling <: SamplingMethod end
struct ImportanceSampling <: SamplingMethod end
struct AcceptRejectSampling <: SamplingMethod end


"""
Returns re-weighted particles according to new parameters.

input: Particles struct P₀, proposal distribution function g, target distribution function f
output: Particles struct
"""
# function weight_f(P₀, g, f)
#     Px = zero(P₀.x); Pv = zero(P₀.v); Pw = zero(P₀.w)
#     Px .= P₀.x; Pv .= P₀.v; Pw .= P₀.w
#     P = ParticleList(Px, Pv, Pw)
#     P.w .*= f(P₀.x, P₀.v) ./ g(P₀.x, P₀.v)
#     return P
# end
