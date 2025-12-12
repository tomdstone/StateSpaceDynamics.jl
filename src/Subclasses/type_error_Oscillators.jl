
# # ! somehow this is erroring, with a type error wrt. T

# mutable struct OscillatorStateModel{
#     T<:Real,M<:AbstractMatrix{T},V<:AbstractVector{T}
# } <: AbstractStateModel{T}
#     a::V
#     ω::V
#     f::Union{V,Nothing}
#     fs::Union{T,Nothing}
#     σ::V
#     A::M
#     Q::M
#     b::V
#     x0::V
#     P0::M
#     # Q_prior::Union{Nothing,IWPrior{T}} = nothing
#     # P0_prior::Union{Nothing,IWPrior{T}} = nothing

# end

# function OscillatorStateModel(;
#     a::V,
#     σ::V,
#     P0::M,
#     fs::Union{T,Nothing}=nothing,
#     f::Union{V,Nothing}=nothing,
#     ω::Union{V,Nothing}=nothing,
#     x0::Union{V,Nothing}=nothing,
# ) where {T<:Real,M<:AbstractMatrix{T},V<:AbstractVector{T}}
#     fs, f, ω = _unify_freq_components(fs, f, ω)
#     A = _build_oscillator_A(a, ω)
#     Q = diagm(repeat(σ, inner=2))
#     b = zeros(T, 2length(a))
#     x0 = zeros(T, 2length(a))

#     return new{T,M,V}(a, ω, f, fs, σ, A, Q, x0, b, P0)
# end
