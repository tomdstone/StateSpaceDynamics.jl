"""
    OscillatorStateModel{T<:Real, M<:AbstractMatrix{T}, V<:AbstractVector{T}}

Implements a specific subclass of Gausssian State Space Model described in Matsuda and
Komaki: https://doi.org/10.1162/NECO_a_00916

Has all the fields present in a GaussianStateModel for ease of use, but makes simplifying
assumptions about certain components. First, the state transition matrix is block diagonal,
with the blocks being 2x2 damped rotation matrices. Second, the process noise covariance
is diagonal, with the noise corresponding to each 2x2 rotation matrix being a multiple
of the identity. Thus the number of latent states is twice the number of oscillators. i.e
`latent_dim = 2 * noscs`

# Oscillator State Model fields

- `a::V`: Damping factors (length `noscs`)
- `f::Union{V,Nothing}`: Rotation frequency of each oscillator (length `noscs`)
- `fs::Union{T,Nothing}`: Sampling frequency
- `ω::V`: Normalized rotation frequencies (length `noscs`)
- `σ::V`: Noise variance of each oscillator (length `noscs`)

# Gaussial State Model Fields

- `A::M`: Transition matrix (size `latent_dim × latent_dim`).
- `Q::M`: Process noise covariance matrix.
- `b::V`: Bias vector (length `latent_dim`).
- `x0::V`: Initial state mean (length `latent_dim`).
- `P0::M`: Initial state covariance (size `latent_dim × latent_dim`).
"""
mutable struct OscillatorStateModel <: AbstractStateModel{Float64}
    a
    ω
    f
    fs
    σ
    A
    Q
    b
    x0
    P0
    function OscillatorStateModel(;
        a,
        σ,
        P0,
        fs=nothing,
        f=nothing,
        ω=nothing,
        x0=nothing,
    )
        fs, f, ω = _unify_freq_components(fs, f, ω)
        A = _build_oscillator_A(a, ω)
        Q = diagm(repeat(σ, inner=2))
        b = zeros(2length(a))
        if isnothing(x0)
            x0 = zeros(2length(a))
        end

        return new(a, ω, f, fs, σ, A, Q, x0, b, P0)
    end
end

# utilities for building oscillator models

function _unify_freq_components(fs_, f_, ω_)

    if isnothing(ω_)
        if isnothing(fs_)
            error(ArgumentError("Must provide either sampling frequency or rotation angles"))
        elseif isnothing(f_)
            error(ArgumentError("Must provide either rotation frequencies or rotation angles"))
        end

        fs = fs_
        f = f_
        ω = f / fs
    elseif isnothing(fs_)
        fs = nothing
        f = nothing
        ω = ω_
    elseif isnothing(f_)
        fs = fs_
        f = ω_ * fs_
        ω = ω_
    else
        @assert f_ ≈ ω_ * fs_ "Rotation angle does not match oscillator frequency"
        fs = fs_
        f = f_
        ω = ω_
    end

    return fs, f, ω
end

function _build_oscillator_block(a0, ω0)
    c = a0 * cospi(2ω0)
    s = a0 * sinpi(2ω0)

    return [c -s; s c]
end

function _build_oscillator_A(a, ω)
    noscs = length(a)

    rows = [
        reduce(
            hcat,
            [i == j ? _build_oscillator_block(a[j], ω[j]) : zeros(2,2) for j in 1:noscs]
        ) for i in 1:noscs
    ]
    A = reduce(vcat, rows)

    return A
end

function Base.show(io::IO, osc::OscillatorStateModel; gap="")
    println(io, gap, "Oscillator State Model:")
    println(io, gap, "-----------------------")

    if length(osc.a) > 4
        println(io, gap, " State Parameters:")
        println(io, gap, "  size(a) = ($(length(osc.a)), )")
        if isnothing(osc.fs)
            println(io, gap, "  size(ω) = ($(length(osc.ω)), )")
        else
            println(io, gap, "  size(f) = ($(length(osc.f)), )")
            println(io, gap, "  fs      = $(round(osc.fs, digits = 3))")
        end
        println(io, gap, "  size(σ) = ($(length(osc.σ)), )")
        println(io, gap, " Initial State:")
        println(io, gap, "  size(x0) = ($(length(osc.x0)), )")
        println(io, gap, "  size(P0) = ($(size(osc.P0,1)), $(size(osc.P0,2)))")
    else
        println(io, gap, " State Parameters:")
        println(io, gap, "  a  = $(round.(osc.a, digits = 4))")
        if isnothing(osc.fs)
            println(io, gap, "  ω  = $(round.(osc.ω, sigdigits = 3))")
        else
            println(io, gap, "  f  = $(round.(osc.f, sigdigits = 3))")
            println(io, gap, "  fs = $(round(osc.fs, digits = 3))")
        end
        println(io, gap, "  σ  = $(round.(osc.σ, sigdigits=3))")
        println(io, gap, " Initial State:")
        println(io, gap, "  x0 = $(round.(osc.x0, digits=2))")
        println(io, gap, "  P0 = $(round.(osc.P0, sigdigits=3))")
    end

    return nothing
end

# converting to a generic Gaussian State Model
function GaussianStateModel(osc::OscillatorStateModel)
    return GaussianStateModel(osc.A, osc.Q, osc.b, osc.x0, osc.P0)
end

# canonical way to observe from an oscillator model
function LinearDynamicalSystem(osc::OscillatorStateModel, R = 1.0)
    noscs = length(osc.a)

    C = repeat([1. 0.], outer = (1,noscs))
    if ! (R isa Matrix)
        R = [R;;]
    end

    gsm = GaussianStateModel(osc)
    gom = GaussianObservationModel(C, R, zeros(1))

    return LinearDynamicalSystem(gsm, gom, 2noscs, 1, [true, true, true, true, false, true])
end
