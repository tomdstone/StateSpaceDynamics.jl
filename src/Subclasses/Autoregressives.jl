mutable struct AutoregressiveStateModel <: AbstractStateModel{Float64}
    coeff
    σ2
    A
    Q
    b
    x0
    P0
    function AutoregressiveStateModel(;
        coeff,
        σ2,
        P0=nothing,
        x0=nothing,
    )
        ord = length(coeff)

        A = vcat(reshape(coeff, (1,:)), hcat(I, zeros(ord - 1)))
        Q = diagm(vcat([σ2], zeros(ord - 1)))

        b  = zeros(ord)
        x0 = something(x0, zeros(ord))
        P0 = something(P0, Q)

        @assert length(x0) == ord       "Initial state vector is wrong size"
        @assert size(P0) == (ord, ord)  "Initial state covariance is wrong size"

        return new(reshape(coeff, :), σ2, A, Q, b, x0, P0)
    end
end

function Base.show(io::IO, arn::AutoregressiveStateModel; gap="")
    println(io, gap, "Autoregressive State Model:")
    println(io, gap, "---------------------------")
    println(io, gap, " State Parameters:")
    println(io, gap, "  order = $(order(arn))")
    println(io, gap, "  coeff = $(round.(arn.coeff, sigdigits=3))")
    println(io, gap, "  σ2    = $(round.(arn.σ2, digits=2))")
    println(io, gap, " Initial State: ")

    if order(arn) < 4
        println(io, gap, "  x0 = $(round.(arn.x0, digits=2))")
        println(io, gap, "  P0 = $(round.(arn.P0, digits=2))")
    else
        println(io, gap, "  size(x0) = ($(length(arn.x0)),)")
        println(io, gap, "  size(P0) = ($(size(arn.P0,1)), $(size(arn.P0,2)))")
    end

    return nothing
end

function order(arn::AutoregressiveStateModel)
    return length(arn.coeff)
end
