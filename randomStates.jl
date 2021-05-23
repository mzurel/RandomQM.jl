using LinearAlgebra

include("classicalEnsembles.jl")

"""
    randomPureState(D, format)

A random pure state of a D-level quantum system.  If format='DM' (default) then
the state returned is represented as a density matrix.  If format='Ket', the
state returned is represented as a vector.  If real=true, then a real state is
returned.
"""
function randomPureState(D::Int; format::String="DM", real::Bool=false)
    if real
        ψ = randn(D)
    else
        ψ = randn(D) + randn(D)im
    end
    ψ /= norm(ψ)

    if format == "Ket"
        return ψ
    elseif format == "DM"
        return ψ * transpose(ψ)
    else
        throw(ArgumentError("Argument format must be either 'Ket' or 'DM' (for density matrix)"))
    end
end

"""
    randomMixedState(D)

A random mixed state of a D-level quantum system represented as a density matrix.
If measure="HS" then the state is distributed uniformly according to the Hilbert-
Schmidt measure and if measure="Bures" then the state is distributed uniformly
according to the Bures measure.  If real=true then a real density matrix is
returned.
[Phys. Lett. A 373 320-324 (2009); J. Phys. A: Math. Gen. 34 7111-7125 (2001); J. Phys. A: 43 055302 (2010).]
"""
function randomMixedState(D::Int; real::Bool=false, measure::String="HS")
    if measure == "HS"
        if real
            X = realGinibre(D, D+1)
        else
            X = complexGinibre(D)
        end
        X = X * transpose(conj.(X))
        return X / tr(X)

    elseif measure == "Bures"
        if real
            X = realGinibre(D, D+1)
            V = circularOrthogonal(D)
            ρ = sqrt((I+V)*adjoint(I+V)) * X * transpose(X) * sqrt((I+adjoint(V))*adjoint(I+adjoint(V)))
            return ρ / tr(ρ)
        else
            X = complexGinibre(D)
            U = circularUnitary(D)
            ρ = (I + U) * X * adjoint(X) * (I + adjoint(U))
            return ρ / tr(ρ)
        end
    else
        throw(ArgumentError("Argument measure must be either 'HS' (for Hilbert-Schmidt) or 'Bures'"))
    end
end
