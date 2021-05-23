using LinearAlgebra

"""
    circularOrthogonal(n)

A n×n real orthogonal matrix distributed according to the circular orthogonal
ensemble. [arXiv:math-ph/0609050]
"""
function circularOrthogonal(n)
    Z = randn((n,n))
    Q, R = qr(Z)
    D = diag(R)
    D = diagm(D ./ abs.(D))
    return Q * D
end

"""
    circularUnitary(n)

A n×n complex unitary matrix distributed according to the circular unitary
ensemble.  I.e. a random unitary distributed uniformly according to the Harr
measure.  [arXiv:math-ph/0609050]
"""
function circularUnitary(n)
    Z = (randn((n,n)) + randn((n,n))im) ./ sqrt(2)
    Q, R = qr(Z)
    D = diag(R)
    D = diagm(D ./ abs.(D))
    return Q * D
end

"""
    GaussianOrthogonal(n)

A n×n real matrix distributed according to the Ginibre ensemble.
"""
function realGinibre(n)
    A = randn((n,n))
    return A
end

"""
    complexGinibre(n)

A n×n complex matrix distributed according to the Ginibre ensemble.
"""
function complexGinibre(n)
    return (randn((n,n)) + randn((n,n))) ./ sqrt(2)
end
