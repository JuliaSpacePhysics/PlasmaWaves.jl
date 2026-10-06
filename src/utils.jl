_transpose(x) = x
_transpose(x::AbstractMatrix) = x'

function _transpose(f, X, args...; dim = 1, kw...)
    in = dim == 1 ? X : X'
    out = f(in, args...; kw...)
    return dim == 1 ? out : map(_transpose, out)
end
