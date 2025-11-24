_transpose(x) = x
_transpose(x::AbstractMatrix) = x'

# Helper function to transpose the input and result for consistency
function _transpose(f, X, args...; dim = nothing, kw...)
    dim = @something dim 1
    in = dim == 1 ? X : X'
    out = f(in, args...; kw...)
    return dim == 1 ? out : map(_transpose, out)
end
