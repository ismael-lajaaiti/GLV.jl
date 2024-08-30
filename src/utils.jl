"""
    offdiag(A)

Dictionnary of off-diagonal elements of a matrix `A`.
Keys are pairs of indices `(i, j)` and values are the corresponding elements of `A`.

```julia
A = [1 2; 3 4]
offdiag(A)
```
"""
function offdiag(A)
    @assert size(A, 1) == size(A, 2)
    S = size(A, 1)
    Dict([(i, j) => A[i, j] for i = 1:S, j = 1:S if i != j])
end
export offdiag
