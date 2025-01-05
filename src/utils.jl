catCol(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
# Simeon Schaub https://discourse.julialang.org/t/adding-rows-to-a-matrix-dynamically/52984
