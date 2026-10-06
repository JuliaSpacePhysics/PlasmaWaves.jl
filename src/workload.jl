function workload(S = Float64)
    X = rand(S, 1000, 3)
    wavpol(X, S(1.0))
    twavpol(X)
    twavpol_svd(X)
    return twavpol(permutedims(X); dim = 2), twavpol_svd(permutedims(X); dim = 2)
end


@setup_workload begin
    @compile_workload begin
        workload()
        workload(Float32)
    end
end
