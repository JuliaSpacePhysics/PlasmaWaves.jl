function workload(S = Float64)
    X = rand(S, 1000, 3)
    return map((Means(), Samson(), Santolik())) do method
        wavpol(X, S(1.0); method)
        twavpol(X; method)
        twavpol(permutedims(X); dim = 2, method)
    end
end


@setup_workload begin
    @compile_workload begin
        workload()
        workload(Float32)
    end
end
