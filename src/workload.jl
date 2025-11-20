function workload(S = Float64)
    X = rand(S, 1000, 3)
    return wavpol(X, S(1.0))
end


@setup_workload begin
    @compile_workload begin
        workload()
        workload(Float32)
    end
end
