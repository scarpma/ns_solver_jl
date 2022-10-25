using Plots
include("./IO.jl")
include("./FiniteDifference.jl")
using .IO
using .FiniteDifference

function test()
    # test d1
    res = Float64[]
    for i=2:12
        push!(res, FiniteDifference.test_d2(i^2))
    end
    p = plot([i^2 for i in 2:12], res, xaxis=:log, yaxis=:log)
    x = [i for i in 2:12]
    plot!(x, x.^(-2.), title="d1 (x1**3 + x2**3)")
    display(p)

    # test d2
    res = Float64[]
    for i=2:12
        push!(res, FiniteDifference.test_d2(i^2))
    end
    p = plot([i^2 for i in 2:12], res, xaxis=:log, yaxis=:log)
    x = [i for i in 2:12]
    plot!(x, x.^(-2.), title="d2 (x1**3 + x2**3)")
    display(p)

    # test dd1
    res = Float64[]
    for i=2:12
        push!(res, FiniteDifference.test_dd1(i^2))
    end
    p = plot([i^2 for i in 2:12], res, xaxis=:log, yaxis=:log)
    x = [i for i in 2:12]
    plot!(x, x.^(-2.), title="dd1 (x1**3 + x2**3)")
    display(p)

    # test dd1
    res = Float64[]
    for i=2:12
        push!(res, FiniteDifference.test_dd2(i^2))
    end
    p = plot([i^2 for i in 2:12], res, xaxis=:log, yaxis=:log)
    x = [i for i in 2:12]
    plot!(x, x.^(-2.), title="dd2 (x1**3 + x2**3)")
    display(p)

end


test()
