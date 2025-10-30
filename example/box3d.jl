using BoundaryIntegral
using Plots

ds = []
for d in [0.5, 0.1, 0.01]
    src = (1 - d, 0.0, 0.0)
    ts = Float64[]
    for n_interface in 1:20
        box = uniform_box3d(n_interface)
        t = 0.0
        for i in 1:box.n
            t += laplace3d_doublelayer(src, box.points[i], box.norms[i]) * box.weights[i]
        end
        push!(ts, t)
    end
    push!(ds, ts)
end

plot(1:20, log10.(abs.(ds[1] .- 1.0)), label="d=0.5")
plot!(1:20, log10.(abs.(ds[2] .- 1.0)), label="d=0.1")
plot!(1:20, log10.(abs.(ds[3] .- 1.0)), label="d=0.01")
savefig("single_box_error.svg")