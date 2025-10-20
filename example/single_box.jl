using BoundaryIntegral
using Plots

ds = []
for d in [1.0, 0.1, 0.01]
    src = (1 - d, 0.0, 0.0)
    ts = Float64[]
    for n_surface in 1:20
        box = uniform_box(n_surface)
        t = 0.0
        for i in 1:box.n
            t += coulomb_En(src, box.points[i], box.norms[i]) * box.weights[i]
        end
        push!(ts, t)
    end
    push!(ds, ts)
end

plot(1:20, log10.(abs.(ds[1] .- 1.0)), label="d=1.0")
plot!(1:20, log10.(abs.(ds[2] .- 1.0)), label="d=0.1")
plot!(1:20, log10.(abs.(ds[3] .- 1.0)), label="d=0.01")
savefig("single_box_error.png")