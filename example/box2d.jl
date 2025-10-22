using BoundaryIntegral
using Plots

dir = @__DIR__

nquads = [2, 4, 8, 16, 32, 64]

ds = []
for d in [0.5, 0.1, 0.01]
    src = (1 - d, 0.0)
    ts = Float64[]
    for n_quad in nquads
        box = uniform_box2d(n_quad)
        t = 0.0
        for panel in box.panels
            for (point, weight) in zip(panel.points, panel.weights)
                t += laplace2d_doublelayer(src, point, panel.normal) * weight
            end
        end
        push!(ts, t)
    end
    push!(ds, ts)
end

scatter(nquads, log10.(abs.(ds[1] .- 1.0)), label="d=0.5")
scatter!(nquads, log10.(abs.(ds[2] .- 1.0)), label="d=0.1")
scatter!(nquads, log10.(abs.(ds[3] .- 1.0)), label="d=0.01")

savefig(joinpath(dir, "box2d_doublelayer_error.svg"))