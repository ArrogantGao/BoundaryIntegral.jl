using BoundaryIntegral
import BoundaryIntegral as BI
using OMEinsum, LinearAlgebra

begin
    eps_box = 2.0
    src = (0.1, 0.1)

    box = BI.dielectric_box2d(4, 64, true, Float64, n_adapt = 2)
    lhs = BI.Lhs_dielectric_box2d(eps_box, box)
    rhs = BI.Rhs_dielectric_box2d(eps_box, box, src)

    x = BI.solve_lu(lhs, rhs);
    err = norm(lhs * x - rhs)

    g = BI.l2d_singlelayer_gi(box, x, 2.0, 32) + 1.0 / eps_box

    # identity check
    DT = BI.laplace2d_DT(box);

    @show ein"i, (ji, j, j) -> "(x, DT, BI.all_weights(box), ones(length(x))), sum(x) / 2

    self_check = ein"(ij, j), j -> "(DT, BI.all_weights(box), x);

    @show g, self_check, sum(x) / 2
end

xs = Float64[]
cs = Float64[]
i = 0
for n_panels in 1:box.n
    panel = box.panels[n_panels]
    for (weight, point) in zip(panel.weights, panel.points)
        i += 1
        if point[2] == 1.0
            push!(xs, point[1])
            push!(cs, x[i])
        end
    end
end

order_xs = sortperm(xs)
sxs = xs[order_xs]
scs = cs[order_xs]

scatter(sxs, scs)

scatter(log10.(1.0 .- sxs[end - 32:end]), log10.(scs[end - 32:end]))