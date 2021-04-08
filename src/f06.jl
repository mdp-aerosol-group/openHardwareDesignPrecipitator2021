using Cairo
using Fontconfig
using Gadfly
using Underscores
using CSV
using DataFrames
using Lazy
using MLStyle
using Interpolations

include("Scripts/precipitator.jl")
const Λ = Precipitator(0.25inches / 2, 0.75inches / 2, 12inches, 298.15, 1e5, 1lpm)
include("Scripts/reshape_fields.jl")

df, fx, fr, fθ = interpolated_field_functions()

function velocityprofile(Λ)
    r = range(Λ.r₁, stop = Λ.r₂, length = 100)
    v = @_ map(ux(Λ, _), r)
    DataFrame(r = 1000 * r, v = v, label = "Analytical")
end

function interpolated_profile(f, x, ψ)
    r = range(Λ.r₁, stop = Λ.r₂, length = 100)
    v = @_ map(f(x, _, ψ), r)
    angle = @> round(ψ * 180.0 / π, digits = 0) Int
    DataFrame(r = 1000 * r, v = v, label = "Interpolated", θ = angle)
end

function panel2()
    colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]
    xs = @>> df x -> x[!, :x] unique sort x -> x[80]
    doo = @>> df filter(x -> (x[:x] == xs)) x ->
        DataFrame(r = 1000 * x[!, :r], v = x[!, :vx], label = "Simulation")
    l1 = @> velocityprofile(Λ) layer(x = :r, y = :v, color = :label, Geom.line)
    l2 = @> interpolated_profile(fx, xs, 0) layer(x = :r, y = :v, color = :label, Geom.line)
    l3 = @> doo layer(x = doo[!, :r], y = doo[!, :v], color = :label, Geom.point)
    xl1 = [2 / 16, 3 / 16, 4 / 16, 5 / 16, 6 / 16]
    xl = ["2/16", "3/16", "4/16", "5/16", "6/16"]
    lfun(x) = xl[argmin(abs.(x ./ inches ./ 1000 .- xl1))]
    xticks =
        @> range(Λ.r₁, stop = Λ.r₂, length = 5) collect x -> round.(1000 * x; digits = 2)
    xlab = round(xs, digits = 3)

    plot(
        l1,
        l3,
        Guide.xlabel("Radial direction (mm)"),
        Guide.ylabel("v (m/s)"),
        Guide.title("x = $xlab (m)"),
        Guide.yticks(ticks = collect(0:0.02:0.1)),
        Guide.xticks(ticks = xticks),
        Guide.colorkey(title = ""),
        Scale.x_continuous(labels = lfun),
        Scale.color_discrete_manual(colors...),
        Coord.Cartesian(xmin = Λ.r₁ * 1000, xmax = Λ.r₂ * 1000),
        Theme(plot_padding = [-5mm, 2mm, 2mm, 2mm]),
    )
end

function panel1(posx)
    angles = -pi:pi/16:0
    df = @_ mapfoldl(interpolated_profile(fx, posx, _), vcat, angles)
    l = @> df layer(x = df[!, :r], y = :v, color = :θ, Geom.line)
    xl1 = [2 / 16, 3 / 16, 4 / 16, 5 / 16, 6 / 16]
    xl = ["2/16", "3/16", "4/16", "5/16", "6/16"]
    lfun(x) = xl[argmin(abs.(x ./ inches ./ 1000 .- xl1))]
    xticks =
        @> range(Λ.r₁, stop = Λ.r₂, length = 5) collect x -> round.(1000 * x; digits = 2)

    plot(
        l,
        Scale.color_continuous(minvalue = -180, maxvalue = 0),
        Guide.xlabel("Radial direction (inches)"),
        Guide.ylabel("v (m/s)"),
        Guide.title("x = $posx (m)"),
        Guide.yticks(ticks = collect(0:0.02:0.14)),
        Guide.xticks(ticks = xticks),
        Scale.x_continuous(labels = lfun),
        Coord.Cartesian(xmin = Λ.r₁ * 1000, xmax = Λ.r₂ * 1000),
        Theme(plot_padding = [0mm, 5mm, 2mm, 2mm]),
    )
end

set_default_plot_size(25cm, 8cm)
p = hstack(panel1(0.015), panel2(), panel1(0.29))
Gadfly.draw(PNG("Figures/f06.png", dpi = 600), p)
