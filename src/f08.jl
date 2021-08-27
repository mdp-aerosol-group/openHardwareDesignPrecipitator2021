using CSV
using DataFrames
using Cairo
using Fontconfig
using Gadfly
using Statistics
using Lazy
using MLStyle
using Underscores
using ProgressMeter
using Interpolations
using DifferentialEquations
using ParameterizedFunctions
using NumericIO

try
    Gadfly.pop_theme()
catch
end;

theTheme = Theme(
    highlight_width = 0.05pt,

    major_label_font_size = 8pt,
    line_width = 1pt,
    key_swatch_color = "black",
    plot_padding = [0.1inch, 0.0inch, 0.1inch, -0.1inch],
    key_label_font_size = 4pt,
    key_title_font_size = 5pt,
    point_size = 2pt,
)
Gadfly.push_theme(theTheme)

include("Scripts/precipitator.jl")
const Λ = Precipitator(0.25inches / 2, 0.75inches / 2, 12inches, 298.15, 1e5, 1lpm)

ODEsimple = @ode_def_bare begin
    dx = ux(Λ, r)
    dr = -zp * Er(Λ, v, r)
end Λ v zp

function trajectorymodel(posr, zp)
    function condition(u, t, integrator)
        terminate = (u[1] > Λ.l) || (u[2] < Λ.r₁) || (u[2] > Λ.r₂)
        out = terminate ? 0.0 : 1.0
    end
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)
    prob = ODEProblem(ODEsimple, [0.0, posr], (0.0, 60), [Λ v zp])
    sol = solve(
        prob,
        AutoTsit5(Rosenbrock23()),
        abstol = 1e-17,
        reltol = 1e-17,
        callback = cb,
    ) # error is here
    return mapfoldl(x -> x, hcat, sol.u)
end

Dp = 100e-9
v = dtov(Λ, Dp)
zp1 = dtoz(Λ, Dp)
zp2 = 0.5 * zp1

pos1 = trajectorymodel(0.75inches / 2.0, zp1)
pos2 = trajectorymodel(0.75inches / 2.0, 0.5 * zp1)

pos3 = trajectorymodel(0.75inches / 3.0, zp1)
pos4 = trajectorymodel(0.75inches / 3.0, 0.5 * zp1)

pos5 = trajectorymodel(0.75inches / 4.8, zp1)
pos6 = trajectorymodel(0.75inches / 4.8, 0.5 * zp1)

df1 = DataFrame(
    x = (Λ.l .- pos1[1, :]) .* 39.3701,
    r = (pos1[2, :]) .* 39.3701,
    zp = formatted(zp1, :ENG, ndigits = 3),
)
df2 = DataFrame(
    x = (Λ.l .- pos2[1, :]) .* 39.3701,
    r = (pos2[2, :]) .* 39.3701,
    zp = formatted(zp2, :ENG, ndigits = 3),
)
df3 = DataFrame(
    x = (Λ.l .- pos3[1, :]) .* 39.3701,
    r = (pos3[2, :]) .* 39.3701,
    zp = formatted(zp1, :ENG, ndigits = 3),
)
df4 = DataFrame(
    x = (Λ.l .- pos4[1, :]) .* 39.3701,
    r = (pos4[2, :]) .* 39.3701,
    zp = formatted(zp2, :ENG, ndigits = 3),
)
df5 = DataFrame(
    x = (Λ.l .- pos5[1, :]) .* 39.3701,
    r = (pos5[2, :]) .* 39.3701,
    zp = formatted(zp1, :ENG, ndigits = 3),
)
df6 = DataFrame(
    x = (Λ.l .- pos6[1, :]) .* 39.3701,
    r = (pos6[2, :]) .* 39.3701,
    zp = formatted(zp2, :ENG, ndigits = 3),
)
PJ_palette_paper = [
    colorant"rgb(0, 0, 139)",
    colorant"rgb(204,204,0)",
    colorant"rgb(0,0,138)",
    colorant"rgb(203,203,0)",
    colorant"rgb(0,0,137)",
    colorant"rgb(202,202,0)",
]
ytick = [0.125, 0.1875, 0.25, 0.3125, 0.375]
labels = Dict(zip(ytick, ["2/16", "3/16", "4/16", "5/16", "6/16"]))
Colorkey = ["xx"]

yticksfn = ["2/16", "3/16", "4/16", "5/16", "6/16"]
cut = plot(
    layer(df1, x = (df1[!, 1] * -1) .+ 12, y = :r, color = :zp, Geom.line),
    layer(df2, x = (df2[!, 1] * -1) .+ 12, y = :r, color = :zp, Geom.line),
    layer(df3, x = (df3[!, 1] * -1) .+ 12, y = :r, color = :zp, Geom.line),
    layer(df4, x = (df4[!, 1] * -1) .+ 12, y = :r, color = :zp, Geom.line),
    layer(df5, x = (df5[!, 1] * -1) .+ 12, y = :r, color = :zp, Geom.line),
    layer(df6, x = (df6[!, 1] * -1) .+ 12, y = :r, color = :zp, Geom.line),
    Guide.xlabel("Stream direction (inches)"),
    Guide.ylabel("Radial direction (inches)", orientation = :vertical),
    Scale.color_discrete_manual(PJ_palette_paper...),
    Guide.xticks(ticks = [0:2:12;]),
    Guide.yticks(ticks = ytick),
    Scale.y_continuous(labels = x -> labels[x]),
    theTheme,
    Guide.colorkey(title = "zp m<sup>2</sup> V<sup>-1</sup> s<sup>-1</sup>"),
    Coord.cartesian(xmin = 0, xmax = 13),
)

img = PNG("Figures/f08.png", dpi = 600, 3.1inch, 2.2inch);
draw(img, cut)
