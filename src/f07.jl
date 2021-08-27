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

try
    Gadfly.pop_theme()
catch
end;
theTheme = Theme(
    highlight_width = 0.05pt,

    major_label_font_size = 8pt,
    default_color = "black",
    point_size = 2pt,
    key_swatch_color = "black",
    line_width = 1pt,
    plot_padding = [0.1inch, 0.1inch, 0.1inch, -0.1inch],
    key_label_font_size = 4pt,
    key_title_font_size = 5pt,
)
Gadfly.push_theme(theTheme)

include("Scripts/precipitator.jl")
const Λ = Precipitator(0.25inches / 2, 0.75inches / 2, 12inches, 298.15, 1e5, 1lpm)
include("Scripts/reshape_fields.jl")

df, fx, fr, fθ = interpolated_field_functions()

Dp = 200e-9
v = dtov(Λ, Dp)
zs = dtoz(Λ, Dp)
zzs = [0.0001; 0.01:0.01:0.1; 0.15:0.05:0.9; 0.91:0.01:1; 1.2] |> collect

ODEsimple = @ode_def_bare begin
    dx = ux(Λ, r)
    dr = -zp * Er(Λ, v, r)
end Λ v zp

function simplemodel(posr, zp)
    function condition(u, t, integrator)
        out = (u[1] > Λ.l) || (u[2] < Λ.r₁) || (u[2] > Λ.r₂) ? 0.0 : 1.0
    end
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)
    prob = ODEProblem(ODEsimple, [0.0, posr], (0.0, 60), [Λ v zp])
    sol = solve(prob, RK4(), abstol = 1e-7, reltol = 1e-7, callback = cb)
    ext = mapfoldl(x -> x, hcat, sol.u)
    T = (ext[1, end] > Λ.l) & (sol.t[end] < 60) ? 1.0 : 0.0
end

rrs = range(Λ.r₁, stop = Λ.r₂, length = 200)
rmid = (rrs[2:end] - rrs[1:end-1]) ./ 2 .+ rrs[1:end-1]
us = map(x -> ux(Λ, x), rmid)

frac = map(1:length(rrs)-1) do i
    pi * us[i] * (rrs[i+1]^2 - rrs[i]^2) ./ 1.6666666e-5
end

f(zzs) = @> map(x -> simplemodel(x, zzs * zs), rmid) (x -> x .* frac) sum
tr = @showprogress map(f, zzs)
df = DataFrame(zzs = zzs, T = tr, label = "Diffeq")
p = plot(
    df,
    x = :zzs,
    y = :T,
    Geom.line,
    Guide.xlabel("z/zˢ"),
    Guide.xticks(ticks = 0:0.2:1.2),
    Guide.yticks(ticks = 0:0.2:1),
    Guide.ylabel("Fraction transmitted (-)"),
    Coord.cartesian(xmin = 0, xmax = 1.2),
    theTheme,
)
set_default_plot_size(3.1inch, 2.2inch)
Gadfly.draw(PNG("Figures/f07.png", dpi = 600), p)
