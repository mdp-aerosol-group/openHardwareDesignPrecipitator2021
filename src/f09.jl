using CSV
using DataFrames
using Cairo
using Fontconfig
using Gadfly
using Compose
using Printf
using NumericIO
using Lazy


try
    Gadfly.pop_theme()
catch
end;
theme = style(
    key_title_font = "PT Sans",
    major_label_font = "PT Sans",
    key_label_font = "PT Sans",
    key_label_font_size = 6pt,
    key_title_font_size = 0pt,
    point_size = 2.0pt,
    major_label_font_size = 8pt,
    line_width = 1pt,
    highlight_width = 6pt,
    plot_padding = [0.1inch, 0.2inch, 0.2inch, -0.1inch],
)
Gadfly.push_theme(theme)

t50 = 1 - (0.5^(1 / 3))
include("Scripts/precipitator.jl")
const Λ = Precipitator(0.003175, 0.007874, 0.130, 298.15, 1e5, 1lpm)
x = range(10, stop = 3000, length = 100)
voltage = repeat(x, 10)[:]
zstar = vtoz(Λ, voltage)

diameter = map(x -> ztod(Λ, x), zstar)
transmission50_1 = map(x -> ztod(Λ, x), zstar .* t50)

const Λ1 = Precipitator(0.003175, 0.007874, 0.215, 298.15, 1e5, 0.5lpm)
zstar05 = vtoz(Λ1, voltage)
diameter05 = map(x -> ztod(Λ1, x), zstar05)
transmission50_05 = map(x -> ztod(Λ1, x), zstar05 .* t50)

f = "1 L/min"
ff = "50% 1 L/min"
f05 = "0.5 L/min"
ff05 = "50% 0.5 L/min"
df1 = DataFrame(
    v = voltage,
    d = diameter * 1e9,
    tr50_1 = transmission50_1 * 1e9,
    Legend = [f for i = 1:size(voltage, 1)],
    Legend2 = [ff for i = 1:size(voltage, 1)],
)
df05 = DataFrame(
    v = voltage,
    d = diameter05 * 1e9,
    tr50_05 = transmission50_05 * 1e9,
    Legend = [f05 for i = 1:size(voltage, 1)],
    Legend2 = [ff05 for i = 1:size(voltage, 1)],
)


PJ_palette_paper =
    [colorant"rgb(139,58,58)" colorant"rgb(233, 150, 122)" colorant"rgb(46, 139, 87)" colorant"rgb(155,205,155)"]

fn1(x) = (x in 1:3) ? string(@sprintf("%d", 10^x)) : ""

gengrid(r) = [vcat(map(x -> x:x:9x, r)...); r[end] * 10]

xticks = log10.(gengrid([10, 10, 100, 1000]))
yticks = log10.(gengrid([10, 10, 100]))

cut = plot(
    layer(df1, x = :v, y = :d, color = :Legend, Geom.line),
    layer(
        df1,
        x = :v,
        y = :tr50_1,
        color = :Legend2,
        Geom.line,
        Theme(line_style = [:dash]),
    ),
    layer(df05, x = :v, y = :d, color = :Legend, Geom.line),
    layer(
        df05,
        x = :v,
        y = :tr50_05,
        color = :Legend2,
        Geom.line,
        Theme(line_style = [:dash]),
    ),
    Guide.xlabel("Voltage (V)"),
    Guide.ylabel("Diameter (nm)", orientation = :vertical),
    Scale.color_discrete_manual(PJ_palette_paper...),
    Scale.x_log10(labels = fn1),
    Scale.y_log10(labels = fn1),
    Guide.xticks(ticks = xticks),
    Guide.yticks(ticks = yticks),
    Guide.colorkey,
    Coord.cartesian(ymax = log10(1000.150), xmax = log10(3000)),
)

img = PNG("Figures/f09.png", dpi = 600, 3.25inch, 2inch);
draw(img, cut)
