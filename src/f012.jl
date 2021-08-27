using CSV
using Cairo
using Fontconfig
using Gadfly
using Printf
using DataFrames
using Compose
using Colors
using ColorSchemes
using Contour

shapes = [
    Gadfly.Shape.utriangle,
    Gadfly.Shape.circle,
    Gadfly.Shape.cross,
    Gadfly.Shape.diamond,
    Gadfly.Shape.octagon,
]

try
    Gadfly.pop_theme()
catch
end;
theTheme = Theme(
    alphas = [0.1],
    discrete_highlight_color = c -> RGBA{Float32}(c.r, c.g, c.b, 1),
    point_shapes = shapes,
    highlight_width = 0.05pt,

    major_label_font_size = 8pt,
    point_size = 2pt,
    key_swatch_color = "black",
    plot_padding = [0.1inch, 0.0inch, 0.1inch, -0.1inch],
    key_label_font_size = 4pt,
    key_title_font_size = 5pt,
)
Gadfly.push_theme(theTheme)

dfall = CSV.read("Data/Pexperiments.csv", DataFrame, header = true)
dfT = CSV.read("Data/transfer.csv", DataFrame, header = true)

layers = []
push!(layers, layer(dfT, x = :x, y = :y, Geom.line, Theme(default_color = "blue")))
push!(
    layers,
    layer(
        dfall,
        x = :ratio,
        y = :normalized,
        color = :legend3,
        shape = :legend2,
        Geom.point,
        Theme(
            alphas = [0.1],
            discrete_highlight_color = c -> RGBA{Float32}(c.r, c.g, c.b, 1),
            point_shapes = shapes,
            highlight_width = 0.4pt,
            point_size = 2.5pt,
        ),
    ),
)

guides = []
push!(guides, Guide.xlabel("z<sub>DMA</sub>/z<sup>s</sup> (-)"))
push!(guides, Guide.ylabel("Fraction transmitted (-)"))
push!(guides, Guide.YTicks(ticks = collect(0:0.2:1.0)))
push!(guides, Guide.XTicks(ticks = collect(0:0.2:1.6)))
push!(guides, Guide.shapekey(title = "Shapekey"))
push!(guides, Guide.colorkey("Colorkey"))

PJ_palette_paper = [
	colorant"rgb(128,0,0)",
	colorant"rgb(0,0,139)",
	colorant"rgb(255,140,0)",
	colorant"rgb(0,100,0)"
]

scales = []
push!(scales, Scale.color_discrete_manual(PJ_palette_paper...))
push!(scales, Scale.x_continuous())

coords = []
push!(coords, Coord.cartesian(ymin = 0, ymax = 1, xmax = 1.6, xmin = 0))

pre = Gadfly.plot(layers..., guides..., scales..., coords..., theTheme)
img = PNG("Figures/f012.png", dpi = 600, 3.1inch, 2.2inch)
draw(img, pre)
