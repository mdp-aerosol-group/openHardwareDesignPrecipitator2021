using Cairo, Fontconfig, Lazy, CSV, Dates, DataFrames, Gadfly
using Printf, Statistics, SpecialFunctions
using Compose, Printf, Colors
using TimeSeries

try
    Gadfly.pop_theme()
catch
end;
theTheme = Theme(
    highlight_width = 0.05pt,
    major_label_font = "PT Sans",
    minor_label_font = "PT Sans",
    key_label_font = "PT Sans",
    major_label_font_size = 8pt,
    line_width = 0.5pt,
    key_swatch_color = "black",
    plot_padding = [0.1inch, 0.1inch, 0.1inch, -0.1inch],
    key_label_font_size = 4pt,
    key_title_font_size = 0pt,
)

Gadfly.push_theme(theTheme)

dftime = CSV.read("Data/timeseries.csv", DataFrame, header = true)

layers = []
push!(layers, layer(dftime, x = :time, y = :cpc1, color = :legend, Geom.line))
push!(layers, layer(dftime, x = :time, y = :cpc2, color = :legend2, Geom.line))

push!(layers, layer(dftime, x = :time, y = :ratios, color = :legend4, Geom.line))
push!(layers, layer(dftime, x = :time, y = :mean, color = :legend3, Geom.line))

guides = []
guides = []
push!(guides, Guide.xlabel("Time steps (s)"))
push!(guides, Guide.ylabel("Concentration counts (cm<sup>-3</sup>)"))

PJ_palette_paper = [
	colorant"rgb(0,0,153)",
	colorant"rgb(0,0,0)", 
	colorant"rgb(0,153,0)", 
	colorant"rgb(153,0,0)"
]
scales = []

push!(scales, Scale.color_discrete_manual(PJ_palette_paper...))
push!(scales, Scale.x_continuous())

pre = Gadfly.plot(layers..., guides..., scales..., theTheme)

img = PNG("Figures/f08s.png", dpi = 600, 3.1inch, 2.2inch)
draw(img, pre)
