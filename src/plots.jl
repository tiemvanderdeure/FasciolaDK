using Makie
# Plot Denmark with bornholm on a inset
function plot_denmark!(ax, ras::Raster; kw...)
    dk = copy(ras[X = 7 .. 12.8])
    bornholm = ras[Extent(X = (14.5, 15.2), Y = (54.8, 55.5))]
    dk[Extent(X = (11.8, 12.5), Y = (56.8, 57.5))] .= parent(bornholm)
    poly!(ax, Rect2(11.85, 56.8,0.75,0.75), color = :transparent, strokecolor = :black, strokewidth = 1)
    plot!(ax,dk; kw...)
end
