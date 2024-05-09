 # heatmap
theme_methy_heatmap <- ggplot2::theme_bw() +
        ggplot2::theme(
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )

heatmap_fill_scale <- ggplot2::scale_fill_gradient2(low = "#CCCC3D", mid = "#53857E", high = "#1933B2", midpoint = 0.5)
heatmap_col_scale <- ggplot2::scale_color_gradient2(low = "#CCCC3D", mid = "#53857E", high = "#1933B2", midpoint = 0.5)

# Old scheme
# heatmap_fill_scale <- scico::scale_fill_scico(palette = "imola", direction = -1)
# heatmap_col_scale <- scico::scale_colour_scico(palette = "imola", direction = -1)

# Future release candidate
# heatmap_fill_scale <- scico::scale_fill_scico(palette = "roma", direction = -1)
# heatmap_col_scale <- scico::scale_colour_scico(palette = "roma", direction = -1)
