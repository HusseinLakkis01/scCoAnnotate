# function for a minimal theme for ggplot
# author: Selin Jessa


theme_min <- function(base_size = 11, base_family = "",
                      border_colour = "grey90",
                      border_size = 1) {
    
    theme_light(base_size = 11, base_family = "") +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = border_colour, size = border_size),
            axis.ticks = element_line(colour = border_colour),
            strip.background = element_rect(fill = NA, colour = NA),
            strip.text.x = element_text(colour = "black", size = rel(1.2)),
            strip.text.y = element_text(colour = "black", size = rel(1.2)),
            title = element_text(size = rel(0.9)),
            axis.text = element_text(colour = "black", size = rel(0.8)),
            axis.title = element_text(colour = "black", size = rel(0.9)),
            legend.title = element_text(colour = "black", size = rel(0.9)),
            legend.key.size = unit(0.9, "lines"),
            legend.text = element_text(size = rel(0.7), colour = "black"),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.background = element_rect(colour = NA, fill = NA)
        )
}
