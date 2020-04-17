# custom ggplot theme


custom_theme <- function() {
  theme_bw() +
    theme(
      axis.text = element_text(size=8),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title = element_text(size =8),
      panel.grid.major = element_line(color = 'gray'),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = '#ffffff'),
      strip.background = element_rect(fill = 'gray', color = 'black', size =0.5),
      strip.text = element_text(face = 'bold', size = 8, color = 'white'),
      legend.position = 'right',
      legend.justification = 'center',
      legend.background = element_blank(),
      panel.border = element_rect(color = 'grey5', fill = NA, size = 0.5))
}

theme_set(custom_theme())
