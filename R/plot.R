# The script visualize the learning process of the MVBLUP model:
# Notes:
# p1: the MVBLUP model iteratively learns the weights for different types of data
# p2: optimal weight changes throughout the iterative process of the MVBLUP model

plot_learning <- function(Training_Accuracy, Training_Weight, MVBLUP_information, output_path) {
  col1 <- colorRampPalette((pal_npg("nrc")(9)))(NP)
  Index_iter <- as.numeric(MVBLUP_information$values[n+1])
  dt_f <- Training_Accuracy
  dt_f <- na.omit(dt_f)
  dt_f$group <- factor(dt_f$group, levels = unique(dt_f$group))
  
  p1 <- ggplot(dt_f, aes(iter, Accuracy, color = group, group = group)) +
    geom_point(size = 1) +
    geom_line(size = 0.3) +
    theme_bw() +
    scale_color_manual(values = col1) +
    ylab("Training accuracy") +
    xlab("Iteration") +
    geom_vline(xintercept = Index_iter) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 12, vjust = 0.5),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_blank()
    )
  
  ggsave(paste0(output_path,"trait",trid,"-Training_Accuracy.png"), p1, width = 8.5, height = 2.5)
  
  dt_w <- Training_Weight
  dt_w <- na.omit(dt_w)
  Parameters <- paste0("View",1:n,"_W")
  dt_w$Parameter <- factor(Parameters[as.integer(dt_w$Parameter)], levels = Parameters)
  
  p2 <- ggplot(dt_w, aes(iter, Weight, color = Parameter, group = Parameter)) +
    geom_point(size = 1) +
    geom_line(size = 0.3) +
    theme_bw() +
    scale_color_manual(values = col1) +
    labs(color = "Parameter") +
    ylab("Optimal weight") +
    xlab("Iteration") +
    geom_vline(xintercept = Index_iter) +
    theme(legend.position = "top",
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, vjust = 0.5),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.direction = 'horizontal') +
    guides(color = guide_legend(nrow = 1, keywidth = 0.5, keyheight = 0.5))
  ggsave(paste0(output_path,"trait",trid,"-Optimal weight.png"), p2, width = 8.5, height = 2.5)
}