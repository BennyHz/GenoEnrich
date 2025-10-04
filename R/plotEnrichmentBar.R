#' Create a bar plot for genomic coverage
#'
#' Shows total state size and overlap with target region as grouped or stacked bars.
#'
#' @param result_df Output from regionEnrich()
#' @param color_scheme Color scheme: "teal-rose", "blue-red", etc., or list(target=..., background=...)
#' @param stacked Whether to use stacked bars (default: FALSE = grouped)
#' @param x_order Character vector of State order
#' @param title Plot title
#' @return ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual scale_y_continuous labs theme_minimal theme
#' @importFrom tidyr pivot_longer
#' @importFrom scales label_comma
plotEnrichmentBar <- function(
    result_df,
    color_scheme = "teal-rose",
    stacked = FALSE,
    x_order = NULL,
    title = "Genomic Coverage by Chromatin State"
) {
  required_cols <- c("State", "Overlap_BP", "Total_State_BP")
  if (!all(required_cols %in% names(result_df))) {
    stop("输入数据必须包含列: ", paste(required_cols, collapse = ", "))
  }

  # === 颜色方案处理（修复 switch 问题）===
  if (is.character(color_scheme) && length(color_scheme) == 1) {
    # 预设名称
    color_settings <- switch(
      color_scheme,
      "blue-red"     = list(target = "#e41a1c", background = "#4575b4"),
      "green-purple" = list(target = "#8e5ea2", background = "#66c2a5"),
      "teal-rose"    = list(target = "#a46b6b", background = "#6b8ca4"),
      "viridis"      = list(target = "#fde725", background = "#440154"),
      "magma"        = list(target = "#fcc52d", background = "#000004"),
      stop("未知的 color_scheme: '", color_scheme, "'. 支持: blue-red, green-purple, teal-rose, viridis, magma")
    )
  } else if (is.list(color_scheme) &&
             all(c("target", "background") %in% names(color_scheme))) {
    # 用户自定义 list
    color_settings <- color_scheme
  } else {
    stop("color_scheme 必须是预设字符串（如 'teal-rose'）或 list(target = '#...', background = '#...')")
  }

  # 准备数据
  plot_data <- result_df[, required_cols]
  plot_data$State <- as.character(plot_data$State)
  plot_data$NonOverlap_BP <- pmax(0, plot_data$Total_State_BP - plot_data$Overlap_BP)

  if (is.null(x_order)) {
    x_order <- sort(unique(plot_data$State))
  }
  x_order <- x_order[x_order %in% unique(plot_data$State)]
  plot_data <- plot_data[plot_data$State %in% x_order, ]

  if (stacked) {
    plot_long <- tidyr::pivot_longer(
      plot_data,
      cols = c(Overlap_BP, NonOverlap_BP),
      names_to = "Category",
      values_to = "Size_BP"
    )
    plot_long$Category <- factor(
      plot_long$Category,
      levels = c("NonOverlap_BP", "Overlap_BP"),
      labels = c("Non-overlap", "Overlap with Target")
    )
    position <- "stack"
    legend_title <- "Region Type"
  } else {
    plot_long <- tidyr::pivot_longer(
      plot_data,
      cols = c(Total_State_BP, Overlap_BP),
      names_to = "Category",
      values_to = "Size_BP"
    )
    plot_long$Category <- factor(
      plot_long$Category,
      levels = c("Total_State_BP", "Overlap_BP"),
      labels = c("Total State Size", "Overlap with Target")
    )
    position <- "dodge"
    legend_title <- "Category"
  }

  # 绘图
  p <- ggplot2::ggplot(plot_long,
                       aes(x = factor(State, levels = x_order),
                           y = Size_BP,
                           fill = Category)) +
    ggplot2::geom_col(
      position = position,
      width = 0.7,
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::scale_fill_manual(
      values = if (stacked) {
        c("Non-overlap" = color_settings$background,
          "Overlap with Target" = color_settings$target)
      } else {
        c("Total State Size" = color_settings$background,
          "Overlap with Target" = color_settings$target)
      },
      name = legend_title
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::labs(
      title = title,
      x = "Chromatin States",
      y = "Genomic Size (bp)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9)
    )

  return(p)
}
