#' Create a bubble plot for enrichment results
#'
#' Visualizes enrichment results with bubble size = enrichment ratio,
#' color = odds ratio (log scale), and significance stars.
#'
#' @param enrichment_df Output from regionEnrich() or batch_regionEnrich()
#' @param domain_mapping Optional mapping from tissue_state to domain_type (named vector or function)
#' @param x_order Character vector of State order (default: sorted)
#' @param y_order Character vector of domain_type order (default: unique in data)
#' @param sort_y Whether to sort y_order alphabetically (default: FALSE)
#' @param sort_x_by Sort X-axis by "Enrichment_Ratio" or "Odds_Ratio" (default: NULL = no sorting)
#' @param sort_x_ref Reference domain_type for sorting (required if sort_x_by is not NULL)
#' @param color_scheme Color scheme: "teal-rose", "blue-red", "green-purple", "viridis", "magma", "plasma", or list(low=..., high=...)
#' @param title Main plot title
#' @param subtitle Subtitle (default: auto-generated)
#' @return ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_gradient2 scale_size_continuous geom_text labs theme_minimal theme
#' @importFrom dplyr mutate case_when
plotEnrichmentBubble <- function(
    enrichment_df,
    domain_mapping = NULL,
    x_order = NULL,
    y_order = NULL,
    sort_y = FALSE,
    sort_x_by = NULL,
    sort_x_ref = NULL,
    color_scheme = "teal-rose",
    title = "Chromatin State Enrichment",
    subtitle = NULL
) {
  required_cols <- c("tissue_state", "State", "Enrichment_Ratio", "Odds_Ratio", "FDR")
  if (!all(required_cols %in% names(enrichment_df))) {
    stop("输入数据必须包含列: ", paste(required_cols, collapse = ", "))
  }

  # 生成 domain_type
  plot_data <- enrichment_df
  if (is.null(domain_mapping)) {
    plot_data$domain_type <- as.character(plot_data$tissue_state)
  } else if (is.function(domain_mapping)) {
    plot_data$domain_type <- domain_mapping(plot_data$tissue_state)
  } else if (is.character(domain_mapping)) {
    plot_data$domain_type <- domain_mapping[as.character(plot_data$tissue_state)]
    unmatched <- is.na(plot_data$domain_type)
    if (any(unmatched)) {
      plot_data$domain_type[unmatched] <- as.character(plot_data$tissue_state[unmatched])
    }
  } else {
    stop("domain_mapping 必须是 NULL、函数或命名字符向量")
  }

  plot_data$domain_type <- as.character(plot_data$domain_type)
  plot_data$State <- as.character(plot_data$State)
  plot_data$Enrichment_Ratio_plot <- ifelse(plot_data$Enrichment_Ratio <= 0, 0.01, plot_data$Enrichment_Ratio)

  # Y 轴顺序
  if (is.null(y_order)) {
    y_order <- unique(plot_data$domain_type)
    if (sort_y) y_order <- sort(y_order)
  } else {
    if (sort_y) y_order <- sort(y_order)
  }
  y_order <- y_order[y_order %in% unique(plot_data$domain_type)]
  plot_data <- plot_data[plot_data$domain_type %in% y_order, ]

  # X 轴排序
  if (!is.null(sort_x_by)) {
    if (!sort_x_by %in% c("Enrichment_Ratio", "Odds_Ratio")) {
      stop("sort_x_by 必须是 'Enrichment_Ratio' 或 'Odds_Ratio'")
    }

    if (is.null(sort_x_ref)) {
      if (length(y_order) == 1) {
        sort_x_ref <- y_order[1]
      } else {
        stop("当 sort_x_by 非 NULL 时，必须指定 sort_x_ref（参考的 domain_type）")
      }
    }

    ref_row <- plot_data[plot_data$domain_type == sort_x_ref, ]
    if (nrow(ref_row) == 0) {
      stop("sort_x_ref '", sort_x_ref, "' 在数据中不存在")
    }

    ref_row <- ref_row[order(-ref_row[[sort_x_by]]), ]
    x_order <- ref_row$State
  } else {
    if (is.null(x_order)) {
      x_order <- sort(unique(plot_data$State))
    }
  }
  x_order <- x_order[x_order %in% unique(plot_data$State)]
  plot_data <- plot_data[plot_data$State %in% x_order & plot_data$domain_type %in% y_order, ]

  # 显著性标记
  plot_data$Significance <- ifelse(plot_data$FDR < 0.001, "***",
                                   ifelse(plot_data$FDR < 0.01, "**",
                                          ifelse(plot_data$FDR < 0.05, "*", "")))

  # === 修复 color_scheme 的 switch 问题 ===
  if (is.character(color_scheme) && length(color_scheme) == 1) {
    color_settings <- switch(
      color_scheme,
      "blue-red"     = list(low = "#4575b4", high = "#d73027"),
      "green-purple" = list(low = "#66c2a5", high = "#8e5ea2"),
      "teal-rose"    = list(low = "#6b8ca4", high = "#a46b6b"),
      "viridis"      = list(low = "#440154", high = "#fde725"),
      "magma"        = list(low = "#000004", high = "#fcc52d"),
      "plasma"       = list(low = "#0d0887", high = "#f0f921"),
      stop("未知 color_scheme: '", color_scheme, "'. 支持: blue-red, green-purple, teal-rose, viridis, magma, plasma")
    )
  } else if (is.list(color_scheme) && length(color_scheme) == 2) {
    # 允许 list(low = ..., high = ...) 或命名 list
    if (!("low" %in% names(color_scheme) && "high" %in% names(color_scheme))) {
      # 如果没命名，按位置取
      color_settings <- list(low = color_scheme[[1]], high = color_scheme[[2]])
    } else {
      color_settings <- color_scheme
    }
  } else {
    stop("color_scheme 必须是预设字符串或 list(low=..., high=...)")
  }

  # 默认副标题
  if (is.null(subtitle)) {
    subtitle <- "Size: enrichment ratio; Color: odds ratio (log scale)\n*** FDR<0.001, ** FDR<0.01, * FDR<0.05"
  }

  # 绘图
  p <- ggplot2::ggplot(plot_data,
                       aes(x = factor(State, levels = x_order),
                           y = factor(domain_type, levels = y_order))) +
    ggplot2::geom_point(
      aes(size = Enrichment_Ratio_plot, fill = Odds_Ratio),
      shape = 21, color = "black", alpha = 0.85, stroke = 0.5
    ) +
    ggplot2::scale_fill_gradient2(
      low = color_settings$low,
      mid = "white",
      high = color_settings$high,
      midpoint = 1,
      name = "Odds Ratio",
      trans = "log",
      na.value = "lightgray"
    ) +
    ggplot2::scale_size_continuous(range = c(2.5, 10), name = "Enrichment\nRatio") +
    ggplot2::geom_text(aes(label = Significance), size = 2.5, fontface = "bold", vjust = 0.8) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Chromatin States",
      y = "Domain Types"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  return(p)
}
