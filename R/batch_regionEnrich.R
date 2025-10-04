#' Batch enrichment analysis for multiple target regions
#'
#' Runs regionEnrich on a list of target GRanges and combines results.
#'
#' @param target_gr_list List of GRanges objects (e.g., list(NAD = nad_gr, interNAD = inter_gr))
#' @param feature_gr GRanges with 'state' column (same for all targets)
#' @param target_labels Character vector of labels (default: names of target_gr_list)
#' @param background_gr Optional background GRanges (passed to regionEnrich)
#' @param min_overlap_bp Minimum overlap bp (passed to regionEnrich)
#' @param sort_by Sorting method (passed to regionEnrich)
#' @param add_log10p Whether to add log10 p-values (passed to regionEnrich)
#' @param sig_method Significance method (passed to regionEnrich)
#' @return Combined data.frame with all results
#' @export
batch_regionEnrich <- function(
    target_gr_list,
    feature_gr,
    target_labels = names(target_gr_list),
    background_gr = NULL,
    min_overlap_bp = 1,
    sort_by = "Enrichment_Ratio",
    add_log10p = TRUE,
    sig_method = "FDR"
) {
  # 参数检查
  if (!is.list(target_gr_list) || length(target_gr_list) == 0)
    stop("target_gr_list 必须是非空列表")

  if (!inherits(feature_gr, "GRanges"))
    stop("feature_gr 必须是 GRanges 对象")

  if (!"state" %in% names(mcols(feature_gr)))
    stop("feature_gr 的 mcols 必须包含 'state' 列")

  if (is.null(target_labels)) {
    target_labels <- paste0("Target_", seq_along(target_gr_list))
  } else {
    if (length(target_labels) != length(target_gr_list))
      stop("target_labels 长度必须与 target_gr_list 一致")
  }

  # 确保所有 target_gr 是 GRanges
  for (i in seq_along(target_gr_list)) {
    if (!inherits(target_gr_list[[i]], "GRanges")) {
      stop("target_gr_list[[", i, "]] 不是 GRanges 对象")
    }
  }

  cat("开始批量富集分析，共", length(target_gr_list), "个目标区域\n")

  # 存储所有结果
  all_results <- list()

  # 循环每个目标区域
  for (i in seq_along(target_gr_list)) {
    label <- target_labels[i]
    target_gr <- target_gr_list[[i]]

    cat("\n>>> 分析目标区域:", label, "\n")

    # 调用你已有的 regionEnrich 函数（不返回 overlap）
    res_df <- regionEnrich(
      target_gr = target_gr,
      feature_gr = feature_gr,
      target_label = label,
      sig_method = sig_method,
      background_gr = background_gr,
      min_overlap_bp = min_overlap_bp,
      return_overlap_gr = FALSE,   # 批量分析通常不需要重叠区域
      sort_by = sort_by,
      add_log10p = add_log10p
    )

    all_results[[i]] <- res_df
  }

  # 合并结果
  final_df <- do.call(rbind, all_results)
  rownames(final_df) <- NULL

  # 按 target_label 和排序字段整理（可选）
  final_df <- final_df[order(final_df$tissue_state, -final_df[[sort_by]]), ]
  rownames(final_df) <- NULL

  cat("\n✅ 批量富集完成！共", nrow(final_df), "行结果\n")
  return(final_df)
}
