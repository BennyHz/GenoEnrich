#' Perform genomic region enrichment analysis
#'
#' Computes enrichment of target regions (e.g., NAD domains) against annotated genomic states.
#' Uses hypergeometric test and reports enrichment ratio, odds ratio, and significance.
#'
#' @param target_gr GRanges object of target regions (e.g., NAD domains)
#' @param feature_gr GRanges object with a 'state' column in mcols (e.g., chromatin states)
#' @param target_label Label for the target region (e.g., "Hep_NAD")
#' @param sig_method Significance method: "FDR" (default) or "pValue"
#' @param background_gr Optional GRanges for background (default: total feature_gr)
#' @param min_overlap_bp Minimum overlap length (in bp) to count as overlap (default: 1)
#' @param return_overlap_gr Whether to return actual overlap GRanges (default: FALSE)
#' @param sort_by Sorting method: "Enrichment_Ratio", "Odds_Ratio", "P_value", "-log10(P_value)", or "-log10(FDR)"
#' @param add_log10p Whether to add -log10(P_value) and -log10(FDR) columns (default: TRUE)
#' @return data.frame with enrichment results, or list(result_df, overlaps) if return_overlap_gr=TRUE
#' @export
#' @importFrom GenomicRanges intersect width mcols
#' @importFrom stats phyper p.adjust
regionEnrich <- function(
    target_gr,
    feature_gr,
    target_label = "Target",
    sig_method = "FDR",
    background_gr = NULL,
    min_overlap_bp = 1,
    return_overlap_gr = FALSE,
    sort_by = "Enrichment_Ratio",
    add_log10p = TRUE
) {
  # 参数检查
  if (!inherits(target_gr, "GRanges"))
    stop("target_gr 必须是 GRanges 对象")
  if (!inherits(feature_gr, "GRanges"))
    stop("feature_gr 必须是 GRanges 对象")
  if (!"state" %in% names(mcols(feature_gr)))
    stop("feature_gr 的 mcols 必须包含 'state' 列")

  if (!sig_method %in% c("pValue", "FDR"))
    stop("sig_method 必须是 'pValue' 或 'FDR'")

  allowed_sort <- c("Enrichment_Ratio", "Odds_Ratio", "P_value", "-log10(P_value)", "-log10(FDR)")
  if (!sort_by %in% allowed_sort)
    stop("sort_by 必须是以下之一: ", paste(allowed_sort, collapse = ", "))

  if (min_overlap_bp < 1) min_overlap_bp <- 1

  # 设置背景区域
  if (is.null(background_gr)) {
    background_gr <- feature_gr
  } else {
    if (!inherits(background_gr, "GRanges"))
      stop("background_gr 必须是 GRanges 对象或 NULL")
  }

  # 提取状态并确保是字符型
  mcols(feature_gr)$state <- as.character(mcols(feature_gr)$state)
  unique_states <- unique(mcols(feature_gr)$state)
  cat("检测到染色质状态：", paste(unique_states, collapse = ", "), "\n")

  # 总长度
  total_target_bp <- sum(width(target_gr))
  total_background_bp <- sum(width(background_gr))
  cat("目标区域总长度:", format(total_target_bp, scientific = TRUE), "bp\n")
  cat("背景区域总长度:", format(total_background_bp, scientific = TRUE), "bp\n")

  # 存储结果和可选的重叠区域
  enrichment_results <- list()
  overlap_list <- list()

  for (state in unique_states) {
    cat("  处理状态:", state, "\n")

    current_state_gr <- feature_gr[mcols(feature_gr)$state == state]
    if (length(current_state_gr) == 0) {
      cat("  警告：状态", state, "无区域，跳过\n")
      next
    }

    # 计算重叠（并过滤最小重叠长度）
    overlap_gr <- intersect(target_gr, current_state_gr)
    if (length(overlap_gr) > 0) {
      overlap_gr <- overlap_gr[width(overlap_gr) >= min_overlap_bp]
    }
    total_overlap_bp <- as.numeric(sum(width(overlap_gr)))

    total_state_bp <- as.numeric(sum(width(current_state_gr)))
    total_target_bp_num <- as.numeric(total_target_bp)
    total_background_bp_num <- as.numeric(total_background_bp)

    # 期望重叠（基于 background_gr）
    expected_bp <- total_target_bp_num * total_state_bp / total_background_bp_num
    enrichment_ratio <- ifelse(expected_bp > 0, total_overlap_bp / expected_bp, 0)
    if (!is.finite(enrichment_ratio)) enrichment_ratio <- 1e6

    # 列联表
    a <- total_overlap_bp
    b <- total_target_bp_num - a
    c <- total_state_bp - a
    d <- total_background_bp_num - total_state_bp - b
    c <- max(0, c)
    d <- max(0, d)

    # Odds Ratio
    if (b > 0 && c > 0) {
      odds_ratio <- (a * d) / (b * c)
    } else if (a > 0 && d > 0 && b == 0 && c == 0) {
      odds_ratio <- Inf
    } else {
      odds_ratio <- 1
    }

    # 超几何检验 p 值
    if (total_state_bp > 0 && total_target_bp_num > 0 && total_overlap_bp >= 0) {
      p_val <- phyper(total_overlap_bp - 1,
                      total_state_bp,
                      total_background_bp_num - total_state_bp,
                      total_target_bp_num,
                      lower.tail = FALSE)
      if (is.na(p_val) || p_val < .Machine$double.xmin) {
        p_val <- .Machine$double.xmin
      }
    } else {
      p_val <- NA
    }

    # 保存结果
    res_row <- data.frame(
      tissue_state = target_label,
      State = state,
      Overlap_BP = total_overlap_bp,
      Total_Target_BP = total_target_bp_num,
      Total_State_BP = total_state_bp,
      Total_Background_BP = total_background_bp_num,
      Expected_BP = expected_bp,
      Enrichment_Ratio = enrichment_ratio,
      Odds_Ratio = ifelse(is.infinite(odds_ratio), 1e308,
                          ifelse(is.nan(odds_ratio), 1, odds_ratio)),
      P_value = p_val,
      stringsAsFactors = FALSE
    )
    enrichment_results[[state]] <- res_row

    # 保存重叠区域（如果需要）
    if (return_overlap_gr) {
      overlap_list[[state]] <- overlap_gr
    }
  }

  if (length(enrichment_results) == 0)
    stop("未生成任何富集结果！")

  result_df <- do.call(rbind, enrichment_results)
  rownames(result_df) <- NULL

  # 显著性处理
  if (sig_method == "FDR") {
    result_df$FDR <- p.adjust(result_df$P_value, method = "fdr")
    result_df$Significance <- result_df$FDR
  } else {
    result_df$FDR <- NA_real_
    result_df$Significance <- result_df$P_value
  }

  # 添加 -log10 列（如果需要）
  if (add_log10p) {
    result_df$log10P <- -log10(pmax(result_df$P_value, .Machine$double.xmin))
    result_df$log10FDR <- ifelse(is.na(result_df$FDR),
                                 NA_real_,
                                 -log10(pmax(result_df$FDR, .Machine$double.xmin)))
  }

  # 处理无穷大的 Odds Ratio
  finite_or <- result_df$Odds_Ratio[is.finite(result_df$Odds_Ratio) & result_df$Odds_Ratio != 1e308]
  if (length(finite_or) > 0) {
    max_or <- max(finite_or)
    result_df$Odds_Ratio[result_df$Odds_Ratio == 1e308] <- max_or * 10
  } else {
    result_df$Odds_Ratio[result_df$Odds_Ratio == 1e308] <- 1e6
  }

  # 排序
  if (sort_by == "Enrichment_Ratio") {
    result_df <- result_df[order(-result_df$Enrichment_Ratio), ]
  } else if (sort_by == "Odds_Ratio") {
    result_df <- result_df[order(-result_df$Odds_Ratio), ]
  } else if (sort_by == "P_value") {
    result_df <- result_df[order(result_df$P_value), ]
  } else if (sort_by == "-log10(P_value)") {
    if (!add_log10p) {
      logp <- -log10(pmax(result_df$P_value, .Machine$double.xmin))
      result_df <- result_df[order(-logp), ]
    } else {
      result_df <- result_df[order(-result_df$log10P), ]
    }
  } else if (sort_by == "-log10(FDR)") {
    if (!add_log10p) {
      logfdr <- ifelse(is.na(result_df$FDR),
                       -Inf,
                       -log10(pmax(result_df$FDR, .Machine$double.xmin)))
      result_df <- result_df[order(-logfdr), ]
    } else {
      result_df <- result_df[order(-result_df$log10FDR), ]
    }
  }

  rownames(result_df) <- NULL

  # 返回结果
  if (return_overlap_gr) {
    return(list(result_df = result_df, overlaps = overlap_list))
  } else {
    return(result_df)
  }
}
