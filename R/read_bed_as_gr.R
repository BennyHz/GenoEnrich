#' Read a 3-column BED file as a GRanges object
#'
#' Reads a BED file (chrom, start, end) and converts it to a GRanges object.
#'
#' @param file Path to the BED file (must have at least 3 columns)
#' @return A GRanges object
#' @export
#' @importFrom GenomicRanges makeGRangesFromDataFrame
read_bed_as_gr <- function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
  colnames(df) <- c("seqnames", "start", "end")
  GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
}
