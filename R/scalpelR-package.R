#' scalpelR: downstream analysis of SCALPEL (nf-core/scalpel) results
#'
#' @keywords internal
"_PACKAGE"

# data.table syntax in the package functions
.datatable.aware <- TRUE

utils::globalVariables(c(
  ".", ".SD", ".N", ":=", "start", "end", "strand", "gene_name", "transcript_name", "o", "te5", "three_end",
  "significant", "gene_tr", "seqnames", "type", "width", "end_CDS", "UTRlength", "features", "cells", "value", "group", "barcode", "position", "y", "d", "umis", "usage", "id", "N", "label", "feature", "xend", "x"
))
