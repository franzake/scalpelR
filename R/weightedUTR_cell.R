#' Calculate Weighted UTR Length per Cell
#'
#' This function computes the weighted UTR length per cell using transcript information from a Seurat object.
#' It requires a gene annotation file with UTR and CDS information to calculate UTR lengths across isoforms.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param annotation_gr A GRanges object containing gene annotation, including UTR and CDS regions.
#' @param group.by The metadata column in the Seurat object to group cells by.
#' @param features A vector of features (genes/transcripts) to subset the Seurat object by. Default is NULL.
#' @param assay The assay to use for expression values, default is "RNA".
#' @param min.counts The minimum number of counts required to consider a transcript for further analysis. Default is 3.
#'
#' @return A data.table with one row per cell: \code{cells}, \code{group}, \code{wUTRlength.cell} (mean over the genes
#'   expressed in the cell of the gene weighted 3' UTR length) and \code{n_genes} (genes used).
#'
#' @details
#' For each gene and cell, the weighted 3' UTR length is the sum over its isoforms of the isoform fraction
#' of the gene counts x the isoform 3' UTR length (isoforms with more than \code{min.counts} counts in the
#' cell). Counts split in several layers (one per sample) are joined; computed on sparse matrices.
#'
#' @examples
#' \dontrun{
#' weightedUTR_cell(seurat_obj = my_seurat_obj, features = NULL,
#'                  annotation_gr = my_annotation, group.by = "Cluster",
#'                  assay = "RNA", min.counts = 3)
#' }
#'
#' @export
weightedUTR_cell <- function(
        seurat_obj,
        annotation_gr,
        group.by,
        features = NULL,
        assay = "RNA",
        min.counts = 3) {
    # Checks
    if (is.null(group.by)) stop("Error! Provide group.by argument.")
    if (!(group.by %in% colnames(seurat_obj[[]]))) stop(paste0("Error! ", group.by, " not found in Seurat object."))

    message("Processing annotation file...(1/3)")
    #- 3' UTR length of each transcript: UTR segments downstream of the CDS end...
    my_annot <- .utr3_length(annotation_gr)

    message("Processing isoform expression...(2/3)")
    #- isoform counts (sparse, layers joined) of the annotated isoforms...
    seurat_obj <- SeuratObject::JoinLayers(seurat_obj, assay = assay)
    counts <- SeuratObject::LayerData(seurat_obj, assay = assay, layer = "counts")
    feats <- data.table::data.table(gene_tr = rownames(counts))
    feats[, c("gene_name", "transcript_name") := data.table::tstrsplit(gene_tr, "***", fixed = TRUE)]
    feats <- merge(feats, my_annot, by = c("gene_name", "transcript_name"), sort = FALSE)
    if (!is.null(features)) feats <- feats[gene_tr %in% features | gene_name %in% features]
    if (nrow(feats) == 0) stop("No isoform of the Seurat object found in the annotation (UTR / CDS).")
    x <- counts[feats$gene_tr, , drop = FALSE]
    x@x[x@x <= min.counts] <- 0
    x <- Matrix::drop0(x)

    #- gene weighted 3' UTR length per cell = sum(isoform fraction x isoform 3' UTR length)...
    genes <- factor(feats$gene_name)
    A <- Matrix::sparseMatrix(i = as.integer(genes), j = seq_along(genes), x = 1, dims = c(nlevels(genes), length(genes)))
    tot <- A %*% x
    wlen <- A %*% (x * feats$UTRlength)
    expressed <- tot > 0
    tot_inv <- tot
    tot_inv@x <- 1 / tot_inv@x
    wgene <- wlen * tot_inv

    message("Integrating Seurat metadata...(3/3)")
    n_genes <- Matrix::colSums(expressed)
    my_expr <- data.table::data.table(
        cells = colnames(seurat_obj),
        group = seurat_obj[[]][, group.by],
        wUTRlength.cell = ifelse(n_genes > 0, Matrix::colSums(wgene) / n_genes, NA_real_),
        n_genes = n_genes)
    return(my_expr)
}
