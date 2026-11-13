#' Calculate Weighted UTR Length per Cell
#'
#' This function computes the weighted 3' UTR length per cell using the isoform quantification of a Seurat object
#' and a gene annotation with UTR and CDS information (3' UTR length of each isoform).
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param annotation_gr A GRanges object containing gene annotation, including UTR and CDS regions.
#' @param group.by The metadata column in the Seurat object to group cells by.
#' @param features A vector of features (genes/transcripts) to subset the Seurat object by. Default is NULL.
#' @param assay The assay to use for expression values, default is "RNA".
#' @param method \code{"centered"} (default) or \code{"length"}, see Details.
#' @param min.counts Isoforms with \code{min.counts} counts or less in a cell are not used in this cell. Default:
#'   0 with \code{method = "centered"}, 3 with \code{method = "length"}.
#'
#' @return A data.table with one row per cell: \code{cells}, \code{group}, \code{wUTRlength.cell} and
#'   \code{n_genes} (genes used in the cell).
#'
#' @details
#' For each gene and cell, the weighted 3' UTR length is the sum over its isoforms of the isoform fraction of the gene
#' counts x the isoform 3' UTR length.
#'
#' \code{method = "centered"}: only genes whose isoforms differ in 3' UTR length are used, and the weighted 3' UTR length
#' of each gene in the cell is centered on the weighted 3' UTR length of the gene in all the cells of the object;
#' \code{wUTRlength.cell} is the mean of these differences over the genes of the cell (nt; > 0: longer 3' UTRs than the
#' average cell, < 0: shorter). It does not depend on which genes are detected in the cell, i.e. on the sequencing depth.
#'
#' \code{method = "length"}: \code{wUTRlength.cell} is the mean over the genes of the cell of their weighted 3' UTR
#' length (nt). Genes with a single 3' UTR length contribute their fixed length: the value depends on the genes detected
#' in the cell, which increase with the sequencing depth.
#'
#' Counts split in several layers (one per sample) are joined; computed on sparse matrices.
#'
#' @examples
#' \dontrun{
#' weightedUTR_cell(seurat_obj = my_seurat_obj, annotation_gr = my_annotation, group.by = "Cluster")
#' }
#'
#' @export
weightedUTR_cell <- function(
        seurat_obj,
        annotation_gr,
        group.by,
        features = NULL,
        assay = "RNA",
        method = "centered",
        min.counts = NULL) {
    # Checks
    if (is.null(group.by)) stop("Error! Provide group.by argument.")
    if (!(group.by %in% colnames(seurat_obj[[]]))) stop(paste0("Error! ", group.by, " not found in Seurat object."))
    if (!(method %in% c("centered", "length"))) stop("Error! method must be 'centered' or 'length'.")
    if (is.null(min.counts)) min.counts <- if (method == "centered") 0 else 3

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
    #- centered: genes whose isoforms differ in 3' UTR length only...
    if (method == "centered") feats <- feats[, if (data.table::uniqueN(UTRlength) > 1) .SD, by = gene_name]
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
    if (method == "centered") {
        #- difference to the weighted 3' UTR length of the gene in all the cells...
        ref <- as.vector(Matrix::rowSums(wlen) / pmax(Matrix::rowSums(tot), 1e-12))
        wgene <- wgene - Matrix::Diagonal(x = ref) %*% (expressed * 1)
    }

    message("Integrating Seurat metadata...(3/3)")
    n_genes <- Matrix::colSums(expressed)
    my_expr <- data.table::data.table(
        cells = colnames(seurat_obj),
        group = seurat_obj[[]][, group.by],
        wUTRlength.cell = ifelse(n_genes > 0, Matrix::colSums(wgene) / n_genes, NA_real_),
        n_genes = n_genes)
    return(my_expr)
}
