#- simulated iDGE Seurat object: 2 samples (one counts layer each), 2 cell types, 3 genes
sim_obj <- function() {
  set.seed(1)
  n <- 100
  feats <- c("Ga***Ga-201", "Ga***Ga-202", "Gb***Gb-201", "Gb***Gb-202", "Gc***Gc-201")
  one <- function(s) {
    ct <- rep(c("A", "B"), each = n / 2)
    m <- rbind(
      rpois(n, ifelse(ct == "A", 16, 4)), rpois(n, ifelse(ct == "A", 4, 16)),
      rpois(n, 10), rpois(n, 10), rpois(n, 5))
    dimnames(m) <- list(feats, paste0("C", seq_len(n), "-1"))
    o <- SeuratObject::CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE), project = s)
    o$orig.ident <- s
    o$Barcodes <- colnames(o)
    o$celltype <- ct
    o
  }
  merge(one("S1"), one("S2"), add.cell.ids = c("S1", "S2"))
}
annots <- function(dir) {
  a <- data.table::data.table(
    seqnames = "chr1", start = c(1000, 1000, 5000, 9000, 20000), end = c(1100, 2100, 5200, 9300, 20100), width = 0,
    strand = "+", startR = 0, endR = 0, gene_id = "g", gene_name = c("Ga", "Ga", "Gb", "Gb", "Gc"),
    transcript_id = "t", transcript_name = c("Ga-201", "Ga-202", "Gb-201", "Gb-202", "Gc-201"), exon_number = 1,
    TPM_perc = 1, transcript_collapseds = "", n_isoforms_in_profile = 1)
  f <- file.path(dir, "chr1.annots.gz")
  data.table::fwrite(a, f, sep = "\t")
  f
}

test_that("FindIsoforms on a Seurat object with layers = pipeline engine", {
  obj <- sim_obj()
  f <- annots(tempdir())
  res <- FindIsoforms(obj, group_by = "celltype", isoform_annotation = f, extended = TRUE)
  m <- SeuratObject::LayerData(SeuratObject::JoinLayers(obj), layer = "counts")
  ref <- scalpelR:::.diu_engine(m, obj$celltype, strata = obj$orig.ident, ends = scalpelR:::.Read_isoformEnds(f))
  ref[, significant := NULL]
  expect_equal(res, ref)
  expect_identical(colnames(res)[1:2], c("gene_name", "transcript_name"))
  ga <- res[gene_name == "Ga"]
  expect_true(all(ga$event_type == "tandem_3UTR") && all(ga$direction == "shortening") && isTRUE(ga$replicated[1]))
  expect_true(all(is.na(res[gene_name == "Gb"]$direction)))
  #- default: genes with a differential isoform usage only (Ga)...
  expect_identical(unique(FindIsoforms(obj, group_by = "celltype", isoform_annotation = f)$gene_name), "Ga")
})

test_that("FindIsoforms permutation test and GTF annotation", {
  obj <- sim_obj()
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1000, 1000, 5000, 9000, 20000), c(1100, 2100, 5200, 9300, 20100)),
    strand = "+", type = "exon", gene_name = c("Ga", "Ga", "Gb", "Gb", "Gc"), transcript_name = c("Ga-201", "Ga-202", "Gb-201", "Gb-202", "Gc-201"))
  res <- FindIsoforms(obj, group_by = "celltype", isoform_annotation = gr, test = "permutation", n_perm = 500, extended = TRUE)
  expect_true(all(res$test == "permutation"))
  expect_true(res[gene_name == "Ga", p.value[1] <= 1 / 50 & event_type[1] == "tandem_3UTR"])
  expect_true(res[gene_name == "Gb", p.value[1] > 0.01])
})

test_that("weightedUTR_cell: sum over isoforms of usage x 3' UTR length, sparse, layers joined", {
  obj <- sim_obj()
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(c(100, 201, 100, 201, 5000, 5101, 9000, 9101), c(200, 300, 200, 1300, 5100, 5300, 9100, 9200)), strand = "+",
    type = rep(c("CDS", "UTR"), 4), gene_name = rep(c("Ga", "Gb"), each = 4), transcript_name = rep(c("Ga-201", "Ga-202", "Gb-201", "Gb-202"), each = 2))
  w <- weightedUTR_cell(obj, gr, group.by = "celltype", min.counts = 0)
  m <- as.matrix(SeuratObject::LayerData(SeuratObject::JoinLayers(obj), layer = "counts"))
  len <- c(100, 1100, 200, 100)
  ga <- (m[1, ] * len[1] + m[2, ] * len[2]) / (m[1, ] + m[2, ])
  gb <- (m[3, ] * len[3] + m[4, ] * len[4]) / (m[3, ] + m[4, ])
  expected <- rowMeans(cbind(ga, gb), na.rm = TRUE)
  expect_equal(unname(w$wUTRlength.cell), unname(expected))
  expect_identical(w$cells, colnames(obj))
})

test_that("plotExpression on a Seurat object with layers", {
  obj <- sim_obj()
  p <- plotExpression(obj, features = c("Ga***Ga-201", "Ga***Ga-202"), group_var = "celltype")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 2 * ncol(obj))
})

test_that("plotUsage: usage per group (and sample), isoforms ordered by 3' end", {
  obj <- sim_obj()
  f <- annots(tempdir())
  p <- plotUsage(obj, gene_in = "Ga", group_by = "celltype", sample_by = "orig.ident", isoform_annotation = f)
  d <- data.table::as.data.table(p$data)
  expect_identical(levels(d$transcript_name), c("Ga-201", "Ga-202"))
  expect_true(all(abs(d[, sum(usage), by = .(sample, group)]$V1 - 1) < 1e-9))
  m <- SeuratObject::LayerData(SeuratObject::JoinLayers(obj), layer = "counts")
  a <- obj$celltype == "A" & obj$orig.ident == "S1"
  expect_equal(d[sample == "S1" & group == "A" & transcript_name == "Ga-201"]$usage, sum(m["Ga***Ga-201", a]) / sum(m[c("Ga***Ga-201", "Ga***Ga-202"), a]))
})

test_that("CoveragePlot: reads of each group from the BAM files, per 1000 cells, strand of the gene", {
  skip_if(Sys.which("samtools") == "", "samtools not available")
  obj <- sim_obj()
  dir <- tempfile()
  dir.create(dir)
  #- 2 reads of cell C1 (type A, S1) at 1001-1050 (+), 1 read of cell C60 (type B, S1) at 1021-1070 with a 10 nt deletion, 1 read on - strand...
  sam <- c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:30000",
    "r1\t0\tchr1\t1001\t255\t50M\t*\t0\t0\t*\t*\tCB:Z:C1-1",
    "r2\t0\tchr1\t1001\t255\t50M\t*\t0\t0\t*\t*\tCB:Z:C1-1",
    "r4\t16\tchr1\t1001\t255\t50M\t*\t0\t0\t*\t*\tCB:Z:C2-1",
    "r3\t0\tchr1\t1021\t255\t20M10D30M\t*\t0\t0\t*\t*\tCB:Z:C60-1")
  writeLines(sam, file.path(dir, "s1.sam"))
  expect_equal(system(paste("samtools view -b -o", file.path(dir, "s1.bam"), file.path(dir, "s1.sam"), "&& samtools index", file.path(dir, "s1.bam"))), 0)
  #- Ga-201 coding, Ga-202 non-coding (no CDS)...
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(900, 1000, 900), c(2100, 1100, 1500)), strand = "+", type = c("exon", "CDS", "exon"),
    gene_name = "Ga", transcript_name = c("Ga-201", "Ga-201", "Ga-202"))
  p <- CoveragePlot(obj, bamfiles = c(S1 = file.path(dir, "s1.bam")), gene_in = "Ga", genome_gr = gr, group_by = "celltype", features = c("Ga-201", "Ga-202"))
  expect_setequal(levels(data.table::as.data.table(p[[2]]$layers[[3]]$data)$transcript_name), c("Ga-201", "Ga-202"))
  d <- data.table::as.data.table(p[[1]]$data)
  n <- table(obj$celltype)
  expect_equal(d[grepl("^A", group) & position == 1010]$depth, 2 / n[["A"]] * 1000)
  expect_equal(d[grepl("^B", group) & position == 1060]$depth, 1 / n[["B"]] * 1000)
  expect_equal(d[grepl("^A", group) & position == 1060]$depth, 0)
  expect_equal(d[grepl("^B", group) & position == 1045]$depth, 1 / n[["B"]] * 1000)
})

test_that("plotShiftHeatmap: relative 3' UTR length per group & sample, and along a trajectory", {
  obj <- sim_obj()
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(c(100, 201, 100, 201, 5000, 5101, 9000, 9101), c(200, 300, 200, 1300, 5100, 5300, 9100, 9200)), strand = "+",
    type = rep(c("CDS", "UTR"), 4), gene_name = rep(c("Ga", "Gb"), each = 4), transcript_name = rep(c("Ga-201", "Ga-202", "Gb-201", "Gb-202"), each = 2))
  p <- plotShiftHeatmap(obj, gr, group_by = "celltype", sample_by = "orig.ident", scale = FALSE)
  d <- attr(p, "shift")
  expect_equal(nrow(d), 2 * 4)
  #- Ga: relative 3' UTR length = (weighted length - 100) / (1100 - 100), shorter in A...
  m <- as.matrix(SeuratObject::LayerData(SeuratObject::JoinLayers(obj), layer = "counts"))
  a <- obj$celltype == "A" & obj$orig.ident == "S1"
  rel <- (sum(m[1, a]) * 100 + sum(m[2, a]) * 1100) / sum(m[1:2, a])
  expect_equal(d[gene_name == "Ga" & group == "A" & sample == "S1"]$rel_3UTR, (rel - 100) / 1000)
  expect_true(all(d[gene_name == "Ga" & group == "A"]$value < 0) && all(d[gene_name == "Ga" & group == "B"]$value > 0))
  expect_true(all(abs(d[, sum(value), by = gene_name]$V1) < 1e-9))
  #- trajectory: A cells first, Ga 3' UTR lengthening along the bins...
  obj$pt <- ifelse(obj$celltype == "A", 0, 1) + stats::runif(ncol(obj))
  p2 <- plotShiftHeatmap(obj, gr, group_by = "celltype", order_by = "pt", n_bins = 4, min_umis = 1)
  d2 <- attr(p2, "shift")[gene_name == "Ga"][order(column)]
  expect_s3_class(p2, "patchwork")
  expect_true(d2$value[1] < 0 && d2$value[4] > 0)
})

test_that("plotShiftHeatmap: clusters of genes", {
  obj <- sim_obj()
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(c(100, 201, 100, 201, 5000, 5101, 9000, 9101), c(200, 300, 200, 1300, 5100, 5300, 9100, 9200)), strand = "+",
    type = rep(c("CDS", "UTR"), 4), gene_name = rep(c("Ga", "Gb"), each = 4), transcript_name = rep(c("Ga-201", "Ga-202", "Gb-201", "Gb-202"), each = 2))
  d <- attr(plotShiftHeatmap(obj, gr, group_by = "celltype", n_clusters = 2), "shift")
  expect_identical(nlevels(d$cluster), 2L)
  expect_identical(data.table::uniqueN(d[, .(gene_name, cluster)]), 2L)
})
