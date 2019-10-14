library("DESeq2")
library("dplyr")
library("tibble")
library("BiocParallel")
register(MulticoreParam(as.integer(snakemake@threads)))

# Meta data
colData <- read.table(snakemake@input[["sampletable"]], stringsAsFactors = FALSE, header = TRUE) %>%
    select(-orig_filenames) %>%
    filter(tissue == snakemake@wildcards[["tissue"]]) %>%
    filter(driver == snakemake@wildcards[["driver"]]) %>%
    mutate(target = as.factor(.$target)) %>%
    mutate(replicate = as.factor(.$replicate)) %>%
    column_to_rownames("samplename")

colData$target <- relevel(colData$target, "dam")

# Count data
cnts <- read.table(snakemake@input[["counts"]], stringsAsFactors = FALSE, header = TRUE) %>%
    rename(id = 1) %>%
    select(id, row.names(colData)) %>%
    filter(rowSums(select(., -id)) > 0) %>%
    column_to_rownames("id")

# Fit model
dds <- DESeqDataSetFromMatrix(countData = cnts, colData = colData, design = ~target)
dds <- DESeq(dds, fitType = "local", parallel = TRUE)
res <- results(dds, alpha = 0.05, parallel = TRUE)
res <- res %>% as.data.frame %>% rownames_to_column("GATC_site")

# Save results
write.table(res, file = snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE)