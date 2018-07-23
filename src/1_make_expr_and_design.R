library(data.table)
library(stringr)

nThread <- 4

local.path <- file.path(getwd(), "data", "Mouse_BMDM_Non_stimulated.csv")
expr.data <- fread(local.path, header = TRUE, sep = "\t", verbose = TRUE, nThread = nThread)

# clean data

expr.data <- expr.data[!is.na(entrezgene_id)]
expr.data <- expr.data[entrezgene_id != ""]
# expr.data$entrezgene_id <- str_replace(expr.data$entrezgene_id, ',.*', '')
# expr.data$entrezgene_id <- str_replace(expr.data$entrezgene_id, 'entrezgene:', '')


# make design matrix

samples <- colnames(expr.data)[c(7:33)]
samples <- str_replace(samples, 'Non-stimulated', 'sample')
times <- str_replace(samples, 'h.*', '')
times <- str_replace(times, '[^0-9]*', '')
reps <- str_sub(samples, -1, -1)
groups <- ifelse(grepl('Mtb.inf', samples), '2', '1') # 1 - not infected, 2 - infected
design <- data.table(sample = samples, rep = reps, time.h = times, group = groups)

# make an expression matrix

colnames(expr.data)[c(7:33)] <- samples

expr.data[, gene_descr := str_remove(short_description, '.*@')]
no_peaks <- expr.data[, .SD[, entrezgene_id[1]], by=gene_descr]

sum.expr.data <- expr.data[, lapply(.SD, sum, na.rm=TRUE), by=gene_descr, .SDcols=samples ]
expr.data <- merge(no_peaks, sum.expr.data, by = 'gene_descr')
rm(no_peaks, sum.expr.data)
colnames(expr.data)[2] <- 'entrezgene_id'

expr.data[, rep:= as.integer(.SD[, length(strsplit(entrezgene_id, ',')[[1]])]),
                               entrezgene_id]
with.duplicates <- expr.data[,cbind(.SD,dup=1:rep),by='entrezgene_id']

single.entrez <- expr.data[, .SD[, strsplit(entrezgene_id, ',')], entrezgene_id]
colnames(single.entrez) <- c('entrezgene_id', 'entrez_single')

identical(with.duplicates$entrezgene_id, single.entrez$entrezgene_id)
with.duplicates[, entrez_single := single.entrez$entrez_single]
 
# dirty hack:
with.duplicates = with.duplicates[!duplicated(with.duplicates$entrez_single),]

expr.data <- as.matrix(with.duplicates[, samples, with = F])
rownames(expr.data) <- str_remove(with.duplicates$entrez_single, 'entrezgene:')

design.df <- as.data.frame(design)
as.numeric(design.df[, 'time.h']) -> design.df[, 'time.h']
as.factor(design.df[, 'rep']) -> design.df[, 'rep']
as.factor(design.df[, 'group']) -> design.df[, 'group']

identical(colnames(expr.data), design.df$sample)
# save everything

save(expr.data, design.df, file = 'out/expr_and_design.RData')
