require(data.table)
require(stringr)

# loading list of significant genes
sign.genes <- readRDS(file = 'out/sign_genes_model_1.RDS')
n.sign <- length(sign.genes)

# load interesting genes
genes.of.interest <- readRDS(file = 'out/genes_in_complexes.RDS')
all.genes <- paste0('X', genes.of.interest)
n.genes <- length(all.genes)

# load gmt object
load('out/gmt.RData')

fisher.res.dt <- rbindlist(lapply(gmt.complexes$genesets, function(geneset){
  geneset.genes <- paste0('X', geneset)
  n.geneset <- length(intersect(all.genes, geneset.genes))
  n.geneset.sign <- length(intersect(sign.genes, geneset.genes))
  m = matrix(c(n.geneset.sign, n.sign, n.geneset, n.genes), ncol=2)
  res <- fisher.test(m)
  data.table(p.value = res$p.value)
}))

fisher.res.dt[, geneset := gmt.complexes$geneset.descriptions]
fisher.res.dt[, p.adj := p.adjust(p.value)]
setcolorder(fisher.res.dt, c("geneset", "p.value", "p.adj"))

path <- file.path(getwd(), 'out', 'fisher_model_1.csv')
write.table(fisher.res.dt, file = path, append = F, 
            quote = F, sep = '\t', col.names = T, row.names = F)
