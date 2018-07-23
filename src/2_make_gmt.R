library(GSA)
library(org.Mm.eg.db)

nThread <- 4

# Complexes
path1 <- file.path(getwd(), 'data', 'protein_complexes.csv')
complexes <- fread(path1, sep = "\t", verbose = TRUE, nThread = nThread)

complex.groups <- complexes[, c("Complex group name", "Genes in complex"), with = F]
colnames(complex.groups) <- c('complex_group', 'genes')

complex.groups.u <- complex.groups[, .SD[, paste(genes, sep="", 
                            collapse=", ")], by = complex_group]

colnames(complex.groups.u) <- colnames(complex.groups)
# complex.groups.u[, rep:= as.integer(.SD[, length(strsplit(genes, c(', |\\|'))[[1]])]), genes]
complex.groups.u <- complex.groups.u[, .SD[, strsplit(genes, ', |\\|')],
                                     by = c('genes', 'complex_group')]
colnames(complex.groups.u) <- c('genes', 'complex_group', 'single_gene')

complex.groups.u$single_gene <- str_remove(complex.groups.u$single_gene, '\\(|\\)|\\+')
# if we doubt it -- we don't use it:
complex.groups.u <- complex.groups.u[!single_gene %like% "\\?"]

simp.converter <- function(x) {
  s1 <- str_sub(x, 1, 1)
  s2 <- tolower(str_sub(x, 2))
  paste0(s1, s2)
}


# convert Symbol to Entrez
mm <- org.Mm.egSYMBOL2EG
my.symbols <- complex.groups.u$single_gene
my.symbols.low <- unlist(lapply(my.symbols, function(s) simp.converter(s)))
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(mm)
# Convert to a list
annot.big <- as.list(mm[mapped_genes])
annot <- lapply(my.symbols.low, function(x){
  if (x %in% names(annot.big)){
    paste(annot.big[[x]])}
  else{
    paste(NA)
  }
})


#identical(complex.groups.u$single_gene, annot$SYMBOL)
complex.groups.u$single_gene <- annot

complex.groups.u <- complex.groups.u[!(single_gene == 'NA')]
complex.groups <- complex.groups.u[, .SD[, paste(single_gene, collapse = '\t')],
                                   by = complex_group]


path4 <- file.path(getwd(), 'out', 'complexes_GSE.txt')
write.table(complex.groups, 
            file = path4, append = F, 
            quote = F, sep = '\t', col.names = F, row.names = T)

#######################################################################
#Genes
path2 <- file.path(getwd(), 'data', 'genes.csv')
genes <- fread(path2, sep = "\t", verbose = TRUE, nThread = nThread)

gene.func <- genes[, c("Entrez gene ID", "Function"), with = F]
colnames(gene.func) <- c('entrezgene_id', 'func')
gene.func <- gene.func[entrezgene_id != '#' & func != '#']

formatted.genes <- gene.func[, .SD[, paste(entrezgene_id, 
                                           sep="", collapse="\t")], by = func]
colnames(formatted.genes) <- c('func', 'genes')
formatted.genes[, rep:= as.integer(.SD[, length(strsplit(func, ', ')[[1]])]), func]
duplicated.func <- formatted.genes[,cbind(.SD,dup=1:rep), by='func']

single.func <- formatted.genes[, .SD[, strsplit(func, ', ')], func]
colnames(single.func) <- c('func', 'single_func')

identical(duplicated.func$func, single.func$func)
duplicated.func[, single_func := single.func$single_func]
formatted.genes <- duplicated.func[, .SD[, paste(unique( strsplit(genes, '\t')[[1]]),
                collapse = '\t')], by = single_func]

path3 <- file.path(getwd(), 'out', 'genes_GSE.txt')
write.table(formatted.genes, 
            file = path3, append = F, 
            quote = F, sep = '\t', col.names = F, row.names = T)

######################################################################

# load the files to make a .gmt object
gmt.genes <- GSA.read.gmt(path3)
gmt.complexes <- GSA.read.gmt(path4)

save(gmt.genes, gmt.complexes, file = 'out/gmt.RData')
