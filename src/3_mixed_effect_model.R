require(data.table)
require(lme4)
require(stringr)
require(optimx)

# load expression & design matrices
load('out/expr_and_design.RData')

# load gmt object
load('out/gmt.RData')

draw.genes <- function(){
  lapply(sign.genes, function(gene){
    # no infection
    lapply(unique(expr.sign[expr.sign$group == '1', ]$rep), function(rep){
      dt <- expr.sign[(expr.sign$group == '1') & (expr.sign$rep == rep), ]
      setorder(dt, time.ranked)
      lines(x = dt$time.ranked, y = dt[, gene], col = 'sienna2')
    })
    
    
    # infected
    lapply(unique(expr.sign[expr.sign$group == '2', ]$rep), function(rep){
      dt <- expr.sign[(expr.sign$group == '2') & (expr.sign$rep == rep), ]
      setorder(dt, time.ranked)
      lines(x = dt$time.ranked, y = dt[, gene], col = 'slateblue4')
    })
  })
  legend("topleft",legend = c('control', 'Mtb infected'), lty=1, 
         col=c('sienna2', 'slateblue4'))
}

# interesting genes for complexes:
genes.of.interest <- unique(unlist(gmt.complexes$genesets))
genes.of.interest <- genes.of.interest[genes.of.interest %in% rownames(expr.data)]

expr.small <- as.data.frame(t(expr.data[genes.of.interest, ]))
identical(rownames(expr.small), design.df$sample)
colnames(expr.small) <- paste0('X', colnames(expr.small))
genes.of.interest <- paste0('X', genes.of.interest)
saveRDS(genes.of.interest, file = 'out/genes_in_complexes.RDS')

design.df$time.ranked <- as.numeric(factor(design.df$time.h))
# write.csv(file='out/design.csv', x=design.df)
expr.small <- cbind(expr.small, design.df)
meta.cols <- colnames(design.df)

model.res.dt <- rbindlist(lapply(genes.of.interest, function(gene){
    f1 <- as.formula(paste(gene, "~ group + (1|time.ranked)"))
    f2 <- as.formula(paste(gene, "~ (1|time.ranked)"))

    model.gene = lmer(f1, data = expr.small, REML = F)
    model.gene.null = lmer(f2, data = expr.small, REML = F)

    p.val.str <-
      as.character(as.data.frame(summary(anova(model.gene, model.gene.null)))$Freq[[55]])
    p.value <- as.numeric(str_remove(p.val.str, '[^0-9]+'))

    stdev <- c(as.data.frame(VarCorr(model.gene))$sdcor[1])
    t.score <- as.data.frame(coef(summary(model.gene))[,"t value"])[2, ]
    data.table(gene = gene, stdev = stdev, t.score = t.score, anova.p.value = p.value)
}))

model.res.dt$p.adj <- p.adjust(model.res.dt$anova.p.value)
path <- file.path(getwd(), 'out', 'model_random_intercept.csv')
write.table(model.res.dt, file = path, append = F, 
            quote = F, sep = '\t', col.names = T, row.names = F)

sign.genes <- model.res.dt[model.res.dt$p.adj < 0.05]$gene
saveRDS(sign.genes, file = 'out/sign_genes_model_1.RDS')

expr.sign <- expr.small[, c(sign.genes, meta.cols)]

plot(0, 0, xlim = c(1, 5), ylim = c(min(expr.sign[, sign.genes]), 
                                    max(expr.sign[, sign.genes])
                                    ), type="n",
     xlab = 'time point', ylab = 'Expression', main = 'Significant genes')

draw.genes()

# drawing mean between the replics (a line for each group) for 4 random genes:

expr.small.dt <- as.data.table(expr.small)

par(mfrow = c(2, 2), mar = c(5, 5, 2, 2))
lapply(sample(sign.genes, 4), function(gene){
  dt <- expr.small.dt[, .SD[, mean(get(gene))], by = c('group', 'time.ranked')]
  plot(dt[group == 1]$time.ranked, dt[group==1]$V1, xlim =c(1, 5), 
       ylim = c(min(dt$V1), max(dt$V1)), type = 'l', xlab = 'time point',
       ylab = paste(gene, 'expression'), col = 'sienna2')
  lines(dt[group == 2]$time.ranked, dt[group==2]$V1, xlim =c(1, 5), 
        ylim = c(min(dt$V1), max(dt$V1)), col = 'slateblue4')
})

# 












