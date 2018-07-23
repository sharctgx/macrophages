library(TcGSA)
library(data.table)
library(stringr)
library(gridExtra)


# load expression & design matrices
load('out/expr_and_design.RData')

# load gmt object
load('out/gmt.RData')

# perform TcGSA

GSE.gene.func <- TcGSA.LR(expr = expr.data, gmt = gmt.genes, design = design.df,
                          subject_name = 'rep', time_name = 'time.h', group_name = 'group',
                          minGSsize = 5)

# expr.data.div <- expr.data/1000
# design.df$time.h.div <- design.df$time.h/10
design.df$time.ranked <- as.numeric(factor(design.df$time.h))

GSE.complexes <- TcGSA.LR(expr = expr.data, gmt = gmt.complexes, design = design.df,
                          subject_name = 'rep', time_name = 'time.ranked', 
                          time_func = "linear", group_name = 'group', minGSsize = 5)

sgnifs <- signifLRT.TcGSA(GSE.complexes, threshold = 0.05, myproc = "BH",
                          nbsimu_pval = 1000, write=FALSE)

expr.df <- as.data.frame(expr.data)
expr.df <- sapply(expr.df, as.numeric)

summary(GSE.complexes)
plotFit.GS(GSE.complexes, expr.df, design.df, subject_name = "rep",
           time_name = 'time.h', colnames_ID = 'sample', plot_type = "Histogram Obs", 
           GeneSetsList = sgnifs$mixedLRTadjRes$GeneSet,
           color = c("genes", "time", "subjects"), marginal_hist = TRUE,
           gg.add = list(theme()))

# trying some DTW

require(dtwclust)

# checking on uninfected sample replicate 1

cols <- design.df[(design.df$rep == 1) & (design.df$group == 1),]$sample
rows <- sample(dim(expr.data)[1], size = 1000)

Clust.uninf.rep1 <- tsclust(series = expr.data[, cols], type = 'partitional', 
                control = partitional_control(pam.precompute = F),
                distance = "dtw_lb", window.size = 1, trace = T, error.check = T)
save(Clust.uninf.rep1 , file = 'clusters_uninf_rep1.RData')

centr.l <- Clust.uninf.rep1@centroids
centroids.m <- do.call(rbind, centr.l)
plot(centroids.m[1, ], ylim = c(min(centroids.m), max(centroids.m)), type="n",
     xlab = 'time point', ylab = 'Expression', main = 'Centroids (Shape-based distance,
     partitional clustering, no infection, replicate 1)')
mapply(lines, centr.l,col=seq_along(centr.l),lty=2)
legend("topleft",legend = 1:length(centr.l),lty=2,
       col=seq_along(centr.l))

expr.data.shape <- t(apply(expr.data[, cols], 1, function(row)
  row/sum(row)))

load(file = 'clusters_uninf_rep1.RData')







