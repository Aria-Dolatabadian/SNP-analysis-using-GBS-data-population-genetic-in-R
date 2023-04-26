
#https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html


library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(igraph)
library(ggplot2)
library(tidyr)
library(phangorn)

#Opening and examining the vcf file

rubi.VCF <- read.vcfR("prubi_gbs.vcf.gz")

pop.data <- read.table("population_data.gbs.txt", sep = "\t", header = TRUE)

all(colnames(rubi.VCF@gt)[-1] == pop.data$AccessID)

#Converting the dataset to a genlight object

gl.rubi <- vcfR2genlight(rubi.VCF)

ploidy(gl.rubi) <- 2

pop(gl.rubi) <- pop.data$State

gl.rubi

#Population genetic analyses for GBS data
#Distance matrices

x.dist <- dist(x)

x.dist <- poppr::bitwise.dist(gl.rubi)

#Distance tree

tree <- aboot(gl.rubi, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

cols <- brewer.pal(n = nPop(gl.rubi), name = "Dark2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.rubi)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
legend('topleft', legend = c("CA","OR","WA"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")

#Minimum spanning networks

library(igraph)

rubi.dist <- bitwise.dist(gl.rubi)
rubi.msn <- poppr.msn(gl.rubi, rubi.dist, showplot = FALSE, include.ties = T)

node.size <- rep(2, times = nInd(gl.rubi))
names(node.size) <- indNames(gl.rubi)
vertex.attributes(rubi.msn$graph)$size <- node.size

set.seed(9)
plot_poppr_msn(gl.rubi, rubi.msn , palette = brewer.pal(n = nPop(gl.rubi), name = "Dark2"), gadj = 70)

#Principal components analysis

rubi.pca <- glPca(gl.rubi, nf = 3)
barplot(100*rubi.pca$eig/sum(rubi.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)


rubi.pca.scores <- as.data.frame(rubi.pca$scores)
rubi.pca.scores$pop <- pop(gl.rubi)

library(ggplot2)
set.seed(9)
p <- ggplot(rubi.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()

p


#DAPC

pnw.dapc <- dapc(gl.rubi, n.pca = 3, n.da = 2)

scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)

#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(pnw.dapc,col = cols, posi = 'top')

dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.rubi)
dapc.results$indNames <- rownames(dapc.results)

# library(reshape2)
# dapc.results <- melt(dapc.results)
library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

head(dapc.results, n = 6)

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")


p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

#Subsetting a vcfR object to 200 random variants

alphabet <- c("A","B","C","D")
alphabet[c(1,2)]

rubi.VCF

rubi.VCF[c(1:10),]

subset.1 <- sample(size = 200, x= c(1:nrow(rubi.VCF)))
subset.2 <- sample(size = 200, x= c(1:nrow(rubi.VCF)))

identical(subset.1, subset.2)

rubi.VCF.sub1 <- rubi.VCF[subset.1,]
rubi.VCF.sub2 <- rubi.VCF[subset.2,]

rubi.VCF.sub1

rubi.VCF.sub2

# Creating a list object to save our subsets in.
rubi.variant.subset <- vector(mode = "list", length = 50)

# Using a for loop to generate 50 subsets of 200 random variants from the rubi.VCF vcfR object.
for (i in 1:50){
  rubi.variant.subset[[i]] <- rubi.VCF[sample(size = 200, x= c(1:nrow(rubi.VCF)))]
}

# Checking we have 50 vcfR objects:
length(rubi.variant.subset)

head(rubi.variant.subset, n=2)

# Creating the GenLight object
rubi.gl.subset <- lapply(rubi.variant.subset, function (x) suppressWarnings(vcfR2genlight(x)))
for (i in 1:length(rubi.gl.subset)){
  ploidy(rubi.gl.subset[[i]]) <- 2
}

# Creating a simple UPGMA tree per object
library(phangorn)
rubi.trees <- lapply(rubi.gl.subset, function (x) upgma(bitwise.dist(x)))
class(rubi.trees) <- "multiPhylo"

# Overlapping the trees
densiTree(rubi.trees, consensus = tree, scaleX = T, show.tip.label = F, alpha = 0.1)
title(xlab = "Proportion of variants different")



