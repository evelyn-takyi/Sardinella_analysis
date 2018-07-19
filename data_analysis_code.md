---
title: "adegenet1"
output: html_document
---

```{r}
library(vcfR)
library(treemap)
library(readr)
require(devtools)
library(readr)
library(ggplot2)
library(reshape2)
library(pegas)
library(ape)  
library(ade4)
library(adegenet)
library(diveRsity)
library(genetics)
library(hierfstat)
library(iterators)
library(poppr)
library(readxl)
library(readr)
library("magrittr")
library(devtools)
library("mmod")
library(readr)
library(pegas)
library(ape)  
library(ade4)
library(genetics)
library(hierfstat)
library(iterators)
library(poppr)
library(readxl)
library(adegenet)
library(pegas)
library(StAMPP)
library(lattice)
library(ggplot2)
library(ape)
library(ggmap)
library(dplyr)
library(maps)
library("devtools")
library(StAMPP)
library(lattice)
library(ape)
library(ggmap)
```



```{r }
AU1 <- read.vcfR("~/Desktop/SNP.AUdp5p05FHWE2A.recode.vcf")
```

## Including Plots

You can also embed plots, for example:

```{r }
genindB <- vcfR2genind(AU1)
strata<- read.table("~/Desktop/popmap", header=TRUE)
strata_df <- data.frame(strata)
strata(genindB) <- strata_df
setPop(genindB) <- ~Population
```

```{r}
#Test the possibility of loci being in HWE 
cats.hwt <- hw.test(genindB, B=1000)

```

```{r}
## testing population structure
fstat(genindB)
wc(genindB)
```

```{r}
matFst <- pairwise.fst(genindB)
matFst
is.euclid(matFst)
```

```{r}
## Principal component analysis
X <- tab(genindB, freq = TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:20], main = "PCA eigenvalues", col = heat.colors(50))
col <- funky(9) 
s.class(pca1$li, pop(genindB),xax=1,yax=2, col=col, axesell=FALSE, cstar=0, cpoint=2, grid=FALSE)


```

```{r}
## DAPC plot
set.seed(20160308) # Setting a seed for a consistent result
grp <- find.clusters(genindB, max.n.clust = 10, n.pca = 25, choose.n.clust = FALSE) 
names(grp)
grp$grp
dapc1 <- dapc(genindB, grp$grp, n.pca = 25, n.da = 2) 
scatter(dapc1) # plot of the group
```

```{r}
### structure-like plot indicating membership probability of each sample
dapcB.results <- as.data.frame(dapc1$posterior)
dapcB.results$pop <- pop(genindB)
dapcB.results$indNames <- rownames(dapcB.results)
library(reshape2)
dapcB.results <- melt(dapcB.results)
colnames(dapcB.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
p <- ggplot(dapcB.results, aes(x=Sample,y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
p

# calculate mean per population
grouped_data <- group_by(dapcB.results, Original_Pop, Assigned_Pop)
data_means <- summarise(grouped_data, mean=mean(Posterior_membership_probability))
# plot means for each original pop colored by assigned pop
pieB <- ggplot(data_means, aes(x=Original_Pop,y=mean, fill=Assigned_Pop))
pieB <- pieB + geom_bar(stat='identity') + coord_polar("y",start=0)
pieB <- pieB + facet_grid(~Original_Pop, scales = "free") + theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid = element_blank(), strip.background = element_blank())
pieB
```

```{r}
## Plot phylogenetic Tree
titi <- genind2genpop(genindB)
Dgen <- dist.genpop(titi, method=2)
tree <- nj(Dgen)
plot.phylo(tree)
```

