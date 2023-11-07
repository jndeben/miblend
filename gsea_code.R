
# Gene Set Enrichment Analysis

```{r, GSEA}
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(tidyverse)
#BiocManager::install("ggridges")
library(ggridges)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("ReactomePA")
library(ReactomePA)
library(fgsea)
library(stats)
```

## Run GSEA for each comparison

```{r, GSEA code 1 , message=FALSE}
# set parameters
compare1 <- D0vsD1
pvalue <- .05

# create list of t values and gene list for the comparison
t_res1 <- na.omit(compare1[,c(1,6)])
t_res1 <- unlist(t_res1$t)
names(t_res1) <- make.names(na.omit(compare1[,1]), unique = T) %>% 
  gsub('^X', '', .)
t_res1 <- t_res1[order(t_res1, decreasing = T)]

Reactome1 <- ReactomePA::gsePathway(t_res1, 
                                    organism = "human",
                                    pvalueCutoff = pvalue,
                                    pAdjustMethod = "BH", 
                                    verbose = TRUE,
                                    eps = 1e-50,
                                    seed = T)
# view summary
datatable(Reactome@result)
```  

```{r, GSEA code 2}
# set parameters
compare2 <- D0vsD2
pvalue <- .05

# create list of t values and gene list for the comparison
t_res2 <- na.omit(compare2[,c(1,6)])
t_res2 <- unlist(t_res2$t)
names(t_res2) <- make.names(na.omit(compare2[,1]), unique = T) %>% 
  gsub('^X', '', .)
t_res2 <- t_res2[order(t_res2, decreasing = T)]

Reactome2 <- ReactomePA::gsePathway(t_res2, 
                                    organism = "human",
                                    pvalueCutoff = pvalue,
                                    pAdjustMethod = "BH", 
                                    verbose = TRUE,
                                    eps = 1e-50,
                                    seed = T)
# view summary
datatable(Reactome@result)
```


```{r, GSEA code 3}
# set parameters
compare3 <- D0vsD3
pvalue <- .05

# create list of t values and gene list for the comparison
t_res3 <- na.omit(compare3[,c(1,6)])
t_res3 <- unlist(t_res3$t)
names(t_res3) <- make.names(na.omit(compare3[,1]), unique = T) %>% 
  gsub('^X', '', .)
t_res3 <- t_res3[order(t_res3, decreasing = T)]

Reactome3 <- ReactomePA::gsePathway(t_res3, 
                                    organism = "human",
                                    pvalueCutoff = pvalue,
                                    pAdjustMethod = "BH", 
                                    verbose = TRUE,
                                    eps = 1e-50,
                                    seed = T)
# view summary
datatable(Reactome@result)
```

```{r, GSEA code 4}
# set parameters
compare4 <- D0vsD4
pvalue <- .05

# create list of t values and gene list for the comparison
t_res4 <- na.omit(compare4[,c(1,6)])
t_res4 <- unlist(t_res4$t)
names(t_res4) <- make.names(na.omit(compare4[,1]), unique = T) %>% 
  gsub('^X', '', .)
t_res4 <- t_res4[order(t_res4, decreasing = T)]

Reactome4 <- ReactomePA::gsePathway(t_res4, 
                                    organism = "human",
                                    pvalueCutoff = pvalue,
                                    pAdjustMethod = "BH", 
                                    verbose = TRUE,
                                    eps = 1e-50,
                                    seed = T)
# view summary
datatable(Reactome@result)
```

```{r, GSEA code 5}
# set parameters
compare5 <- D0vsD5
pvalue <- .05

# create list of t values and gene list for the comparison
t_res5 <- na.omit(compare5[,c(1,6)])
t_res5 <- unlist(t_res5$t)
names(t_res5) <- make.names(na.omit(compare5[,1]), unique = T) %>% 
  gsub('^X', '', .)
t_res5 <- t_res5[order(t_res5, decreasing = T)]

Reactome5 <- ReactomePA::gsePathway(t_res5, 
                                    organism = "human",
                                    pvalueCutoff = pvalue,
                                    pAdjustMethod = "BH", 
                                    verbose = TRUE,
                                    eps = 1e-50,
                                    seed = T)
# view summary
datatable(Reactome@result)
```

```{r, GSEA code 6}
# set parameters
compare6 <- D0vsD6
pvalue <- .05

# create list of t values and gene list for the comparison
t_res6 <- na.omit(compare6[,c(1,6)])
t_res6 <- unlist(t_res6$t)
names(t_res6) <- make.names(na.omit(compare6[,1]), unique = T) %>% 
  gsub('^X', '', .)
t_res6 <- t_res6[order(t_res6, decreasing = T)]

Reactome <- ReactomePA::gsePathway(t_res6, 
                                   organism = "human",
                                   pvalueCutoff = pvalue,
                                   pAdjustMethod = "BH", 
                                   verbose = TRUE,
                                   eps = 1e-50,
                                   seed = T)
# view summary
datatable(Reactome@result)
```


```{r, GSEA code 7}
# set parameters
compare7 <- D0vsD7
pvalue <- .05

# create list of t values and gene list for the comparison
t_res7 <- na.omit(compare7[,c(1,6)])
t_res7 <- unlist(t_res7$t)
names(t_res7) <- make.names(na.omit(compare7[,1]), unique = T) %>% 
  gsub('^X', '', .)
t_res7 <- t_res7[order(t_res7, decreasing = T)]

Reactome <- ReactomePA::gsePathway(t_res7, 
                                   organism = "human",
                                   pvalueCutoff = pvalue,
                                   pAdjustMethod = "BH", 
                                   verbose = TRUE,
                                   eps = 1e-50,
                                   seed = T)
# view summary
datatable(Reactome@result)
```
