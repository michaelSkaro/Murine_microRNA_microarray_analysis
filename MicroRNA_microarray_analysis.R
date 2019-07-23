# clear environment
rm(list=ls())

# get the libraries
library(Biobase)
library(oligoClasses)
library(knitr)
library(BiocStyle)
library(oligo)
library(geneplotter)
library(ggplot2)
library(dplyr)
library(LSD)
library(gplots)
library(RColorBrewer)
library(ArrayExpress)
library(arrayQualityMetrics)
library(stringr)
library(matrixStats)
library(topGO)
library(genefilter)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
library(pheatmap)
library(mvtnorm)
library(DAAG)
library(multcomp)
library(limma)
library(ReactomePA)
library(clusterProfiler)
library(openxlsx)
library(devtools)
library(biomaRt)
library(ArrayExpress)
#library(EnrichmentBrowser)
set.seed(777)
raw_data_dir <- file.path(getwd(), "rawDataMAWorkflow")


# get my GSEs
library(GEOquery)
my.gse <- "GSE116836"
##get published data and metadata
##this step is slow because of download and reformatting
##create directory for files downloaded from GEO
if(!file.exists("geo_downloads")) dir.create("geo_downloads")
if(!file.exists("results"))  dir.create("results", recursive=TRUE)
##get data from GEO
my.geo.gse <- getGEO(GEO=my.gse, filename=NULL, destdir="./geo_downloads", GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)


# check the type of bject tht was downloaded 

class(my.geo.gse)

# check the platform you want to analyze
names(my.geo.gse)

# get the one i want. In this case there is only one option however in larger studies, data might be submitted all at once and this is the 
# step to decide which platform we need. 
my.geo.gse <- my.geo.gse[[1]]


# the data is now in the form of an expression set. notice the class change in the object.

class(my.geo.gse)


# lets make sure the colnames match the documentation

colnames(pData(my.geo.gse))

# the expression data is in a S4 object. Therefore we have to call the correct facet of the data
# this will be accomplished by using the Expr() command

View(head(exprs(my.geo.gse)))


summary(exprs(my.geo.gse))
head(exprs(my.geo.gse))

pData(my.geo.gse)$data_processing[1]

# [1] Normalization of raw Ct values from qRT-PCR was performed using Exiqon's recommended software, GenEx, according to instructions provided: http://www.biomcc.com/genex-download.html
# Levels: Normalization of raw Ct values from qRT-PCR was performed using Exiqon's recommended software, GenEx, according to instructions provided: http://www.biomcc.com/genex-download.html

# We don't have the license to re-normalize or redo this analysis. We will have to take the result above to be true for now


# Begin the DE analysis.


dat <- exprs(my.geo.gse)
dim(dat)

# We will now group the data into the three appropriate groups. These groups are the mammary tissues that developed 
# tumors, other tissue locations that formed tumors, mice that were at risk and exposed to the ATXN and mice that 
# were not exposed and were true normals

design.mat <- data.table::fread("Design_matrix_microRNA.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(design.mat) <- c("geo_accession_id","sample_UUID","tissue_type","tumor_type","tumor_status")

# make the distringuishing column the factor and make a reference


design.mat$tumor_type <- as.character(design.mat$tumor_type)
design.mat$tumor_type[is.na(design.mat$tumor_type)] <- "healthy_serum"

design.mat$sample.levels <- as.factor(design.mat$tumor_type)
design.mat$sample.levels <- relevel(design.mat$sample.levels, ref="healthy_serum")


table(design.mat$sample.levels)

# Visualize the density plots of the matrix. Show data has been normalized
level.pal <- brewer.pal(4, "Dark2")
level.cols <- level.pal[unname(design.mat$sample.levels)]

plotDensities(dat, legend=F, col=level.cols , main="Arrays Normalized")
legend("topright", legend=levels(design.mat$sample.levels), fill=level.pal)


# boxplot

boxplot(dat, las=2, names=design.mat$sample.levels, outline=F, col=level.cols, main="Arrays Normalized")


# with the data showing the variability it is and the description of the normalization techniques,
# I am either going to assume the data is not normalized, the data is not a microarray or was improperly 
# normalized. I will finish the analysis and conduct the work as if it is correct, as shown in line 84 and 85
# then I will go back and download the raw files. I will write my own software to normalize the data as a qRT-PCR experiment
# or normalize it as a microarray experiment and compare results. While these will not be the reccomended settings set but 
# Exiqon, I will get some solid results pertaining to the identity of the data and if any DGE be relied on. 



# conducting DGE analysis, clustering, pathway enrichment and visualization


# PCA plot

plotMDS(dat, labels=design.mat$sample.levels, gene.selection="common", col = level.cols ,main="MDS Plot to Compare Replicates")

# this result clearly shows the differences in this data exists within the disease state and not the tissue type


# clustering samples on euclidian distance. This will require some dimensionality reduction. As we want to avoid
# highly expressed genes dominating the clustering we will use the normalized data and conduct z-score values on
# the number of SD each value in the sample is away from the sample mean. samples with similar vairances and covariances
# will cluster highly together and those that differ because of disease state or tissue type. 


# hierarchical clustering

cluster.dat <- dat
gene.mean <- apply(cluster.dat, 1, mean)
gene.sd <- apply(cluster.dat, 1, sd)
cluster.dat <- sweep(cluster.dat, 1, gene.mean, "-")
cluster.dat <- sweep(cluster.dat, 1, gene.sd, "/")

my.dist <- dist(t(cluster.dat), method="euclidean")
my.hclust <- hclust(my.dist, method="average")
my.hclust$labels <- design.mat$sample.levels
plot(my.hclust, cex=0.75, main="Comparison of Biological Replicates", xlab="Euclidean Distance")


# clustering based on pearson correlation

my.cor <- cor(dat)
my.dist <- as.dist(1 - my.cor)
my.hclust <- hclust(my.dist, method="average")
my.hclust$labels <- design.mat$sample.levels
plot(my.hclust, cex=0.75, main="Comparison of Biological Replicates")


# heatmaps

library(ComplexHeatmap)

ComplexHeatmap::Heatmap(my.cor, name = "Correlation_of_expresssion", split =2, show_column_dend = FALSE)

ComplexHeatmap::Heatmap(dat, name = "normalized expression", show_row_names = FALSE)



# linear modeling

# new designmatrix with binary operators

design.mat$sample.levels <- str_replace_all(design.mat$sample.levels,"[[:punct:]]","_")
design.mat$sample.levels <- str_replace_all(design.mat$sample.levels," ","_")
design.mat$sample.levels <- str_replace_all(design.mat$sample.levels,"__","_")


my.design <- model.matrix(~0 + sample.levels, design.mat)

colnames(my.design) <- c("at_risk", "healthy_serum", "Liver_cancer", "mammary_tumor")
colnames(my.design)
rownames(my.design) <- design.mat$sample.levels


# linear model fitting

my.fit <- lmFit(dat, my.design)
write.table(my.fit$coefficients, file=paste0("results/",my.gse,"_Limma_Coeff.txt"), sep="\t", quote=F)



# saved here,everything up until this line is okay**

# make contrast groups
my.contrasts <- makeContrasts(HSvLC = healthy_serum - Liver_cancer, HSvMT= healthy_serum - mammary_tumor, 
                              HSvAR = healthy_serum - at_risk, levels = my.design)

contrast.fits <- sapply(colnames(my.contrasts), function(x)(contrasts.fit(my.fit, contrasts=my.contrasts[, x])))
length(contrast.fits)

colnames(contrast.fits)

contrast.fits


# Bayes work: in this sction we will test the true differences in expression by contructing a bayesian network. These tests are
# used in favor of multiple ttests to correct for multiple testing. We will use an adjusted proportion of 5% of the genes expected
# expected to be DE. For larger studies I will reduce this to a proportion of 1% to reduce the FDR.



mr.bayes <- eBayes(my.fit, proportion=0.01, stdev.coef.lim=c(0.1,4), trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))
contrast.tts <- lapply(mr.bayes, function(x)(topTable(x, adjust="BH", number=length(x$coefficients), sort.by="none")))







# Bayesian modeling

contrast.ebs <- lapply(contrast.fits, function(x)(eBayes(x, proportion=0.1, trend=FALSE, robust=FALSE)))
contrast.tts <- lapply(contrast.ebs, function(x)(topTable(x, adjust="BH", number=length(x$coefficients), sort.by="none")))
contrast.tests <- lapply(contrast.ebs, function(x)(decideTests(x, method="separate", adjust.method="BH", p.value=0.05, lfc=0)))







# Visualization 

# generate plots(heatmap,PCA, Volcano) for all of the sequencings results for each cancer type.
library(EnhancedVolcano)
library(pheatmap)
library(EnhancedVolcano)


res <- data.table::fread("DE_Liver_vs_healthy.txt", header = TRUE)

  # volcano
  
EnhancedVolcano(res,
                lab = res$ID,
                x = 'logFC',
                y = 'P.Value',
                xlim = c(-8, 8),
                title = 'Liver Cancer vs. Healthy Controls',
                pCutoff = 0.05,
                FCcutoff = .5,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0)
  
res <- data.table::fread("DE_mammary_vs_healthy.txt", header = TRUE)

# volcano

EnhancedVolcano(res,
                lab = res$ID,
                x = 'logFC',
                y = 'P.Value',
                xlim = c(-8, 8),
                title = 'Mammary Cancer vs. Healthy Controls',
                pCutoff = 0.05,
                FCcutoff = .5,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0)

res <- data.table::fread("DE_at_risk_vs_healthy.txt", header = TRUE)

# volcano

EnhancedVolcano(res,
                lab = res$ID,
                x = 'logFC',
                y = 'P.Value',
                xlim = c(-8, 8),
                title = 'At Risk vs. Healthy Controls',
                pCutoff = 0.05,
                FCcutoff = .5,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0)


res <- data.table::fread("DE_Liver_vs_at_risk.txt", header = TRUE)

# volcano

EnhancedVolcano(res,
                lab = res$ID,
                x = 'logFC',
                y = 'P.Value',
                xlim = c(-8, 8),
                title = 'Liver cancer vs. At risk',
                pCutoff = 0.05,
                FCcutoff = .5,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0)

res <- data.table::fread("DE_Mammary_vs_at_risk.txt", header = TRUE)

# volcano

EnhancedVolcano(res,
                lab = res$ID,
                x = 'logFC',
                y = 'P.Value',
                xlim = c(-8, 8),
                title = 'Mammary cancer vs. At risk',
                pCutoff = 0.05,
                FCcutoff = .5,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0)

res <- data.table::fread("DE_Liver_vs_Mammary.txt", header = TRUE)

# volcano

EnhancedVolcano(res,
                lab = res$ID,
                x = 'logFC',
                y = 'P.Value',
                xlim = c(-8, 8),
                title = 'Liver cancer vs. Mammary cancer',
                pCutoff = 0.05,
                FCcutoff = .5,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0)

# Gene associations drug interactions and disease drivers
# in this section we will use the multiMIR package to investigate the 
# gene interactions of DE microRNAs. Further we will use this data as 
# input data into enriched pathway investigation

#BiocManager::install("multiMiR")
library(multiMiR)
library(edgeR)

# # No significantly DE in this section
res <- data.table::fread("DE_Liver_vs_healthy.txt", header = TRUE)
res <- res[res$adj.P.Val < 0.05,]
res <- res$miRNA_ID

multimir_results <- get_multimir(org     = 'mmu',
                                 mirna   = res,
                                 table   = 'validated',
                                 summary = TRUE)

# no join in this section. # No significantly DE in this section

# #################
# # done
res <- data.table::fread("DE_mammary_vs_healthy.txt", header = TRUE)
res <- res[res$adj.P.Val < 0.05,]
#res <- res$miRNA_ID

multimir_results <- get_multimir(org     = 'mmu',
                                 mirna   = res,
                                 table   = 'validated',
                                 summary = TRUE)
Mammary_vs_Healthy_DE_microRNAs<- multimir_results@data

colnames(Mammary_vs_Healthy_DE_microRNAs)

Mammary_vs_Healthy_DE_microRNAs <- left_join(
  x = Mammary_vs_Healthy_DE_microRNAs,
  y = res,
  by = c("mature_mirna_id" = "miRNA_ID")
)

write.csv(Mammary_vs_Healthy_DE_microRNAs, file = "Mammary_vs_Healthy_DE_microRNAs.csv")

###################
# # no join in this section. # No significantly DE in this section
res <- data.table::fread("DE_Mammary_vs_at_risk.txt", header = TRUE)
res <- res[res$adj.P.Val < 0.05,]
res <- res$miRNA_ID

multimir_results <- get_multimir(org     = 'mmu',
                                 mirna   = res,
                                 table   = 'validated',
                                 summary = TRUE)
dat<- multimir_results@data

dat <- left_join(
  x = dat,
  y = res,
  by = c("mature_mirna_id" = "miRNA_ID")
)



write.csv(dat, file = "Mammary_vs_at_risk_DE_microRNAs.csv")

###################
# no sigDE
# no join in this section. # No significantly DE in this section
# res <- data.table::fread("DE_Liver_vs_Mammary.txt", header = TRUE)
# res <- res[res$adj.P.Val < 0.05,]
# res <- res$miRNA_ID
#
# multimir_results <- get_multimir(org     = 'mmu',
#                                  mirna   = res,
#                                  table   = 'validated',
#                                  summary = TRUE)
# dat<- multimir_results@data
#
# write.csv(dat, file = "Mammary_vs_Liver_DE_microRNAs.csv")

###################
res <- data.table::fread("DE_Liver_vs_at_risk.txt", header = TRUE)
res <- res[res$adj.P.Val < 0.05,]
res <- res$miRNA_ID

multimir_results <- get_multimir(org     = 'mmu',
                                 mirna   = res,
                                 table   = 'validated',
                                 summary = TRUE)
dat<- multimir_results@data

dat <- left_join(
  x = dat, 
  y = res,
  by = c("mature_mirna_id" = "miRNA_ID")
)


write.csv(dat, file = "Liver_vs_at_risk_DE_microRNAs.csv")

###################
res <- data.table::fread("DE_at_risk_vs_healthy.txt", header = TRUE)
res <- res[res$adj.P.Val < 0.05,]
res <- res$miRNA_ID

multimir_results <- get_multimir(org     = 'mmu',
                                 mirna   = res,
                                 table   = 'validated',
                                 summary = TRUE)
dat<- multimir_results@data

dat <- left_join(
  x = dat, 
  y = res,
  by = c("mature_mirna_id" = "miRNA_ID")
)





write.csv(dat, file = "at_risk_vs_healthy_DE_microRNAs.csv")


# Ranked data by the pval and fold change, left join by the enriched 
# microRNAs in the multmir data s3 object. 

# use the candidate microRNAs to construct GSEA relying on pval and predicted target gene sets 
# we will use cluster profiler here


library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(DOSE)
# in this section we will use the the clusterprofiler package to provide enriched pathway information
# in an effort to reduce our candidates

dat <- data.table::fread("Mammary_vs_Healthy_DE_microRNAs.csv", header = TRUE, stringsAsFactors = FALSE)
dat$V1 <- NULL


egoMF<- enrichGO(gene = dat$target_symbol, OrgDb = org.Mm.eg.db , keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH")
egoBP<- enrichGO(gene = dat$target_symbol, OrgDb = org.Mm.eg.db , keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH")
egoCC<- enrichGO(gene = dat$target_symbol, OrgDb = org.Mm.eg.db , keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH")

egoMF <- simplify(egoMF)
#egoBP <- simplify(egoBP)
egoCC <- simplify(egoCC)

#save res
write.csv(egoMF@result, "Dysregulated_pathways_molecular_functions_Mammary_vs_healthy.csv")
write.csv(egoBP@result, "Dysregulated_pathways_biological_functions_Mammary_vs_healthy.csv")
write.csv(egoCC@result, "Dysregulated_pathways_cellular_compartments_Mammary_vs_healthy.csv")


p1 <- enrichplot::dotplot(egoMF, x = "Count", showCategory =40, title = "Mammary tumors vs. Healthy controls")
p <- clusterProfiler::dotplot(egoBP, showCategory =20)
p <- clusterProfiler::dotplot(egoCC, showCategory =20)

heatplot(egoMF, foldChange=2^(dat$logFC))


# library(DOSE)
# data(geneList)
# dat.foo <- names(geneList)[abs(geneList) > 2]
# dat.foo.MF <- enrichGO(dat.foo, org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH")
# clusterProfiler::cnetplot(dat.foo.MF, showCategory = 20)
# # as the CNE plot is simply another version of the upsetplot according to the document vignette,
# # this object should be in the proper input form for the upsetplot, correct?
# upsetplot(dat.foo.MF)

# compare to TCGA data to find similar expression signatures in our mammary tumors 

# subset BRCA patients to ones over expressing ENPP2, then subset their expression to only 
# the microRNAs in the human microNOME. use differenital expression analysis to investigte 
# up and downregulated pathways in these patients. finally compare them to the autotaxin mice.

# :)

library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)




projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

#annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]
#tumor.samples <- tumor.samples[tumor.samples$caseID %in% normal.samples$caseID,]


proj <- projects[2]

df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv")) %>%
  as_tibble() %>%
  column_to_rownames(var = "Ensembl")

coldata.t <- tumor.samples[tumor.samples$project == proj,]
coldata.n <- normal.samples[normal.samples$project == proj,]


# subset the samples to breast cancer tumor and normals

df.exp.t <- df.exp %>%
  dplyr::select(coldata.t$barcode)

df.exp.n <- df.exp %>%
  dplyr::select(coldata.n$barcode)

# subset to ones that have significantly higher expression
# of the ENPP2 compared to normals, note:: should be about 25%
# of the cohort

# maybe not?

rownames(df.exp.t) <- substr(rownames(df.exp.t), 1,15)
rownames(df.exp.n) <- substr(rownames(df.exp.n), 1,15)

ENPP2 <- "ENSG00000136960"

df.exp.t.ENPP2 <- df.exp.t[rownames(df.exp.t) == ENPP2,]
df.exp.n.ENPP2 <- df.exp.n[rownames(df.exp.n) == ENPP2,]

df.exp.t.ENPP2 <- as.data.frame(t(df.exp.t.ENPP2))
df.exp.n.ENPP2 <- as.data.frame(t(df.exp.n.ENPP2))

df.exp.t.ENPP2 <- df.exp.t.ENPP2 %>%
  rownames_to_column("barcode")

df.exp.n.ENPP2 <- df.exp.n.ENPP2 %>%
  rownames_to_column("barcode")

normcutoff<- summary(df.exp.n.ENPP2$ENSG00000136960)

# subset the data for only those that are significantly greater than the mean

df.exp.t.ENPP2 <- df.exp.t.ENPP2 %>%
  filter(ENSG00000136960 > normcutoff[2] )


# convert the microRNAs to the ENSEMBL gene IDs (2655 uniq) using bioMaRt

library(biomaRt)
HSA_microRNAs <- data.table::fread("HSA_microRNAs.txt", header = FALSE)
colnames(HSA_microRNAs) <- "microRNAs"
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
converts<- getBM(attributes = c("mirbase_id","ensembl_gene_id"), filters = "mirbase_id", values = HSA_microRNAs$microRNAs ,mart = ensembl)


# subset the expression of the patients to the microRNAs

df.exp.t <- df.exp.t %>%
  dplyr::select(df.exp.t.ENPP2$barcode)

df.exp.t <- df.exp.t[rownames(df.exp.t) %in% converts$ensembl_gene_id,]


df.exp.n <- df.exp.n %>%
  dplyr::select(df.exp.n.ENPP2$barcode)

df.exp.n <- df.exp.n[rownames(df.exp.n) %in% converts$ensembl_gene_id,]


# conduct Differential expression of the microRNAs in the autotaxin tumors vs the healthy patients

df.exp.f <- cbind(df.exp.n, df.exp.t)

coldata.t <- coldata.t[coldata.t$barcode %in% colnames(df.exp.t),]
coldata.n <- coldata.n[coldata.n$barcode %in% colnames(df.exp.n),]


coldata <- rbind(coldata.n, coldata.t)
rownames(coldata) <- coldata$barcode

# df.exp <- df.exp[ ,colnames(df.exp) %in% coldata$barcode]

coldata$sample_type <- gsub(" ", "_", x = coldata$sample_type)


dds <- DESeqDataSetFromMatrix(countData = df.exp.f, colData = coldata, design = ~ sample_type)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


dds$sample_type <- relevel(dds$sample_type, ref = "Solid_Tissue_Normal")


dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
resOrdered <- as.data.frame(resOrdered)
res <- as.data.frame(res)

resOrdered <- resOrdered[complete.cases(resOrdered), ]

resOrdered <- as.data.frame(resOrdered) %>%
  rownames_to_column("Ensembl")

write.csv(resOrdered, file = str_glue("~/storage/Metastatic_Organo_Tropism/microRNA_characterization/Differential_expression_RT-PCR_exiqon_micrRNA_microarray/{proj}_microRNA_DE.csv"))


# find the enriched gene regulations in the population

View(resOrdered)

resOrdered <- resOrdered[resOrdered$padj <= 1.0e-4,]

resOrdered.miRNA <- left_join(
  x = resOrdered, 
  y = converts,
  by = c("Ensembl" = "ensembl_gene_id")
)

multimir_results <- get_multimir(org     = 'hsa',
                                 mirna   = resOrdered.miRNA$mirbase_id,
                                 table   = 'validated',
                                 summary = TRUE)
miRNA.dat<- multimir_results@data

# miRNA.dat.final <- left_join(
#   x = miRNA.dat, 
#   y = resOrdered.miRNA,
#   by = c("mature_mirna_id" = "mirbase_id")
#   
#   
# )

# use cluster profiler to show the enriched pathways in the autaxin enriched patients

library(org.Hs.eg.db)
resOrdered <- resOrdered[resOrdered$padj <= 1.0e-4,]  

egoMF<- enrichGO(gene = miRNA.dat$target_ensembl, OrgDb = org.Hs.eg.db , keyType = "ENSEMBL", ont = "MF", pAdjustMethod = "BH")
egoBP<- enrichGO(gene = miRNA.dat$target_ensembl, OrgDb = org.Hs.eg.db , keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH")
egoCC<- enrichGO(gene = miRNA.dat$target_ensembl, OrgDb = org.Hs.eg.db , keyType = "ENSEMBL", ont = "CC", pAdjustMethod = "BH")

egoMF <- simplify(egoMF)
egoBP <- simplify(egoBP)
egoCC <- simplify(egoCC)

#save res
write.csv(egoMF@result, "Dysregulated_pathways_molecular_functions_BRCA_vs_Normal.csv")
write.csv(egoBP@result, "Dysregulated_pathways_biological_functions_BRCA_vs_Normal.csv")
write.csv(egoCC@result, "Dysregulated_pathways_cellular_compartments_BRCA_vs_Normal.csv")




p2 <- enrichplot::dotplot(egoMF, x = "Count", showCategory =40, title = "TCGA-BRCA_vs_Normal")
p <- clusterProfiler::dotplot(egoBP, showCategory =20)
p <- clusterProfiler::dotplot(egoCC, showCategory =20)




# venn diagram describing enriched pathways

# library(VennDiagram)
# 
# TCGA_BRCA <- data.table::fread("Dysregulated_pathways_molecular_functions_BRCA_vs_Normal.csv", header = TRUE)
# ATXN_Mammary_tumors <- data.table::fread("Dysregulated_pathways_molecular_functions_Mammary_vs_healthy.csv", header = TRUE)
# 
# TCGA_BRCA <- TCGA_BRCA[TCGA_BRCA$pvalue <= 0.05,]
# ATXN_Mammary_tumors <- ATXN_Mammary_tumors[ATXN_Mammary_tumors$pvalue <= 0.05,]
# 
# TCGA_BRCA <- TCGA_BRCA$Description
# ATXN_Mammary_tumors <- ATXN_Mammary_tumors$Description
# 
# x <- list(TCGA_BRCA = TCGA_BRCA, ATXN_Mammary_tumors = ATXN_Mammary_tumors)
# 
# 
# v0 <- venn.diagram( x, filename=NULL,
#                     fill = c("red", "blue"), 
#                     alpha = 0.50,col = "transparent")
# 
# grid.draw(v0)




session_info()


# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                       
# version  R version 3.5.1 (2018-07-02)
# os       Debian GNU/Linux 9 (stretch)
# system   x86_64, linux-gnu           
# ui       RStudio                     
# language (EN)                        
# collate  en_US.UTF-8                 
# ctype    en_US.UTF-8                 
# tz       Etc/UTC                     
# date     2019-06-28                  
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                        * version   date       lib source        
# affxparser                       1.54.0    2018-10-30 [1] Bioconductor  
# affy                           * 1.60.0    2018-10-30 [1] Bioconductor  
# affyio                           1.52.0    2018-10-30 [1] Bioconductor  
# annotate                       * 1.60.1    2019-03-07 [1] Bioconductor  
# AnnotationDbi                  * 1.44.0    2018-10-30 [1] Bioconductor  
# AnnotationHub                    2.14.5    2019-03-14 [1] Bioconductor  
# aroma.light                      3.12.0    2018-10-30 [1] Bioconductor  
# ArrayExpress                   * 1.42.0    2018-10-30 [1] Bioconductor  
# assertthat                       0.2.1     2019-03-21 [1] CRAN (R 3.5.1)
# backports                        1.1.4     2019-04-10 [1] CRAN (R 3.5.1)
# Biobase                        * 2.42.0    2018-10-30 [1] Bioconductor  
# BiocGenerics                   * 0.28.0    2018-10-30 [1] Bioconductor  
# BiocManager                      1.30.4    2018-11-13 [1] CRAN (R 3.5.1)
# BiocParallel                     1.16.6    2019-02-10 [1] Bioconductor  
# BiocStyle                      * 2.10.0    2018-10-30 [1] Bioconductor  
# biomaRt                        * 2.38.0    2018-10-30 [1] Bioconductor  
# Biostrings                     * 2.50.2    2019-01-03 [1] Bioconductor  
# bit                              1.1-14    2018-05-29 [1] CRAN (R 3.5.1)
# bit64                            0.9-7     2017-05-08 [1] CRAN (R 3.5.1)
# bitops                           1.0-6     2013-08-17 [1] CRAN (R 3.5.1)
# blob                             1.1.1     2018-03-25 [1] CRAN (R 3.5.1)
# broom                            0.5.2     2019-04-07 [1] CRAN (R 3.5.1)
# callr                            3.2.0     2019-03-15 [1] CRAN (R 3.5.1)
# caTools                          1.17.1.2  2019-03-06 [1] CRAN (R 3.5.1)
# checkmate                        1.9.3     2019-05-03 [1] CRAN (R 3.5.1)
# circlize                         0.4.6     2019-04-03 [1] CRAN (R 3.5.1)
# cli                              1.1.0     2019-03-19 [1] CRAN (R 3.5.1)
# cluster                          2.0.7-1   2018-04-13 [2] CRAN (R 3.5.1)
# clusterProfiler                * 3.10.1    2018-12-20 [1] Bioconductor  
# cmprsk                           2.2-7     2014-06-17 [1] CRAN (R 3.5.1)
# codetools                        0.2-15    2016-10-05 [2] CRAN (R 3.5.1)
# colorspace                       1.4-1     2019-03-18 [1] CRAN (R 3.5.1)
# ComplexHeatmap                 * 1.20.0    2018-10-30 [1] Bioconductor  
# ConsensusClusterPlus             1.46.0    2018-10-30 [1] Bioconductor  
# cowplot                          0.9.4     2019-01-08 [1] CRAN (R 3.5.1)
# crayon                           1.3.4     2017-09-16 [1] CRAN (R 3.5.1)
# curl                             3.3       2019-01-10 [1] CRAN (R 3.5.1)
# DAAG                           * 1.22.1    2019-03-02 [1] CRAN (R 3.5.1)
# data.table                       1.12.2    2019-04-07 [1] CRAN (R 3.5.1)
# DBI                            * 1.0.0     2018-05-02 [1] CRAN (R 3.5.1)
# DelayedArray                     0.8.0     2018-10-30 [1] Bioconductor  
# desc                             1.2.0     2018-05-01 [1] CRAN (R 3.5.1)
# DESeq                            1.34.1    2019-01-04 [1] Bioconductor  
# devtools                       * 2.0.2     2019-04-08 [1] CRAN (R 3.5.1)
# digest                           0.6.19    2019-05-20 [1] CRAN (R 3.5.1)
# DNAcopy                          1.56.0    2018-10-30 [1] Bioconductor  
# DO.db                            2.9       2018-11-20 [1] Bioconductor  
# doParallel                       1.0.14    2018-09-24 [1] CRAN (R 3.5.1)
# DOSE                             3.8.2     2019-01-14 [1] Bioconductor  
# downloader                       0.4       2015-07-09 [1] CRAN (R 3.5.1)
# dplyr                          * 0.8.1     2019-05-14 [1] CRAN (R 3.5.1)
# EDASeq                           2.16.3    2019-01-10 [1] Bioconductor  
# edgeR                            3.24.3    2019-01-02 [1] Bioconductor  
# enrichplot                       1.2.0     2018-10-30 [1] Bioconductor  
# europepmc                        0.3       2018-04-20 [1] CRAN (R 3.5.1)
# evaluate                         0.14      2019-05-28 [1] CRAN (R 3.5.1)
# ExperimentHub                    1.8.0     2018-10-30 [1] Bioconductor  
# farver                           1.1.0     2018-11-20 [1] CRAN (R 3.5.1)
# fastmatch                        1.1-0     2017-01-28 [1] CRAN (R 3.5.1)
# ff                               2.2-14    2018-05-15 [1] CRAN (R 3.5.1)
# fgsea                            1.8.0     2018-10-30 [1] Bioconductor  
# foreach                          1.4.4     2017-12-12 [1] CRAN (R 3.5.1)
# fs                               1.3.1     2019-05-06 [1] CRAN (R 3.5.1)
# gdata                            2.18.0    2017-06-06 [1] CRAN (R 3.5.1)
# genefilter                     * 1.64.0    2018-10-30 [1] Bioconductor  
# geneplotter                    * 1.60.0    2018-10-30 [1] Bioconductor  
# generics                         0.0.2     2018-11-29 [1] CRAN (R 3.5.1)
# GenomeInfoDb                     1.18.2    2019-02-12 [1] Bioconductor  
# GenomeInfoDbData                 1.2.0     2018-11-20 [1] Bioconductor  
# GenomicAlignments                1.18.1    2019-01-04 [1] Bioconductor  
# GenomicFeatures                  1.34.8    2019-04-10 [1] Bioconductor  
# GenomicRanges                    1.34.0    2018-10-30 [1] Bioconductor  
# GEOquery                       * 2.50.5    2018-12-22 [1] Bioconductor  
# GetoptLong                       0.1.7     2018-06-10 [1] CRAN (R 3.5.1)
# ggforce                          0.2.2     2019-04-23 [1] CRAN (R 3.5.1)
# ggplot2                        * 3.1.1     2019-04-07 [1] CRAN (R 3.5.1)
# ggplotify                        0.0.3     2018-08-03 [1] CRAN (R 3.5.1)
# ggpubr                           0.2       2018-11-15 [1] CRAN (R 3.5.1)
# ggraph                           1.0.2     2018-07-07 [1] CRAN (R 3.5.1)
# ggrepel                          0.8.1     2019-05-07 [1] CRAN (R 3.5.1)
# ggridges                         0.5.1     2018-09-27 [1] CRAN (R 3.5.1)
# ggthemes                         4.2.0     2019-05-13 [1] CRAN (R 3.5.1)
# GlobalOptions                    0.1.0     2018-06-09 [1] CRAN (R 3.5.1)
# glue                             1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
# GO.db                          * 3.7.0     2019-03-11 [1] Bioconductor  
# GOSemSim                         2.8.0     2018-10-30 [1] Bioconductor  
# gplots                         * 3.0.1.1   2019-01-27 [1] CRAN (R 3.5.1)
# graph                          * 1.60.0    2018-10-30 [1] Bioconductor  
# graphite                         1.28.2    2019-01-18 [1] Bioconductor  
# gridExtra                        2.3       2017-09-09 [1] CRAN (R 3.5.1)
# gridGraphics                     0.4-1     2019-05-20 [1] CRAN (R 3.5.1)
# gtable                           0.3.0     2019-03-25 [1] CRAN (R 3.5.1)
# gtools                           3.8.1     2018-06-26 [1] CRAN (R 3.5.1)
# highr                            0.8       2019-03-20 [1] CRAN (R 3.5.1)
# hms                              0.4.2     2018-03-10 [1] CRAN (R 3.5.1)
# htmltools                        0.3.6     2017-04-28 [1] CRAN (R 3.5.1)
# httpuv                           1.5.1     2019-04-05 [1] CRAN (R 3.5.1)
# httr                             1.4.0     2018-12-11 [1] CRAN (R 3.5.1)
# hugene10sttranscriptcluster.db * 8.7.0     2019-06-13 [1] Bioconductor  
# hwriter                          1.3.2     2014-09-10 [1] CRAN (R 3.5.1)
# igraph                           1.2.4.1   2019-04-22 [1] CRAN (R 3.5.1)
# interactiveDisplayBase           1.20.0    2018-10-30 [1] Bioconductor  
# IRanges                        * 2.16.0    2018-10-30 [1] Bioconductor  
# iterators                        1.0.10    2018-07-13 [1] CRAN (R 3.5.1)
# jsonlite                         1.6       2018-12-07 [1] CRAN (R 3.5.1)
# KernSmooth                       2.23-15   2015-06-29 [2] CRAN (R 3.5.1)
# km.ci                            0.5-2     2009-08-30 [1] CRAN (R 3.5.1)
# KMsurv                           0.1-5     2012-12-03 [1] CRAN (R 3.5.1)
# knitr                          * 1.23      2019-05-18 [1] CRAN (R 3.5.1)
# later                            0.8.0     2019-02-11 [1] CRAN (R 3.5.1)
# lattice                        * 0.20-38   2018-11-04 [2] CRAN (R 3.5.1)
# latticeExtra                     0.6-28    2016-02-09 [1] CRAN (R 3.5.1)
# lazyeval                         0.2.2     2019-03-15 [1] CRAN (R 3.5.1)
# limma                          * 3.38.3    2018-12-02 [1] Bioconductor  
# locfit                           1.5-9.1   2013-04-20 [1] CRAN (R 3.5.1)
# LSD                            * 4.0-0     2018-01-26 [1] CRAN (R 3.5.1)
# magrittr                         1.5       2014-11-22 [1] CRAN (R 3.5.1)
# MASS                           * 7.3-51.1  2018-11-01 [2] CRAN (R 3.5.1)
# matlab                           1.0.2     2014-06-24 [1] CRAN (R 3.5.1)
# Matrix                           1.2-15    2018-11-01 [2] CRAN (R 3.5.1)
# matrixStats                    * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
# memoise                          1.1.0     2017-04-21 [1] CRAN (R 3.5.1)
# mgcv                             1.8-25    2018-10-26 [2] CRAN (R 3.5.1)
# mime                             0.7       2019-06-11 [1] CRAN (R 3.5.1)
# multcomp                       * 1.4-10    2019-03-05 [1] CRAN (R 3.5.1)
# munsell                          0.5.0     2018-06-12 [1] CRAN (R 3.5.1)
# mvtnorm                        * 1.0-10    2019-03-05 [1] CRAN (R 3.5.1)
# nlme                             3.1-137   2018-04-07 [2] CRAN (R 3.5.1)
# oligo                          * 1.46.0    2018-10-30 [1] Bioconductor  
# oligoClasses                   * 1.44.0    2018-10-30 [1] Bioconductor  
# openxlsx                       * 4.1.0.1   2019-05-28 [1] CRAN (R 3.5.1)
# org.Hs.eg.db                   * 3.7.0     2019-01-16 [1] Bioconductor  
# pd.hugene.1.0.st.v1            * 3.14.1    2019-06-13 [1] Bioconductor  
# pheatmap                       * 1.0.12    2019-01-04 [1] CRAN (R 3.5.1)
# pillar                           1.4.1     2019-05-28 [1] CRAN (R 3.5.1)
# pkgbuild                         1.0.3     2019-03-20 [1] CRAN (R 3.5.1)
# pkgconfig                        2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
# pkgload                          1.0.2     2018-10-29 [1] CRAN (R 3.5.1)
# plyr                             1.8.4     2016-06-08 [1] CRAN (R 3.5.1)
# polyclip                         1.10-0    2019-03-14 [1] CRAN (R 3.5.1)
# preprocessCore                   1.44.0    2018-10-30 [1] Bioconductor  
# prettyunits                      1.0.2     2015-07-13 [1] CRAN (R 3.5.1)
# processx                         3.3.1     2019-05-08 [1] CRAN (R 3.5.1)
# progress                         1.2.2     2019-05-16 [1] CRAN (R 3.5.1)
# promises                         1.0.1     2018-04-13 [1] CRAN (R 3.5.1)
# ps                               1.3.0     2018-12-21 [1] CRAN (R 3.5.1)
# purrr                            0.3.2     2019-03-15 [1] CRAN (R 3.5.1)
# qvalue                           2.14.1    2019-01-10 [1] Bioconductor  
# R.methodsS3                      1.7.1     2016-02-16 [1] CRAN (R 3.5.1)
# R.oo                             1.22.0    2018-04-22 [1] CRAN (R 3.5.1)
# R.utils                          2.9.0     2019-06-13 [1] CRAN (R 3.5.1)
# R6                               2.4.0     2019-02-14 [1] CRAN (R 3.5.1)
# randomForest                     4.6-14    2018-03-25 [1] CRAN (R 3.5.1)
# rappdirs                         0.3.1     2016-03-28 [1] CRAN (R 3.5.1)
# RColorBrewer                   * 1.1-2     2014-12-07 [1] CRAN (R 3.5.1)
# Rcpp                             1.0.1     2019-03-17 [1] CRAN (R 3.5.1)
# RCurl                            1.95-4.12 2019-03-04 [1] CRAN (R 3.5.1)
# reactome.db                      1.66.0    2019-01-23 [1] Bioconductor  
# ReactomePA                     * 1.26.0    2018-10-30 [1] Bioconductor  
# readr                            1.3.1     2018-12-21 [1] CRAN (R 3.5.1)
# remotes                          2.0.4     2019-04-10 [1] CRAN (R 3.5.1)
# reshape2                         1.4.3     2017-12-11 [1] CRAN (R 3.5.1)
# rjson                            0.2.20    2018-06-08 [1] CRAN (R 3.5.1)
# rlang                            0.3.4     2019-04-07 [1] CRAN (R 3.5.1)
# rmarkdown                        1.13      2019-05-22 [1] CRAN (R 3.5.1)
# rprojroot                        1.3-2     2018-01-03 [1] CRAN (R 3.5.1)
# Rsamtools                        1.34.1    2019-01-31 [1] Bioconductor  
# RSQLite                        * 2.1.1     2018-05-06 [1] CRAN (R 3.5.1)
# rstudioapi                       0.10      2019-03-19 [1] CRAN (R 3.5.1)
# rtracklayer                      1.42.2    2019-03-01 [1] Bioconductor  
# rvcheck                          0.1.3     2018-12-06 [1] CRAN (R 3.5.1)
# rvest                            0.3.4     2019-05-15 [1] CRAN (R 3.5.1)
# S4Vectors                      * 0.20.1    2018-11-09 [1] Bioconductor  
# sandwich                         2.5-1     2019-04-06 [1] CRAN (R 3.5.1)
# scales                           1.0.0     2018-08-09 [1] CRAN (R 3.5.1)
# selectr                          0.4-1     2018-04-06 [1] CRAN (R 3.5.1)
# sesame                           1.0.0     2018-10-30 [1] Bioconductor  
# sesameData                       1.0.0     2018-11-01 [1] Bioconductor  
# sessioninfo                      1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
# shape                            1.4.4     2018-02-07 [1] CRAN (R 3.5.1)
# shiny                            1.3.2     2019-04-22 [1] CRAN (R 3.5.1)
# ShortRead                        1.40.0    2018-10-30 [1] Bioconductor  
# SparseM                        * 1.77      2017-04-23 [1] CRAN (R 3.5.1)
# stringi                          1.4.3     2019-03-12 [1] CRAN (R 3.5.1)
# stringr                        * 1.4.0     2019-02-10 [1] CRAN (R 3.5.1)
# SummarizedExperiment             1.12.0    2018-10-30 [1] Bioconductor  
# survival                       * 2.44-1.1  2019-04-01 [1] CRAN (R 3.5.1)
# survminer                        0.4.4     2019-05-21 [1] CRAN (R 3.5.1)
# survMisc                         0.5.5     2018-07-05 [1] CRAN (R 3.5.1)
# sva                              3.30.1    2019-01-04 [1] Bioconductor  
# TCGAbiolinks                   * 2.10.5    2019-03-20 [1] Bioconductor  
# TH.data                        * 1.0-10    2019-01-21 [1] CRAN (R 3.5.1)
# tibble                           2.1.3     2019-06-06 [1] CRAN (R 3.5.1)
# tidyr                            0.8.3     2019-03-01 [1] CRAN (R 3.5.1)
# tidyselect                       0.2.5     2018-10-11 [1] CRAN (R 3.5.1)
# topGO                          * 2.34.0    2018-10-30 [1] Bioconductor  
# triebeard                        0.3.0     2016-08-04 [1] CRAN (R 3.5.1)
# tweenr                           1.0.1     2018-12-14 [1] CRAN (R 3.5.1)
# UpSetR                           1.4.0     2019-05-22 [1] CRAN (R 3.5.1)
# urltools                         1.7.3     2019-04-14 [1] CRAN (R 3.5.1)
# usethis                        * 1.5.0     2019-04-07 [1] CRAN (R 3.5.1)
# viridis                          0.5.1     2018-03-29 [1] CRAN (R 3.5.1)
# viridisLite                      0.3.0     2018-02-01 [1] CRAN (R 3.5.1)
# wheatmap                         0.1.0     2018-03-15 [1] CRAN (R 3.5.1)
# withr                            2.1.2     2018-03-15 [1] CRAN (R 3.5.1)
# xfun                             0.7       2019-05-14 [1] CRAN (R 3.5.1)
# XML                            * 3.98-1.20 2019-06-06 [1] CRAN (R 3.5.1)
# xml2                             1.2.0     2018-01-24 [1] CRAN (R 3.5.1)
# xtable                           1.8-4     2019-04-21 [1] CRAN (R 3.5.1)
# XVector                        * 0.22.0    2018-10-30 [1] Bioconductor  
# yaml                             2.2.0     2018-07-25 [1] CRAN (R 3.5.1)
# zip                              2.0.2     2019-05-13 [1] CRAN (R 3.5.1)
# zlibbioc                         1.28.0    2018-10-30 [1] Bioconductor  
# zoo                              1.8-6     2019-05-28 [1] CRAN (R 3.5.1)
# 
# [1] /usr/local/lib/R/site-library
# [2] /usr/local/lib/R/library

