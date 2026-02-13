rm(list = ls())

# Setup ####
path <- .libPaths()
newpath <- "C:/Users/albze08/Desktop/R/mQTL"
.libPaths(newpath)


pack_R <- c("dplyr","ggplot2","ggrepel","umap", "edgeR",
            "RColorBrewer", "pheatmap", "tidyverse",
            "igraph", "biomaRt", 
            "org.Hs.eg.db", "ggpubr", 
            "clusterProfiler",
            "igraph", "readxl")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

setwd("C:/Users/albze08/Desktop/postDoc/PBMC/")

set.seed(1)


# Load ####
source('C:/Users/albze08/Desktop/postDoc/functions/my_fisher_test.R')

# Annotation of individuals
anno <- read.csv("data/Wellness_barcodes.txt", sep="\t", header=T)
anno <- anno[anno$Sample.type=="Helblod",]

#Proteomics
protein_data <- read.csv("data/wellness_norm_final_794_wj_alternative_anno.txt", sep="\t", header=T)
protein_data$subject_id <- anno$Subject[match(protein_data$sample, anno$Wellness.id)]
protein_data$id <- paste0(protein_data$subject_id, ":", protein_data$visit)
rownames(protein_data) <- protein_data$id
protein_data <- subset(protein_data, select=-c(iid, sample, visit, id, subject_id))

#HPA annotation
hpa.cytokine <- read.table("HPA/cytokines/cytokine.tsv", sep="\t", header=T)
hpa.cytokine <- hpa.cytokine[grep("secreted", hpa.cytokine$Protein.class), ]
hpa <- read.table("HPA/all/proteinatlas.tsv", sep="\t", header=T)
cytokine <- protein_data[, colnames(protein_data) %in% hpa.cytokine$Gene]

# metadata
metadata <- read.csv("data/S3_Wellness_visit1_6_Clin_200109.txt", sep="\t", header=T)
metadata$visit <- gsub("Visit ", "", metadata$VisitName)
metadata$VisitName <- NULL
#metadata <- metadata[- which(metadata$Study=="T2D"), ]
#metadata$visit <- as.character(metadata$VisitName)
metadata$subject_id <- gsub("1-", "", metadata$subject)
metadata$id <- paste0(metadata$subject_id, ":", metadata$visit)
rownames(metadata) <- metadata$id

# Filter metadata
clinical <- metadata
rownames(clinical) <- metadata$id
clinical$Gender <- ifelse(clinical$Gender=="f", 0, 1)
clinical <- subset(clinical, select=c(Gender, Height, Cap_Gluc, SBP, DBP, Weight, Waist, Hip, Bioimp_fat, Bioimp_muscle, Bioimp_bone, PhysicalActivity, Age_at_Visit, BMI, 
                                      ProBNP, TNT, CRP, ALAT, GGT, HDL, Chol, LDL, TG, HbA1c, Urate, Gluc, CystC, Crea, ApoB.apoA1, ApoB, ApoA1))


#RNA-seq S3WP
rna_s3wp <- read.table("data/wellness_PBMC_v16_norm.txt", sep="\t", header = T)

#Gather lifestyle
lifestyle <- read.table("C:/Users/albze08/Desktop/postDoc/genome-metabolome/lifestyle/wellness_lifestyle.txt", header=T)
str <- strsplit(lifestyle$iid, "_") %>% unlist()
lifestyle$wellness.id <- str[seq(1,length(str),2)]
lifestyle$visit <- str[seq(2,length(str),2)] %>% as.character()
lifestyle$subject <- anno$Subject[match(lifestyle$wellness.id, anno$Wellness.id)]
lifestyle$id <- paste0(lifestyle$subject, ":", lifestyle$visit)
rownames(lifestyle) <- lifestyle$id
lifestyle <- subset(lifestyle, select=-c(iid, wellness.id, subject, visit, id)) %>% na.omit()

#auto-antibody
autoanti <- read.table("C:/Users/albze08/Desktop/postDoc/genome-protein/data/WELLNESS/Wellness/Data/rawdata/original.autoantibody.score.txt",
                       header=T)

rownames(autoanti) <- paste0(autoanti$subject, ":", autoanti$visit)
autoanti.plate <- autoanti$assay_plate
autoanti <- subset(autoanti, select=-c(SampleID,subject,barcode,visit,assay_plate))

#Cytof S3WP
cytof <- read.table("data/original.cytof.txt", sep="\t", header = T)
rownames(cytof) <- cytof$SampleID
cytof <- subset(cytof, select=-SampleID)
rownames(cytof) <- gsub("_", ":", rownames(cytof))

#Annotation CyTOF files
# macro.anno <- read.csv("data/cell_pop_codes.txt", sep="\t")
# macro.anno.my <- macro.anno[match(colnames(cytof), macro.anno$Populations),]
# macro.anno.my$Canonical_Pop[grep("CD4pos", macro.anno.my$Populations)] <- "CD4pos"
# macro.anno.my$Canonical_Pop[grep("CD8pos", macro.anno.my$Populations)] <- "CD8pos"
# macro.anno.my$Canonical_Pop[grep("CD56pos", macro.anno.my$Populations)] <- "CD56pos"
# macro.anno.my$Canonical_Pop[grep("Naive_B_cells", macro.anno.my$Populations)] <- "Naive_B_cells"
# macro.anno.my$Canonical_Pop[grep("CD4pos_naive", macro.anno.my$Populations)] <- "CD4pos_naive"
# macro.anno.my$Canonical_Pop[grep("CD8pos_naive", macro.anno.my$Populations)] <- "CD8pos_naive"
# rownames(macro.anno.my) <- macro.anno.my$Populations
# macro.anno.my <- subset(macro.anno.my, select=-Populations)
# colnames(macro.anno.my) <- c("family")

# Re-format RNA-seq ####
gene <- unique(rna_s3wp$ensg_id)
Ngene <- length(gene)
ind <- unique(rna_s3wp$subject)
Nind <- length(ind)

# Reformat Gene expression
sample <- paste0(rna_s3wp$subject, ":", rna_s3wp$visit)
rna_s3wp$sample <- sample

rna.long <- rna_s3wp[, c("sample", "ensg_id", "NX")]
rna <- pivot_wider(rna.long, names_from = ensg_id, values_from = NX) %>% as.data.frame()
rownames(rna) <- rna$sample
rna <- rna[,-1]
rownames(rna) <- gsub("1-", "", rownames(rna))
rownames(rna) <- gsub("V", "", rownames(rna))

rm(rna.long)

#remove genes with zero variance
v <- apply(rna, 2, sd)
rna <- rna[, v>0]

#remove lowly expressed genes
filter.genes <- filterByExpr(t(rna), min.count = 5, min.total.count = 10, large.n = 10, min.prop = 0.7)
rna <- rna[,filter.genes]

# symbol ID
ensg <- colnames(rna)
ensg.symbol <- data.frame(ensg=ensg , symbol=mapIds(org.Hs.eg.db, keys = ensg, keytype = "ENSEMBL", column="SYMBOL"))

colnames(rna) <- ensg.symbol$symbol[match(colnames(rna), ensg.symbol$ensg)]
rna <- rna[,!is.na(colnames(rna))]
ensg.symbol <- ensg.symbol %>% na.omit()

#log
rna.log <- log2(rna+1)


# Group populations into categories ####

# Natural Killer (NK) Cells
nk_cells <- c(
  "CD56neg_CD57neg_HLADRneg", "CD56neg_CD57neg_HLADRpos", "CD56neg_CD57pos_HLADRneg",
  "CD56neg_CD57pos_HLADRpos", "CD56_bright", "CD56dim_CD57neg_HLADRneg",
  "CD56dim_CD57neg_HLADRpos", "CD56dim_CD57pos_HLADRneg", "CD56dim_CD57pos_HLADRpos"
)

# Monocytes
monocytes_classical <- c("CD16neg_Classical_monocytes")
monocytes_intermediate <- c("Intermediate_Monocytes")
monocytes_nonclassical <- c("Non_classical_monocytes")

# Dendritic Cells (DCs)
mDC <- "mDC"
pDC <- "pDC"

# Basophils
basophils <- c("basophils")

# Innate Lymphoid Cells (ILCs)
ilcs <- c("Non_NK_ILC")

# Naive B Cells
naive_b_cells <- c("Naive_B_cells")

# B Cells
b_cells <- c(
  "Unswitched_memory_B_cells",
  "CD27neg_switched_memory_B_cells_CD22neg_CD24neg", "CD27neg_switched_memory_B_cells_CD22neg_CD24pos",
  "CD27neg_switched_memory_B_cells_CD22pos_CD24neg", "CD27neg_switched_memory_B_cells_CD22pos_CD24pos",
  "plasmablasts", "CD27pos_switched_memory_B_cells_CD22neg_CD24neg",
  "CD27pos_switched_memory_B_cells_CD22neg_CD24pos", "CD27pos_switched_memory_B_cells_CD22pos_CD24neg",
  "CD27pos_switched_memory_B_cells_CD22pos_CD24pos"
)

# Naive T Cells
naive_cd4 <- "CD4pos_naive"
naive_cd8 <- "CD8pos_naive"

# Central Memory (T_CM)
central_memory_cd4 <- c(
  "CD4pos_central_memory_CD57pos_CD39pos", "CD4pos_central_memory_CD57pos_CD39neg",
  "CD4pos_central_memory_CD57neg_CD39neg", "CD4pos_central_memory_CD57neg_CD39pos")

central_memory_cd8 <- c(
  "CD8pos_Central_memory_CD57neg_CD39neg", "CD8pos_Central_memory_CD57neg_CD39pos",
  "CD8pos_Central_memory_CD57pos_CD39neg", "CD8pos_Central_memory_CD57pos_CD39pos"
)

# Effector Memory (T_EM)
effector_memory_cd4 <- c(
  "CD4pos_Effector_memory_CD57pos_CD39neg", "CD4pos_Effector_memory_CD57neg_CD39pos",
  "CD4pos_Effector_memory_CD57neg_CD39neg", "CD4pos_Effector_memory_CD57pos_CD39pos")

effector_memory_cd8 <- c(
  "CD8pos_Effector_memory_CD57pos_CD39pos", "CD8pos_Effector_memory_CD57neg_CD39neg",
  "CD8pos_Effector_memory_CD57neg_CD39pos", "CD8pos_Effector_memory_CD57pos_CD39neg"
)

# TEMRA
temra_cd4 <- c(
  "CD4pos_TEMRA_CD57pos_CD39neg", "CD4pos_TEMRA_CD57neg_CD39pos",
  "CD4pos_TEMRA_CD57neg_CD39neg", "CD4pos_TEMRA_CD57pos_CD39pos"
)

temra_cd8 <- c(
  "CD8pos_TEMRA_CD57pos_CD39pos", "CD8pos_TEMRA_CD57pos_CD39neg",
  "CD8pos_TEMRA_CD57neg_CD39pos", "CD8pos_TEMRA_CD57neg_CD39neg"
)

# Combine into a list
pop.group <- list(
  NK_Cells = nk_cells,
  Monocytes_classical = monocytes_classical,
  Monocytes_intermediate = monocytes_intermediate,
  Monocytes_nonclassical = monocytes_nonclassical,
  Plasmacytoid_DC = pDC,
  Myeloid_DC = mDC,
  Basophils = basophils,
  ILCs = ilcs,
  Naive_B_Cells = naive_b_cells,
  B_Cells = b_cells,
  Naive_CD4 = naive_cd4,
  Naive_CD8 = naive_cd8,
  Central_Memory_CD4 = central_memory_cd4,
  Central_Memory_CD8 = central_memory_cd8,
  Effector_Memory_CD4 = effector_memory_cd4,
  Effector_Memory_CD8 = effector_memory_cd8,
  TEMRA_CD4 = temra_cd4,
  TEMRA_CD8 = temra_cd8
)
pop.ord <- names(pop.group)


cytof.nonnegative <- cytof
cytof.nonnegative[cytof.nonnegative<0] <- 2.2e-16
for (n in 1:nrow(cytof.nonnegative)){
  cytof.nonnegative[n,] <- cytof.nonnegative[n,] / sum(cytof.nonnegative[n,])
}

cytof.group <- matrix(0, nrow(cytof), length(pop.group)) %>% as.data.frame()
rownames(cytof.group) <- rownames(cytof)
colnames(cytof.group) <- names(pop.group)
for (n in 1:ncol(cytof.group)){
  if (length(pop.group[[n]])==1){
    cytof.group[,n] <- cytof.nonnegative[,pop.group[[n]]]
  } else {
    cytof.group[,n] <- rowSums(cytof.nonnegative[,pop.group[[n]]])
  }
}

macro.anno.my <- setNames(rep("", ncol(cytof)), colnames(cytof))
for (pop in names(pop.group)){
  macro.anno.my[names(macro.anno.my) %in% pop.group[[pop]]] <- pop
}

# Extract GWAS results of specific immune families ####
pval_loose <- 1e-6

#calculate number independent pop, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6693916/
cov.cytof <- cov(cytof.nonnegative)
eigen.cytof <- eigen(cov.cytof, symmetric = T)$values
indep.pop <- (sum(eigen.cytof)^2)/sum(eigen.cytof^2)
pval.cytof <- 5e-8/indep.pop

# gather cQTLs
folder_cQTL <- "cQTL-mean/"
var <- grep("linear", list.files(paste0(folder_cQTL, "conditional")), value=T)
var <-  gsub("\\..*", "", var)
cQTL_files <- grep(paste(var, collapse="|"), list.files(folder_cQTL), value=T)
cQTL_files <- cQTL_files[grep("linear", cQTL_files)]
df_summary <- data.frame()
for (n in 1:length(cQTL_files)){
  f <- read.table(paste0(folder_cQTL, cQTL_files[n]), header=F) 
  colnames(f) <- c("CHR", "BP", "SNP", paste0("V",1:11), "P", "V12")
  f <- f[f$P<pval_loose,]
  if (nrow(f)>0){
    f$group <- gsub("out\\.", "", cQTL_files[n])
    f$group <- gsub("\\.glm\\.linear", "", f$group)
    
    df_summary <- rbind(df_summary, f)
  }
}
df_summary$family <- macro.anno.my[match(df_summary$group, rownames(macro.anno.my))]
GWAS.cytof <- df_summary
GWAS.cytof$id <- paste0(GWAS.cytof$CHR, ":", GWAS.cytof$BP)
GWAS.cytof$BP <- as.numeric(GWAS.cytof$BP)

#LD
folder_cQTL <- "cQTL-mean/"
cQTL_files <- grep("clumps", list.files(paste0(folder_cQTL, "LD")), value=T)
df_summary.cqtl.clump <- data.frame()
for (n in 1:length(cQTL_files)){
  var <- gsub("LD\\.out\\.", "", cQTL_files[n])
  var <- gsub("\\.glm\\.linear\\.clumps", "", var)
  
  f <- read.table(paste0(folder_cQTL, "LD/", cQTL_files[n]), header=F)
  colnames(f) <- c("CHR", "BP", "SNP", "P", paste0("V",1:6))
  f <- f[f$P<pval.cytof,]
  if (nrow(f)>0){
    f$group <- var
    df_summary.cqtl.clump <- rbind(df_summary.cqtl.clump, f[which.min(f$P),c("SNP", "P", "group")])
  }
  if (nrow(f)>1){
    conditional_file <- grep(paste0(var, "\\."), list.files(paste0(folder_cQTL, "conditional")), value=T)
    conditional_file <- grep("linear", conditional_file, value=T)
    if (length(conditional_file)==1){
      f_conditional <- read.table(paste0(folder_cQTL, "conditional/", conditional_file), header=F)
      colnames(f_conditional) <- c("CHR", "BP", "SNP", paste0("V", 1:11), "P", "ERRCODE")
      f_conditional$group <- var
      
      idx <- which(f_conditional$P<0.1) #this is the conditional p-value, i.e. when doing linear model that includes the dominant variant
      f_conditional <- f_conditional[idx, ]
      f_conditional$P <- f$P[match(f_conditional$SNP, f$SNP)] #store the actual p-value from GWAS, not the conditional
      df_summary.cqtl.clump <- rbind(df_summary.cqtl.clump, f_conditional[,c("SNP", "P", "group")])
    } else {
      print(paste0("[WARNING] Number of files found for ", var, " : ", length(conditional_file)))
    }
  }
}
GWAS.cytof.LD <- df_summary.cqtl.clump
GWAS.cytof.LD$family <- macro.anno.my[match(GWAS.cytof.LD$group, rownames(macro.anno.my))]

str <- strsplit(GWAS.cytof.LD$SNP, "\\:") %>% unlist()
GWAS.cytof.LD$CHR <- str[seq(1,length(str),4)] %>% as.numeric()
GWAS.cytof.LD$BP <- str[seq(2,length(str),4)] %>% as.numeric()
GWAS.cytof.LD$id <- paste0(GWAS.cytof.LD$CHR, ":", GWAS.cytof.LD$BP)

#saveRDS(GWAS.cytof.LD, "GWAS.cytof.LD.RDS")


# Extract GWAS results of RNA-seq ####

#cov.rna <- cov(rna.log)
#eigen.rna <- eigen(cov.rna)$values
#indep.gene <- (sum(eigen.rna)^2)/sum(eigen.rna^2)
indep.gene <- 15.96234 #value obtained with 6354 genes after filtering
pval.rna <- 5e-8/indep.gene

# Gather eQTLs
folder_eQTL <- "eQTL-mean/"
var <- grep("linear", list.files(paste0(folder_eQTL, "conditional")), value=T)
var <-  gsub("\\..*", "", var)
eQTL_files <- grep(paste(var, collapse="|"), list.files(folder_eQTL), value=T)
eQTL_files <- eQTL_files[grep("linear", eQTL_files)]
df_summary <- data.frame()
for (n in 1:length(eQTL_files)){
  if (file.info(paste0(folder_eQTL, eQTL_files[n]))$size>0){
    f <- read.table(paste0(folder_eQTL, eQTL_files[n]), header=F)
    colnames(f) <- c("CHR", "BP", "SNP", paste0("V",1:11), "P", "V12")
    f <- f[f$P<pval_loose,]
    if (nrow(f)>0){
      f$group <- gsub("out\\.", "", eQTL_files[n])
      f$group <- gsub("\\.glm\\.linear", "", f$group)
      
      df_summary <- rbind(df_summary, f)
    }
  }
}
GWAS.rna <- df_summary
GWAS.rna$symbol <- GWAS.rna$group

#add symbol
# ensg <- GWAS.rna$group
# ensg.symbol <- data.frame(ensg=ensg , symbol=mapIds(org.Hs.eg.db, keys = ensg, keytype = "ENSEMBL", column="SYMBOL"))
# GWAS.rna$symbol <- ensg.symbol$symbol[match(GWAS.rna$group, ensg.symbol$ensg)]

#LD
folder_eQTL <- "eQTL-mean/"
eQTL_files <- grep("clumps", list.files(paste0(folder_eQTL, "LD")), value=T)
df_summary.eQTL.clump <- data.frame()
for (n in 1:length(eQTL_files)){
  var <- gsub("LD\\.out\\.", "", eQTL_files[n])
  var <- gsub("\\.glm\\.linear\\.clumps", "", var)
  
  f <- read.table(paste0(folder_eQTL, "LD/", eQTL_files[n]), header=F)
  colnames(f) <- c("CHR", "BP", "SNP", "P", paste0("V",1:6))
  f <- f[f$P<pval.rna,]
  if (nrow(f)>0){
    f$group <- var
    df_summary.eQTL.clump <- rbind(df_summary.eQTL.clump, f[which.min(f$P),c("SNP", "P", "group")])
  }
  if (nrow(f)>1){
    conditional_file <- grep(paste0(var, "\\."), list.files(paste0(folder_eQTL, "conditional")), value=T)
    conditional_file <- grep("linear", conditional_file, value=T)
    if (length(conditional_file)==1){
      f_conditional <- read.table(paste0(folder_eQTL, "conditional/", conditional_file), header=F)
      colnames(f_conditional) <- c("CHR", "BP", "SNP", paste0("V", 1:11), "P", "ERRCODE")
      f_conditional$group <- var
      
      idx <- which(f_conditional$P<0.1) #this is the conditional p-value, i.e. when doing linear model that includes the dominant variant
      f_conditional <- f_conditional[idx, ]
      f_conditional$P <- f$P[match(f_conditional$SNP, f$SNP)] #store the actual p-value from GWAS, not the conditional
      df_summary.eQTL.clump <- rbind(df_summary.eQTL.clump, f_conditional[,c("SNP", "P", "group")])
    } else {
      print(paste0("[WARNING] Number of files found for ", var, " : ", length(conditional_file)))
    }
  }
}
GWAS.rna.LD <- df_summary.eQTL.clump
GWAS.rna.LD$symbol <- GWAS.rna.LD$group
#saveRDS(GWAS.rna.LD, "GWAS.rna.LD.RDS")

# Common genetics between PBMC frequencies and gene expression ####
cell.gene <- readRDS("Ensemble_ILR-306025.RDS") %>% bind_rows() %>% mutate(pval.bon=pval*n()) %>% filter(pval.bon<0.05) %>% filter(gene %in% colnames(rna))
colnames(cell.gene)[1] <- "pop"
shared_snp <- intersect(GWAS.cytof.LD %>% pull(SNP), 
                        GWAS.rna %>% filter(P<5e-8/indep.gene) %>% pull(SNP)) 
df <- data.frame(snp=shared_snp, gene="", population="", family="")
for (n in 1:nrow(df)){
  g <- GWAS.rna %>% filter(P<5e-8/indep.gene) %>% filter(SNP == df$snp[n]) %>% pull(group) %>% unique()
  df$gene[n] <- g %>% paste(collapse="; ")
  df$family[n] <- GWAS.cytof %>% filter(P<5e-8/indep.pop) %>% filter(SNP == df$snp[n]) %>% pull(family) %>% unique() %>% paste(collapse="; ")
  df$population[n] <- GWAS.cytof %>% filter(P<5e-8/indep.pop) %>% filter(SNP == df$snp[n]) %>% pull(group) %>% unique() %>% paste(collapse="; ")
  df$population.network[n] <- cell.gene %>% filter(gene %in% g) %>% pull(pop) %>% unique() %>% paste(collapse="; ")
}
df$gene %>% paste(collapse="; ") %>% strsplit("; ") %>% pluck(1) %>% unique()


# Put GWAS of cytof, rna and cytokine in the same plot ####
MAX_SNP <- 10e3

#Prepare GWAS CyTOF
df_plot_cytof <- GWAS.cytof[order(GWAS.cytof$P, decreasing=F), c("CHR", "BP", "P", "group", "family", "SNP")]
df_plot_cytof <- df_plot_cytof[!df_plot_cytof$CHR %in% c("PAR1", "PAR2"),]
df_plot_cytof$CHR <- as.numeric(df_plot_cytof$CHR)

df_plot_cytof$symbol <- df_plot_cytof$group
df_plot_cytof$pop <- df_plot_cytof$group
df_plot_cytof <- df_plot_cytof[, c("CHR", "BP", "P", "pop", "family", "symbol", "SNP")]
N_cqtl <- nrow(df_plot_cytof)

# prepare GWAS RNA-seq
df_plot_rna <- GWAS.rna[order(GWAS.rna$P, decreasing=F)[1:MAX_SNP], c("CHR", "BP", "P", "symbol", "SNP")]
df_plot_rna <- df_plot_rna[!df_plot_rna$CHR %in% c("PAR1", "PAR2"),]
df_plot_rna$CHR <- as.numeric(df_plot_rna$CHR)

df_plot_rna$pop <- ""
df_plot_rna$family <- ""
df_plot_rna <- df_plot_rna[, c("CHR", "BP", "P", "pop", "family", "symbol", "SNP")]
N_eqtl <- nrow(df_plot_rna)

#merge
df_plot <- rbind(df_plot_cytof, df_plot_rna)
idx.cytof <- 1:N_cqtl
idx.rna <- (N_cqtl+1):(N_cqtl+N_eqtl)

chr.df <- data.frame(CHR=1:22)
chr.df$BP_max <- 0
chr.df$BP_mean <- 0
df_plot$BP_cum <- 0
max.previous.chr <- 0

stopifnot(is.numeric(df_plot$CHR))
for (n in 1:22){ #chromosomes
  idx <- which(df_plot$CHR==n)
  
  #new BP in plot df
  df_plot$BP_cum[idx] <- df_plot$BP[idx] + max.previous.chr 
  
  #chromosome charac
  chr.df$BP_max[n] <- max(df_plot$BP_cum[idx])
  chr.df$BP_mean[n] <- mean(df_plot$BP_cum[idx])
  
  #update baseline
  max.previous.chr <- max(df_plot$BP_cum[idx])
}
df_plot$color <- ifelse(df_plot$CHR %% 2 ==0, "even", "odd")
df_plot_cytof <- df_plot[idx.cytof, c("CHR", "BP", "P", "pop", "family", "symbol", "BP_cum", "color")]
df_plot_rna <- df_plot[idx.rna, c("CHR", "BP", "P", "symbol", "SNP", "BP_cum", "color")]


# Plot CyTOF, RNA-seq and cytokine GWAS ####
pop.palette <- readRDS("pop.palette.RDS")
group.color <- pop.palette$color
names(group.color) <- pop.palette$family
df_plot_cytof$symbol <- ifelse(df_plot_cytof$P<pval.cytof, df_plot_cytof$symbol, "")

df_plot_cytof$color <- ifelse(df_plot_cytof$P<pval.cytof, df_plot_cytof$family, df_plot_cytof$color)
p1 <- ggplot(df_plot_cytof, aes(x=BP_cum, y=-log10(P), color=color)) +
  geom_point() + 
  geom_hline(yintercept=-log10(pval.cytof), linetype=2) +
  theme_bw() +
  geom_text_repel(aes(label=symbol), cex=2, color="black") +
  scale_x_continuous(label = chr.df$CHR, breaks = chr.df$BP_mean) +
  scale_color_manual(values=c(group.color, even="gray77", odd="gray23")) + 
  scale_size_continuous(range = c(0.5,3)) + 
  labs(x = NULL, y = "-log<sub>10</sub>(p)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0, size=8, vjust=0.5),
        legend.position="none",
        text=element_text(size=12, family="Arial")) +
  ylab("-Log10 P-value")
#p1

df_plot_rna$symbol <- ifelse(df_plot_rna$P<pval.rna, df_plot_rna$symbol, "")
p2 <- ggplot(df_plot_rna, aes(x=BP_cum, y=-log10(P), color=color)) +
  geom_point() +
  geom_hline(yintercept=-log10(pval.rna), linetype=2) +
  scale_color_manual(values=c(even="gray77", odd="gray23")) + 
  geom_text_repel(aes(label=symbol), cex=2, color="black") +
  theme_bw() +
  scale_x_continuous(label = chr.df$CHR, breaks = chr.df$BP_mean) +
  scale_size_continuous(range = c(0.5,3)) + 
  labs(x = NULL, y = "-log<sub>10</sub>(p)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0, size=8, vjust=0.5),
        legend.position="none",
        text=element_text(size=12, family="Arial")) +
  ylab("-Log10 P-value")
#p2

p <- ggarrange(p1,p2, ncol=1)
p
cairo_pdf("GWAS-cytof-rna.pdf", width=6, height=6.6)
print(p)
dev.off()


# Some pies about SNPs ####
group.color <- pop.palette$color
names(group.color) <- pop.palette$family
df.plot <- table(GWAS.cytof.LD$family) %>% as.data.frame()
df.plot$Prc <- df.plot$Freq/sum(df.plot$Freq)

df.plot$Var1 <- factor(df.plot$Var1, levels=df.plot$Var1[order(df.plot$Freq, decreasing=F)])
p <- ggplot(df.plot, aes(x="", y=Freq, fill=Var1)) + 
  geom_bar(stat="identity", width=1, color="white") + 
  geom_text(aes(label = paste0(round(Prc*100), "% ", "n=", Freq)), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values=group.color) +
  coord_polar("y", start=0) + 
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        size=16, family="Arial") +
  ggtitle("Immune populations") +
  xlab("") + ylab("")
p
cairo_pdf("pie-populations.pdf", width=3, height=3)
print(p)
dev.off()


# Plot overlap between PBMC omics

# Venn between omics
library("ggvenn")

data <- list(
  ctQTL = GWAS.cytof.LD$SNP %>% unique(),
  eQTL = GWAS.rna.LD$SNP %>% unique(),
  pQTL = GWAS.cytokine.LD$SNP %>% unique()
)

pdf("QTL-venn.pdf", width=3, height=3)
ggvenn(data, c("ctQTL", "eQTL", "pQTL"))
dev.off()


# Load Genotype ####
library("vcfR")
vcf <- read.vcfR("C:/Users/albze08/Desktop/postDoc/wellness/Wellness/Data/WGS/wellness_101_subject_QC.vcf")
snp <- vcf@fix %>% as.data.frame()
geno <- vcf@gt %>% as.data.frame()
geno <- geno[,2:ncol(geno)]
colnames(geno) <- anno$Subject[match(colnames(geno), paste0("P", anno$Barcode))]


# Plot examples ####
df.plot <- data.frame(snp=intersect(GWAS.lmm.LD$SNP, GWAS.lmm.eqtl$SNP)) 
df.plot$pop <- GWAS.lmm.LD$group[match(df.plot$snp, GWAS.lmm.LD$SNP)]
df.plot$gene <- GWAS.lmm.eqtl$group[match(df.plot$snp, GWAS.lmm.eqtl$SNP)]

for (n in 1:nrow(df.plot)){
  idx <- match(df.plot$snp[n], snp$ID)
  
  df.plot.pop <- data.frame(sample=rownames(cytof), variable=cytof.nonnegative[,df.plot$pop[n]], omic="CyTOF")
  df.plot.pop$rank <- 100*(df.plot.pop$variable %>% rank())/nrow(df.plot.pop)
  df.plot.pop$ind <- gsub("\\:.*", "", df.plot.pop$sample)
  
  df.plot.pop$geno <- geno[idx, df.plot.pop$ind] %>% as.character()
  
  df.plot.gene <- data.frame(sample=rownames(rna.log), variable=rna.log[,df.plot$gene[n]], omic="RNA-seq")
  df.plot.gene$rank <- 100*(df.plot.gene$variable %>% rank())/nrow(df.plot.gene)
  df.plot.gene$ind <- gsub("\\:.*", "", df.plot.gene$sample)
  
  df.plot.gene$geno <- geno[idx, df.plot.gene$ind] %>% as.character()
  
  df.plot.n <- rbind(df.plot.pop, df.plot.gene) %>% na.omit()
  
  print(paste0(df.plot$pop[n], "; ", df.plot$gene[n],"; ", df.plot$snp[n], ";"))
  print(table(df.plot.pop$geno))
  print(table(df.plot.gene$geno))
  
  p <- ggplot(df.plot.n, aes(x=geno, y=rank, fill=omic)) + 
    geom_boxplot() + 
    #geom_beeswarm(dodge.width=0.8,  alpha=0.6) + 
    theme_classic() + xlab("Genotype") + ylab("Percentile in the cohort") +
    ggtitle(paste0(df.plot$pop[n], "; ", df.plot$gene[n],"; ", df.plot$snp[n]))
  
cairo_pdf(paste0("pdfs/", df.plot$pop[n], "; ", df.plot$gene[n], ".pdf"), height=3, width=4)
  print(p)
  dev.off()
}

# Write GWAS Table ####
tab1 <- GWAS.cytof.LD[1:3]
colnames(tab1)[3] <- "population"

tab2 <- GWAS.rna.LD[,c(1,2,4)]
colnames(tab2)[3] <- "gene"

writexl::write_xlsx(list(tab1,tab2), "Supplementary Table 4_GWAS.xlsx")


#Load GWAS ####
GWAS.cytof.LD <- readRDS("GWAS.cytof.LD.RDS") %>% mutate(id=paste0(CHR, ":", BP))
GWAS.rna.LD <- readRDS("GWAS.rna.LD.RDS") %>% mutate(id=paste0(CHR, ":", BP))
GWAS.cytof.main.LD <- GWAS.cytof.LD
GWAS.cytof.main.LD$group <- GWAS.cytof.main.LD$family


# function to lift ####

lift_fun <- function(gwas.df){
  library("GenomicRanges")
  library("rtracklayer")
  chain <- import.chain("liftOver/hg19ToHg38.over.chain/hg19ToHg38.over.chain")
  
  gr <- GRanges(
    seqnames = paste0("chr", gwas.df$CHR),
    ranges = IRanges(start = gwas.df$BP, width = 1),
    strand = "*"
  )
  lifted <- liftOver(gr, chain)
  
  mapped <- elementNROWS(lifted) == 1  # Only 1:1 mappings
  gwas.df.mapped <- gwas.df[mapped, ]
  
  lifted_gr <- unlist(lifted[mapped])
  
  gwas.df.mapped$CHR <- gsub("chr", "", as.character(seqnames(lifted_gr)))
  gwas.df.mapped$BP <- start(lifted_gr)
  gwas.df.mapped$id <- paste0(gwas.df.mapped$CHR, ":", gwas.df.mapped$BP)
  
  return(gwas.df.mapped)
}

closest_snp_fun <- function(df1, df2){
  df.out <- df1[df1$CHR %in% df2$CHR,]
  df.out$distance <- NA
  df.out$CHR_other <- NA
  df.out$BP_other <- NA
  df.out$idx <- NA
  for (n in 1:nrow(df.out)){
    idx <- which(df2$CHR == df1$CHR[n])
    d <- abs(df1$BP[n] - df2$BP[idx])
    idx.min <- which.min(d)
    
    df.out$idx[n] <- idx[idx.min]
    df.out$distance[n] <- d[idx.min]
    df.out$CHR_other[n] <- df2$CHR[idx[idx.min]]
    df.out$BP_other[n] <- df2$BP[idx[idx.min]]
  }
  return(df.out %>% arrange(distance))
}

# Roederer ####
meta.roederer <- read_excel("meta-GWAS/Roederer/NIHMS668942-supplement-11.xlsx")
meta.roederer <- meta.roederer[2:nrow(meta.roederer),]
colnames(meta.roederer) <- meta.roederer[1,]
meta.roederer <- meta.roederer[2:nrow(meta.roederer),]

gwas.roederer <- read_excel("meta-GWAS/Roederer/NIHMS668942-supplement-13.xlsx")
gwas.roederer <- gwas.roederer[2:nrow(gwas.roederer),]
colnames(gwas.roederer) <- gwas.roederer[1,]
gwas.roederer <- gwas.roederer[2:nrow(gwas.roederer),]

colnames(gwas.roederer)[colnames(gwas.roederer)=="Chromosome"] <- "CHR"
colnames(gwas.roederer)[colnames(gwas.roederer)=="Position"] <- "BP"

gwas.roederer$`P Value` <- as.numeric(gwas.roederer$`P Value`)
gwas.roederer$BP <- as.numeric(gwas.roederer$BP)

gwas.roederer.lift <- lift_fun(gwas.roederer)


GWAS.summary.roederer <- closest_snp_fun(GWAS.cytof.LD, gwas.roederer.lift) %>%
  mutate(trait=gwas.roederer.lift$`Cell-type`[idx])

GWAS.summary.roederer$celltype <- meta.roederer$`Cell Subset`[match(GWAS.summary.roederer$trait, meta.roederer$`SPEL Trait ID`)]


# Astle ####
gwas.astle <- read_excel("meta-GWAS/Astle/mmc3.xlsx")
colnames(gwas.astle) <- gwas.astle[1,]
gwas.astle <- gwas.astle[2:nrow(gwas.astle),]

gwas.astle$CHR <- gwas.astle$`Chr (GRCh37)`
gwas.astle$BP <- gwas.astle$`BP (GRCh37)` %>% as.numeric()

gwas.astle.lift <- lift_fun(gwas.astle)

GWAS.summary.astle <- closest_snp_fun(GWAS.cytof.LD, gwas.astle.lift) %>% 
  mutate(trait=gwas.astle.lift$`Blood Index Classes Associated with Locus`[idx])


# Orru #### 
GWAS.orru <- read_excel("Orru/NIHMS1739370-supplement-Supp_table.xlsx", sheet="Supplementary Table 3A")
GWAS.orru <- GWAS.orru[!is.na(GWAS.orru[,2]), ]
colnames(GWAS.orru) <- GWAS.orru[1,]
GWAS.orru <- GWAS.orru[2:nrow(GWAS.orru),]

str <- strsplit(GWAS.orru$`Top variant (chr:pos in GRC37-hg19)`, ":") %>% unlist()
GWAS.orru$CHR <- str[seq(1,length(str),by=2)] 
GWAS.orru$BP <- str[seq(2,length(str),by=2)] %>% as.numeric()
GWAS.orru$id <- paste0(GWAS.orru$CHR, ":", GWAS.orru$BP)

# GWAS.orru.full <- data.frame()
# for (n in 1:nrow(GWAS.orru)){
#   id.vec <- strsplit(GWAS.orru$`Reported credible set list of SNPs`[n], ",")[[1]]
#   GWAS.orru.full <- rbind(GWAS.orru.full, 
#                           data.frame(Trait=GWAS.orru$Trait[n], id=id.vec))
# }
# str <- strsplit(GWAS.orru.full$id, ":") %>% unlist()
# GWAS.orru.full$CHR <- str[seq(1,length(str),by=2)] %>% as.numeric()
# GWAS.orru.full$BP <- str[seq(2,length(str),by=2)] %>% as.numeric()
# GWAS.orru <- GWAS.orru.full

GWAS.orru.lift <- lift_fun(GWAS.orru)

GWAS.summary.orru <- closest_snp_fun(GWAS.cytof.LD, GWAS.orru.lift) %>%
  mutate(trait=GWAS.orru.lift$`Statistical trait name`[idx])


# Ferreira et al., https://genepi.qimr.edu.au/staff/manuelf/gwas_results/main.html ####
gwas.folder <- list.dirs("meta-GWAS/QIMR") %>% setdiff("meta-GWAS/QIMR")
GWAS.summary.ferreira <- data.frame()
for (f in gwas.folder){
  file_name <- list.files(f)
  #load
  gwas <- read.table(paste0(f, "/", file_name), header=T) %>% filter(P<1e-3)
  gwas$id <- paste0(gwas$CHR, ":", gwas$BP)
  
  gwas.mapped <- lift_fun(gwas)
  
  closest_matches <- closest_snp_fun(GWAS.cytof.LD, gwas.mapped) %>% 
    mutate(trait=file_name) %>% 
    filter(distance<1e4)
  
  GWAS.summary.ferreira <- rbind(GWAS.summary.ferreira, closest_matches)
}
GWAS.summary.ferreira <- GWAS.summary.ferreira %>% arrange(distance)


#create SNP id of same format in my data
str <- strsplit(GWAS.cytof.LD$SNP, "\\:") %>% unlist()
GWAS.cytof.LD$CHR <- str[seq(1,length(str),4)]
GWAS.cytof.LD$BP <- str[seq(2,length(str),4)]
GWAS.cytof.LD$A1 <- str[seq(3,length(str),4)]
GWAS.cytof.LD$id <- paste0(GWAS.cytof.LD$CHR, ":", GWAS.cytof.LD$BP)

str <- strsplit(GWAS.cytof$SNP, "\\:") %>% unlist()
GWAS.cytof$CHR <- str[seq(1,length(str),4)]
GWAS.cytof$BP <- str[seq(2,length(str),4)]
GWAS.cytof$A1 <- str[seq(3,length(str),4)]
GWAS.cytof$id <- paste0(GWAS.cytof$CHR, ":", GWAS.cytof$BP)



# DICE ####
gwas.folder <- list.files("meta-GWAS/DICE/") 

DICE.summary <- data.frame()
for (f in gwas.folder){

  gwas <- read.table(paste0("meta-GWAS/DICE/", f), sep="\t") 
  colnames(gwas) <- c("CHR", "BP", "rs", "A1", "A2", "remove", "remove2", "str")
  str <- strsplit(gwas$str, ";") %>% unlist()
  gwas$ensg <- str[seq(1,length(str),by=4)]
  gwas$symbol <- str[seq(2,length(str),by=4)]
  gwas$pval <- str[seq(3,length(str),by=4)]
  gwas$beta <- str[seq(4,length(str),by=4)]
  
  gwas <- select(gwas, -remove, -remove2, -str)
  
  gwas$CHR <- gsub("chr", "", gwas$CHR) 
  gwas$pval <- gsub("Pvalue=", "", gwas$pval) %>% as.numeric()
  
  gwas <- gwas %>% filter(pval<5e-8)
  
  closest_matches <- GWAS.cytof.LD %>%
    rowwise() %>%
    mutate(
      # all SNPs on same chromosome in df2
      candidates = list(gwas$BP[gwas$CHR == CHR]),
      # only compute if there are candidates
      idx = if (length(candidates) > 0) which.min(abs(BP - candidates)) else NA_integer_,
      closest_BP = if (!is.na(idx)) candidates[idx] else NA_integer_,
      distance   = if (!is.na(idx)) abs(BP - closest_BP) else NA_integer_
    ) %>%
    ungroup() %>%
    select(-candidates, -idx) %>%
    select(id, group, closest_BP, distance) %>% 
    mutate(trait=gsub("\\.vcf", "", f)) %>% 
    filter(distance<1e3)
  
  DICE.summary <- rbind(DICE.summary, closest_matches)
  
}

# Specific traits from Orru ####
gwas <- read.table("meta-GWAS/ebi-a-GCST90001407.vcf/ebi-a-GCST90001407.vcf")
gwas$id <- paste0(gwas$V1, ":", gwas$V2)
str <- strsplit(gwas$V9, ":") %>% unlist()
gwas$ES <- str[seq(1,length(str), 5)]
gwas$SE <- str[seq(2,length(str), 5)]
gwas$LP <- str[seq(3,length(str), 5)]
gwas$AF <- str[seq(4,length(str), 5)]
gwas$ID <- str[seq(5,length(str), 5)]

gwas.ebi <- read.vcfR("meta-GWAS/ebi-a-GCST90001436.vcf/ebi-a-GCST90001436.vcf")
chrom <- getCHROM(gwas.ebi)
pos   <- getPOS(gwas.ebi)
snp.ebi <- data.frame(snp=paste0(chrom, ":", pos))
snp.ebi$logP <- extract.gt(gwas.ebi, element = "LP", as.numeric = TRUE)

# Download summary data for eQTL Meta-GWAS ####
tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

tissue_label <- c("B cell", "memory B cell",
                  "CD4+ T cell", "CD4+ CTL cell", "CD4+ memory T cell", "CD4+ TCM cell", "CD4+ TEM cell", 
                  "CD8+ T cell", "CD8+ TCM cell", "CD8+ TEM cell",
                  "CD16+ monocyte", "monocyte",
                  "CD56+ NK cell", "NK cell", 
                  "dendritic cell", "dnT cell", "plasmacytoid dendritic cell")

for (tis in tissue_label){
  print(tis)
  idx <- which(tabix_paths$tissue_label==tis & tabix_paths$quant_method=="ge" & tabix_paths$condition_label=="naive")
  df <- tabix_paths[idx,]
  for (k in 1:nrow(df)){
    file_name <- df$ftp_path[k] 
    file_name <- paste0(sub("/[^/]*$", "", file_name), "/", df$dataset_id[k], ".permuted.tsv.gz")
    
    dest_file <- paste0(tis, "-", df$sample_group[k])
    download.file(file_name,
                  destfile = paste0("C:/Users/albze08/Desktop/postDoc/wellness-rna/meta-GWAS/eQTL_catalog/", dest_file, ".gz"),
                  mode="wb")
  }
}
# download.file("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000499/QTD000499.permuted.tsv.gz",
#               destfile = paste0("C:/Users/albze08/Desktop/postDoc/wellness-rna/meta-GWAS/eQTL_catalog/", dest_file, ".gz"),
#               mode="wb")


# eQTL Meta-GWAS on memory B cells, CD8 effector, NK ####
gwas.folder <- list.dirs("meta-GWAS/eQTL_catalog") %>% setdiff("meta-GWAS/eQTL_catalog")
gwas.list <- list()
df.out <- data.frame()
for (f in gwas.folder){
  file_name <- list.files(f)
  gwas <- read.table(paste0(f, "/", file_name), header=T) %>% filter(pvalue<5e-8)
  #SNP id
  str <- strsplit(gwas$variant, "_") %>% unlist()
  gwas$snp <- paste0(str[seq(1,length(str), by=4)], ":",
                     str[seq(2,length(str), by=4)], ":",
                     str[seq(3,length(str), by=4)], ":",
                     str[seq(4,length(str), by=4)])
  gwas$snp <- gsub("chr", "", gwas$snp)
  
  gwas$symbol <- ensg.symbol$symbol[match(gwas$molecular_trait_object_id, ensg.symbol$ensg)]
  gwas <- na.omit(gwas)
  
  #store
  idx <- which(gwas$snp %in% GWAS.rna$SNP)
  gwas.list[[file_name]] <- intersect(GWAS.rna$SNP, gwas$snp)
  
  df <- gwas[idx, c("snp", "symbol")]
  df.out <- rbind(df.out, data.frame(df, celltype=file_name))
}
  
df.out$symbol_my <- GWAS.rna$group[match(df.out$snp, GWAS.rna$SNP)]


# cell.gene <- readRDS("Ensemble_ILR-306025.RDS") %>% bind_rows() %>% 
#   mutate(pval.bon=pval*n()) %>% filter(pval.bon<0.05) %>% filter(gene %in% colnames(rna))
cell.gene <- readRDS("Ensemble_ILR-306025.RDS") %>% bind_rows() %>% 
  mutate(pval.BH=p.adjust(pval, method="BH")) %>% filter(pval.BH<0.05) %>% filter(gene %in% colnames(rna))
colnames(cell.gene)[1] <- "pop"
df.out$celltype_my <- lapply(df.out$symbol_my, function(x){cell.gene %>% filter(gene==x) %>% pull(pop) %>% paste(collapse=";")}) %>% unlist()

  
  

# Big linear model ####
var.group <- data.frame(var=c("subject_id",
                              "Age_at_Visit", 
                              "Gender", 
                              "Height", "Weight", "Waist", "Hip", "BMI", "Bioimp_fat", "Bioimp_muscle", "Bioimp_bone", 
                              "SBP", "DBP", 
                              "Cap_Gluc", "Gluc", "HbA1c", 
                              "Chol", "LDL", "HDL", "TG", "ApoB", "ApoA1", "ApoB.apoA1",   
                              "CRP", 
                              "ALAT", "GGT", 
                              "ProBNP", "TNT", 
                              "Urate", 
                              "CystC", "Crea",
                              "Smoking", 
                              "MaritalStatus_change", "Housing_change", "Employment_change",
                              "Stress",
                              "PhysicalActivity", "SedentaryTime_hours",
                              "PerceivedHealth",
                              "Common_cold_influenzae",
                              "NSAID_painmed", "Bp_med", "Lipid_med", "Antibiotics_med"),
                        group=c("subject_id",
                                "Age", 
                                "Sex", 
                                "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", 
                                "Blood pressure", "Blood pressure", 
                                "Glucose homeostasis", "Glucose homeostasis", "Glucose homeostasis", 
                                "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile",   
                                "CRP", 
                                "Organ biomarker", "Organ biomarker", 
                                "Organ biomarker", "Organ biomarker", 
                                "Organ biomarker", 
                                "Organ biomarker", "Organ biomarker",
                                "Lifestyle", 
                                "Lifestyle", "Lifestyle", "Lifestyle",
                                "Lifestyle",
                                "Lifestyle", "Lifestyle",
                                "Lifestyle",
                                "Lifestyle",
                                "Lifestyle", "Lifestyle", "Lifestyle", "Lifestyle"))


var.group.SNP <- data.frame(var=unique(c(GWAS.cytof.main.LD$SNP, GWAS.cytof.LD$SNP, GWAS.rna.LD$SNP)),
                            group="Genetics")
var.group.SNP$var <- paste0("SNP:", var.group.SNP$var)
var.group.SNP$var <- gsub(":", "\\.", var.group.SNP$var)

var.group <- rbind(var.group, var.group.SNP)

lm_fun <- function(X,Y, GWAS.summary){
  common.samples <- intersect(rownames(X), rownames(Y))
  stopifnot(length(common.samples)>0)
  #stopifnot(any(GWAS.summary$group %in% colnames(Y)))
  
  X <- X[common.samples,] %>% as.data.frame()
  Y <- Y[common.samples,] %>% as.data.frame()
  
  df.out <- vector(mode="list", ncol(Y))
  names(df.out) <- colnames(Y)
  for (n in 1:ncol(Y)){
    print(paste0(n, "/", ncol(Y)))
    
    #gather genotype of significant SNPs
    snp_n <- GWAS.summary$SNP[GWAS.summary$group == colnames(Y)[n]]
    if (length(snp_n)>0){
      idx <- match(snp_n, snp$ID)
      geno_n <- geno[idx,] %>% t() %>% as.data.frame()
      colnames(geno_n) <- paste0("SNP:", snp_n)
      colnames(geno_n) <- gsub(":", "\\.", colnames(geno_n) ) #formula transforms them into . anyway
      
      subject <- gsub("\\:.*", "", rownames(X))
      X.geno <- geno_n[subject,] %>% as.data.frame()
      X.geno[X.geno=="0/0"] <- 0
      X.geno[X.geno=="0/1"] <- 1
      X.geno[X.geno=="1/1"] <- 2
      
      for (k in 1:ncol(X.geno)){
        X.geno[,k] <- as.numeric(X.geno[,k])
      }
      if (length(snp_n)==1){
        colnames(X.geno) <- paste0("SNP.", snp_n)
      }
      df <- data.frame(y = Y[,n] %>% as.numeric(), cbind(X, X.geno)) 
    } else {
      df <- data.frame(y = Y[,n] %>% as.numeric(), X) 
    }
    
    #lm
    lmFit <- lm(formula(df), df)
    lm.df <- summary(lmFit)$coefficient %>% as.data.frame()
    lm.df$var <- rownames(lm.df)
    
    #ANOVA
    af <- anova(lmFit) %>% as.data.frame()
    af$var <- rownames(af)
    afss <- af$"Sum Sq"
    af$PctExp <- afss/sum(afss)*100
    
    #merge
    lm.df <- merge(lm.df, af, by.x="var", by.y="var") #%>% bind_rows(data.frame(var="subject_id", PctExp=af["subject_id", "PctExp"]))
    
    #variable category
    lm.df <- merge(lm.df, var.group, by.x="var", by.y="var")
    
    #store
    df.out[[n]] <- lm.df
  }
  return(df.out)
}

color_palette <- c("darkseagreen",	"#A30059",	"#9e5a28ff",		
                   "tomato",
                   "#4a6fe3",	"lightblue",	"steelblue" ,
                   "bisque",	"#1CE6FF",
                   "goldenrod", "seagreen")
names(color_palette) <- c("Genetics", "Sex", "Age", 
                          "Blood pressure", 
                          "Body composition", "Lipid profile","Glucose homeostasis",
                          "Organ biomarker", "CRP",
                          "Lifestyle", "subject_id")

pdf("varExpl.legend.pdf", width=4, height=2)
ggplot(color_df <- data.frame(Group = names(color_palette),Color = color_palette), aes(x = Group, y = 1, fill = Color)) +
  geom_tile() +
  scale_fill_identity() + # Use exact colors from the palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
  ) 
dev.off()

plot_varExpl <- function(df.out, Nmax=100){
  
  #plot max Nmax
  if (length(df.out)>Nmax){
    tot.pctexp <- rep(0, length(df.out))
    for (n in 1:length(df.out)){
      tot.pctexp[n] <- sum(df.out[[n]]$PctExp)
    }
    df.out <- df.out[order(tot.pctexp, decreasing=T)[1:Nmax]]
  }
  # df <- df.plot %>% group_by(by=x) %>% summarise(sumPct=sum(PctExp)) %>% as.data.frame()
  # x.keep <- df$by[order(df$sumPct, decreasing=T)[1:Nmax]]
  # df.plot <- df.plot[df.plot$x %in% x.keep,]
  
  #gather
  df.plot <- data.frame()
  for (n in 1:length(df.out)){
    df <- df.out[[n]]
    group.uniq <- unique(df$group)
    for (k in 1:length(group.uniq)){
      idx <- which(df$group==group.uniq[k])
      df.plot <- rbind(df.plot,
                       data.frame(x=names(df.out)[n], group=group.uniq[k], PctExp=sum(df$PctExp[idx])))
    }
  }
  
  #order x based on total variance explained
  df <- df.plot %>% na.omit() %>% group_by(by=x) %>% summarise(sumPct=sum(PctExp)) %>% as.data.frame() 
  df.plot$x <- factor(df.plot$x, df$by[order(df$sumPct, decreasing=T)])
  
  #order groups based on PctExp
  df <- df.plot %>% group_by(by=group) %>% summarise(sumPct=sum(PctExp)) %>% as.data.frame()
  df.plot$group <- factor(df.plot$group, levels=df$by[order(df$sumPct, decreasing=F)])
  
  #plot
  p <- ggplot(df.plot, aes(x=x, y=PctExp, fill=group)) +
    geom_bar(stat="identity", position="stack", color="white") +
    scale_fill_manual(values=color_palette) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0)) +
    theme(text=element_text(size=8),
          legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Variance explained") + xlab("")
  
  
  #Plot total variance explained
  df$meanPct <- df$sumPct/length(unique(df.plot$x))
  
  
  return(list(p=p, df=df.plot, df2=df))  
}


# Variance explained of individual immune frequencies ####
common.samples.lifestyle <- intersect(rownames(clinical), rownames(lifestyle))
X <- cbind(clinical[common.samples.lifestyle,], lifestyle[common.samples.lifestyle,]) %>% scale() %>% as.data.frame()
Y <- cytof[common.samples.lifestyle, ] %>% scale() %>% as.data.frame()
lm.cytof.all <- lm_fun(X=X, 
                       Y=Y, 
                       GWAS.summary=GWAS.cytof.LD)
#saveRDS(lm.cytof.all, "lm.cytof.all.RDS")

out.varExpl <- plot_varExpl(lm.cytof.all)

df.plot <- data.frame(out.varExpl$df)
df.plot <- merge(df.plot, macro.anno.my, by.x="x", by.y="row.names")


pop.palette <- c("NK_cells"="turquoise",  "monocytes"="sienna2", "Dendritic_cells"="lemonchiffon1", "basophils"="goldenrod2", "Other"="gray78", 
                 "Naive_B_cells"="dodgerblue2", "B_cells"="dodgerblue4", "CD4pos_naive" ="#B2DF8A",   
                 "CD8pos_naive"="#FB9A99", "CD4pos"="#33A02C", "CD8pos"="#E31A1C")
color_palette <- c(color_palette, pop.palette)


p <- ggplot(df.plot, aes(x=x, y=PctExp, fill=group)) +
  geom_bar(stat="identity", position="stack", color="white") +
  geom_tile(aes(x=x, y=-3, height=3, fill=family)) +
  scale_fill_manual(values=color_palette) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0)) +
  theme(text=element_text(size=8),
        legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Variance explained") + xlab("") + ggtitle("Immune frequencies")
p
pdf("varExpl-pop.pdf", height=7, width=7)
print(p)
dev.off()

#Extract all associations ####
df.lm.all.cytof <- data.frame()
for (n in 1:length(lm.cytof.all)){
  print(paste0(n, "/", length(lm.cytof.all)))
  df.lm.all.cytof <- rbind(df.lm.all.cytof,
                           lm.cytof.all[[n]] %>% mutate(adj.pval=p.adjust(`Pr(>|t|)`, method="BH")) %>% mutate(gene=names(lm.cytof.all)[n]) )
}
df.lm.all.cytof %>% filter(adj.pval<0.05) %>% pull(var) %>% table() %>% sort(decreasing=T) %>% head(n=10)
df.lm.all.cytof %>% filter(adj.pval<0.05) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% head(n=10)

# Variance explained for major immune frequencies ####
GWAS.cytof.group.LD <- GWAS.cytof.LD
for (g in names(pop.group)){
  GWAS.cytof.group.LD$group[GWAS.cytof.group.LD$group %in% pop.group[[g]]] <- g
}

common.samples.lifestyle <- intersect(rownames(cytof), rownames(clinical)) %>% intersect(rownames(lifestyle))
X <- cbind(clinical[common.samples.lifestyle,], lifestyle[common.samples.lifestyle,]) %>% scale() %>% as.data.frame()
Y <- cytof.group[common.samples.lifestyle,] %>% scale() %>% as.data.frame()
lm.cytof.out <- lm_fun(X=X, 
                       Y=Y, 
                       GWAS.summary=GWAS.cytof.group.LD)

out.varExpl <- plot_varExpl(lm.cytof.out)
out.varExpl$p

pdf("varExpl-pop.group.pdf", height=3, width=4)
print(out.varExpl$p)
dev.off()

#saveRDS(lm.cytof.out, "lm.cytof.RDS")
lm.cytof.out <- readRDS("lm.cytof.RDS")

# Variance explained for gene modules ####


gene.module <- cell.gene %>% group_by(pop) %>% summarise(g=list(gene)) 
gene.module <- setNames(gene.module$g, gene.module$pop)
MIN_SIZE_MODULE <- 50
for (cl in names(gene.module)){
  for (cl2 in names(gene.module)){
    if (cl != cl2){
      g1 <- setdiff(gene.module[[cl]], gene.module[[cl2]])
      g2 <- setdiff(gene.module[[cl2]], gene.module[[cl]])
      gene.module[[cl]] <- g1
      gene.module[[cl2]] <- g2
    }
  }
}
gene.module <- gene.module[lapply(gene.module, length) %>% unlist() >=MIN_SIZE_MODULE]

common.samples.lifestyle <- intersect(rownames(clinical), rownames(lifestyle))
df.plot <- data.frame()
#pdf(paste0("VarExpl-all modules", ".pdf"), height=4, width=7)
for (m in names(gene.module)){
  print(m)
  X <- cbind(clinical[common.samples.lifestyle,], lifestyle[common.samples.lifestyle,]) %>% scale() %>% as.data.frame()
  
  Y <- rna.log[common.samples.lifestyle, gene.module[[m]]] %>% scale() %>% as.data.frame()
  lm.out <- lm_fun(X=X, 
                   Y=Y, 
                   GWAS.summary=GWAS.rna.LD)
  
  out.varExpl <- plot_varExpl(lm.out)
  p <- out.varExpl$p + ggtitle(m)
  
  df.plot <- rbind(df.plot, data.frame(out.varExpl$df2, module=m))
  
  print(p)
}
#dev.off()

df <- df.plot %>% group_by(by) %>% summarise(totPct=sum(sumPct)) %>% as.data.frame()
df.plot$by <- factor(df.plot$by, levels=df$by[order(df$totPct, decreasing=F)])
p <- ggplot(df.plot, aes(x=module, y=meanPct, fill=by)) +
  geom_bar(stat="identity", position="stack", width=0.8, color="gray28") +
  scale_fill_manual(values=color_palette) + 
  ylim(c(0,100)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0)) +
  theme(text=element_text(size=10, family="Arial"),
        legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Variance explained") + xlab("")
p
pdf("varExpl-tot.pdf", height=5, width=4)
print(p)
dev.off()


# Genes explained variance ####
common.samples.lifestyle <- intersect(rownames(clinical), rownames(lifestyle))
X <- cbind(clinical[common.samples.lifestyle,], lifestyle[common.samples.lifestyle,]) %>% scale() %>% as.data.frame()
Y <- rna.log[common.samples.lifestyle, ] %>% scale() %>% as.data.frame()
lm.out.rna <- lm_fun(X=X, 
                     Y=Y, 
                     GWAS.summary=GWAS.rna.LD)

saveRDS(lm.out.rna, "lm.rna.RDS")
lm.out.rna <- readRDS("lm.rna.RDS")


# Number of significant associations
df.lm.all <- data.frame()
for (n in 1:length(lm.out.rna)){
  print(paste0(n, "/", length(lm.out.rna)))
  df.lm.all <- rbind(df.lm.all,
                     lm.out.rna[[n]] %>% mutate(adj.pval=p.adjust(`Pr(>|t|)`, method="BH")) %>% filter(adj.pval<0.05) %>% mutate(gene=names(lm.out.rna)[n]) )
}

# Enrichment of associations in each immune module ####
source("C:/Users/albze08/Desktop/postDoc/functions/my_fisher_test.R")
group.uniq <- unique(df.lm.all$group)
df.fisher <- data.frame()
for (pop in names(gene.module)){
  for (g in group.uniq){
    
    x <- df.lm.all %>% filter(adj.pval<0.05) %>% filter(group==g) %>% pull(gene)
    df <- my_fisher_test(gene.module[[pop]], colnames(rna.log), x)
    
    df.fisher <- rbind(df.fisher, data.frame(df, ntot=length(gene.module[[pop]]), group=g, module=pop))
  }
}
df.fisher$adj.pval <- p.adjust(df.fisher$p_value, method="BH")
group.ord <- df.lm.all %>% filter(adj.pval<1) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% names()
df.fisher$group <- factor(df.fisher$group, group.uniq %>% rev())

#cairo_pdf("ass-module-enrich.pdf", height=4, width=5)
ggplot(df.fisher %>% filter(adj.pval<1), aes(x=group, y=module)) +
  geom_tile(aes(fill=-log10(adj.pval))) +
  #scale_fill_gradientn(colors=c("darkred", "tomato", "gray", "gray88"), values=c(0,0.05, 0.1,1), breaks=0.05, labels="0.05") +
  theme_classic() +
  xlab("") + ylab("")
#dev.off()

### OLD ###################
#Variance explained of each immune module ####

# Load GWAS results
# GWAS.cytof.LD <- readRDS("GWAS.cytof.LD.RDS")
# GWAS.rna.LD <- readRDS("GWAS.rna.LD.RDS")
# GWAS.rna.LD$group <- ensg.symbol$symbol[match(GWAS.rna.LD$group, ensg.symbol$ensg)]
# GWAS.rna.LD <- na.omit(GWAS.rna.LD)
# GWAS.cytokine.LD <- readRDS("GWAS.cytokine.LD.RDS")
# GWAS.autoanti.LD <- readRDS("GWAS.autoanti.LD.RDS")

GWAS.cytof.main.LD <- GWAS.cytof.LD
GWAS.cytof.main.LD$group <- GWAS.cytof.main.LD$family

# Genes associated with each immune module
immune.cl <- c(`Central_Memory_CD4`= "Module 1", `TEMRA_CD4`= "Module 1",
               `B_Cells`= "Module 2", `Naive_B_Cells`= "Module 2", 
               `Myeloid_DC`="Module 3", `Monocytes_classical`= "Module 3",
               `Monocytes_nonclassical`="Module 4", `Plasmacytoid_DC`="Module 4", 
               `Naive_CD8`="Module 5", `Effector_Memory_CD4`= "Module 5", `Naive_CD4`= "Module 5",
               `Central_Memory_CD8`= "Module 6",`Basophils`= "Module 6", 
               `TEMRA_CD8`= "Module 7", `NK_Cells`= "Module 7", `Monocytes_intermediate`= "Module 7",
               `ILCs`= "Module x", `Effector_Memory_CD8`= "Module x")

cl.uniq <- unique(immune.cl)


cell.gene.sub <- readRDS("ILR-cytof.group-310125.RDS")$cell.gene.sub
gene.module <- vector(mode="list", length(cl.uniq))
names(gene.module) <- cl.uniq
MIN_SIZE_MODULE <- 50
for (cl in cl.uniq){
  pop <- names(immune.cl)[immune.cl==cl]
  gene.module[[cl]] <- cell.gene.sub$gene[cell.gene.sub$pop %in% pop]  
}
for (cl in cl.uniq){
  for (cl2 in cl.uniq){
    if (cl != cl2){
      g1 <- setdiff(gene.module[[cl]], gene.module[[cl2]])
      g2 <- setdiff(gene.module[[cl2]], gene.module[[cl]])
      gene.module[[cl]] <- g1
      gene.module[[cl2]] <- g2
    }
  }
}
gene.module <- gene.module[lapply(gene.module, length) %>% unlist() >=MIN_SIZE_MODULE]


# Gene - population ####
pop.uniq <- unique(cell.gene.sub$pop)
gene.pop <- vector(mode="list", length(pop.uniq))
names(gene.pop) <- pop.uniq
MIN_SIZE_MODULE <- 50
for (pop in pop.uniq){
  gene.pop[[pop]] <- cell.gene.sub$gene[cell.gene.sub$pop == pop]  
}
for (pop in pop.uniq){
  for (pop2 in pop.uniq){
    if (pop != pop2){
      g1 <- setdiff(gene.pop[[pop]], gene.pop[[pop2]])
      g2 <- setdiff(gene.pop[[pop2]], gene.pop[[pop]])
      gene.pop[[pop]] <- g1
      gene.pop[[pop2]] <- g2
    }
  }
}
gene.pop <- gene.pop[lapply(gene.pop, length) %>% unlist() >= MIN_SIZE_MODULE]


# Enrichment of eQTLs among immune populations in the network ####
g <- gene.module[[2]]
my_fisher_test(g, colnames(rna.log), unique(GWAS.rna.LD$group))

lapply(gene.pop, function(g){my_fisher_test(g, colnames(rna.log), unique(GWAS.rna.LD$group))})

df.fisher <- data.frame()
for (pop in names(gene.pop)){
  g <- gene.pop[[pop]]
  f <- my_fisher_test(g, colnames(rna.log), unique(GWAS.rna.LD$group))
  df.fisher <- rbind(df.fisher, 
                     data.frame(f, pop=pop, n.pop=length(g)))
  
}
df.fisher$adj.pval <- p.adjust(df.fisher$p_value, method="BH")
df.fisher <- df.fisher %>% filter(adj.pval<0.05)





# Big linear model
var.group <- data.frame(var=c("subject_id",
                              "Age_at_Visit", 
                              "Gender", 
                              "Height", "Weight", "Waist", "Hip", "BMI", "Bioimp_fat", "Bioimp_muscle", "Bioimp_bone", 
                              "SBP", "DBP", 
                              "Cap_Gluc", "Gluc", "HbA1c", 
                              "Chol", "LDL", "HDL", "TG", "ApoB", "ApoA1", "ApoB.apoA1",   
                              "CRP", 
                              "ALAT", "GGT", 
                              "ProBNP", "TNT", 
                              "Urate", 
                              "CystC", "Crea",
                              "Smoking", 
                              "MaritalStatus_change", "Housing_change", "Employment_change",
                              "Stress",
                              "PhysicalActivity", "SedentaryTime_hours",
                              "PerceivedHealth",
                              "Common_cold_influenzae",
                              "NSAID_painmed", "Bp_med", "Lipid_med", "Antibiotics_med"),
                        group=c("subject_id",
                                "Age", 
                                "Sex", 
                                "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", 
                                "Blood pressure", "Blood pressure", 
                                "Glucose homeostasis", "Glucose homeostasis", "Glucose homeostasis", 
                                "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile",   
                                "CRP", 
                                "Organ biomarker", "Organ biomarker", 
                                "Organ biomarker", "Organ biomarker", 
                                "Organ biomarker", 
                                "Organ biomarker", "Organ biomarker",
                                "Lifestyle", 
                                "Lifestyle", "Lifestyle", "Lifestyle",
                                "Lifestyle",
                                "Lifestyle", "Lifestyle",
                                "Lifestyle",
                                "Lifestyle",
                                "Lifestyle", "Lifestyle", "Lifestyle", "Lifestyle"))
var.group.SNP <- data.frame(var=unique(c(GWAS.cytof.main.LD$SNP, GWAS.cytof.LD$SNP, GWAS.rna.LD$SNP)),
                            group="Genetics")
var.group.SNP$var <- paste0("SNP:", var.group.SNP$var)
var.group.SNP$var <- gsub(":", "\\.", var.group.SNP$var)

var.group <- rbind(var.group, var.group.SNP)

lm_fun <- function(X,Y, GWAS.summary){
  common.samples <- intersect(rownames(X), rownames(Y))
  stopifnot(length(common.samples)>0)
  #stopifnot(any(GWAS.summary$group %in% colnames(Y)))
  
  X <- X[common.samples,] %>% as.data.frame()
  Y <- Y[common.samples,] %>% as.data.frame()
  
  df.out <- vector(mode="list", ncol(Y))
  names(df.out) <- colnames(Y)
  for (n in 1:ncol(Y)){
    print(paste0(n, "/", ncol(Y)))
    
    #gather genotype of significant SNPs
    snp_n <- GWAS.summary$SNP[GWAS.summary$group == colnames(Y)[n]]
    if (length(snp_n)>0){
      idx <- match(snp_n, snp$ID)
      geno_n <- geno[idx,] %>% t() %>% as.data.frame()
      colnames(geno_n) <- paste0("SNP:", snp_n)
      colnames(geno_n) <- gsub(":", "\\.", colnames(geno_n) ) #formula transforms them into . anyway
      
      subject <- gsub("\\:.*", "", rownames(X))
      X.geno <- geno_n[subject,] %>% as.data.frame()
      X.geno[X.geno=="0/0"] <- 0
      X.geno[X.geno=="0/1"] <- 1
      X.geno[X.geno=="1/1"] <- 2
      
      for (k in 1:ncol(X.geno)){
        X.geno[,k] <- as.numeric(X.geno[,k])
      }
      if (length(snp_n)==1){
        colnames(X.geno) <- paste0("SNP.", snp_n)
      }
      df <- data.frame(y = Y[,n] %>% as.numeric(), cbind(X, X.geno)) 
    } else {
      df <- data.frame(y = Y[,n] %>% as.numeric(), X) 
    }
    
    #lm
    lmFit <- lm(formula(df), df)
    lm.df <- summary(lmFit)$coefficient %>% as.data.frame()
    lm.df$var <- rownames(lm.df)
    
    #ANOVA
    af <- anova(lmFit) %>% as.data.frame()
    af$var <- rownames(af)
    afss <- af$"Sum Sq"
    af$PctExp <- afss/sum(afss)*100
    
    #merge
    lm.df <- merge(lm.df, af, by.x="var", by.y="var") #%>% bind_rows(data.frame(var="subject_id", PctExp=af["subject_id", "PctExp"]))
    
    #variable category
    lm.df <- merge(lm.df, var.group, by.x="var", by.y="var")
    
    #store
    df.out[[n]] <- lm.df
  }
  return(df.out)
}

color_palette <- c("darkseagreen",	"#A30059",	"#9e5a28ff",		
                   "tomato",
                   "#4a6fe3",	"lightblue",	"steelblue" ,
                   "bisque",	"#1CE6FF",
                   "goldenrod", "seagreen")
names(color_palette) <- c("Genetics", "Sex", "Age", 
                          "Blood pressure", 
                          "Body composition", "Lipid profile","Glucose homeostasis",
                          "Organ biomarker", "CRP",
                          "Lifestyle", "subject_id")

pdf("varExpl.legend.pdf", width=4, height=2)
ggplot(color_df <- data.frame(Group = names(color_palette),Color = color_palette), aes(x = Group, y = 1, fill = Color)) +
  geom_tile() +
  scale_fill_identity() + # Use exact colors from the palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
  ) 
dev.off()

plot_varExpl <- function(df.out, Nmax=100){
  
  #plot max Nmax
  if (length(df.out)>Nmax){
    tot.pctexp <- rep(0, length(df.out))
    for (n in 1:length(df.out)){
      tot.pctexp[n] <- sum(df.out[[n]]$PctExp)
    }
    df.out <- df.out[order(tot.pctexp, decreasing=T)[1:Nmax]]
  }
  # df <- df.plot %>% group_by(by=x) %>% summarise(sumPct=sum(PctExp)) %>% as.data.frame()
  # x.keep <- df$by[order(df$sumPct, decreasing=T)[1:Nmax]]
  # df.plot <- df.plot[df.plot$x %in% x.keep,]
  
  #gather
  df.plot <- data.frame()
  for (n in 1:length(df.out)){
    df <- df.out[[n]]
    group.uniq <- unique(df$group)
    for (k in 1:length(group.uniq)){
      idx <- which(df$group==group.uniq[k])
      df.plot <- rbind(df.plot,
                       data.frame(x=names(df.out)[n], group=group.uniq[k], PctExp=sum(df$PctExp[idx])))
    }
  }
  
  #order x based on total variance explained
  df <- df.plot %>% na.omit() %>% group_by(by=x) %>% summarise(sumPct=sum(PctExp)) %>% as.data.frame() 
  df.plot$x <- factor(df.plot$x, df$by[order(df$sumPct, decreasing=T)])
  
  #order groups based on PctExp
  df <- df.plot %>% group_by(by=group) %>% summarise(sumPct=sum(PctExp)) %>% as.data.frame()
  df.plot$group <- factor(df.plot$group, levels=df$by[order(df$sumPct, decreasing=F)])
  
  #plot
  p <- ggplot(df.plot, aes(x=x, y=PctExp, fill=group)) +
    geom_bar(stat="identity", position="stack", color="white") +
    scale_fill_manual(values=color_palette) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0)) +
    theme(text=element_text(size=8),
          legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Variance explained") + xlab("")
  
  
  #Plot total variance explained
  df$meanPct <- df$sumPct/length(unique(df.plot$x))
  
  
  return(list(p=p, df=df.plot, df2=df))  
}


# Variance explained of individual immune frequencies
common.samples.lifestyle <- intersect(rownames(clinical), rownames(lifestyle))
X <- cbind(clinical[common.samples.lifestyle,], lifestyle[common.samples.lifestyle,]) %>% scale() %>% as.data.frame()
Y <- cytof[common.samples.lifestyle, ] %>% scale() %>% as.data.frame()
lm.cytof.all <- lm_fun(X=X, 
                 Y=Y, 
                 GWAS.summary=GWAS.cytof.LD)
#saveRDS(lm.cytof.all, "lm.cytof.all.RDS")

out.varExpl <- plot_varExpl(lm.cytof.all)

df.plot <- data.frame(out.varExpl$df)
df.plot <- merge(df.plot, macro.anno.my, by.x="x", by.y="row.names")


pop.palette <- c("NK_cells"="turquoise",  "monocytes"="sienna2", "Dendritic_cells"="lemonchiffon1", "basophils"="goldenrod2", "Other"="gray78", 
                 "Naive_B_cells"="dodgerblue2", "B_cells"="dodgerblue4", "CD4pos_naive" ="#B2DF8A",   
                 "CD8pos_naive"="#FB9A99", "CD4pos"="#33A02C", "CD8pos"="#E31A1C")
color_palette <- c(color_palette, pop.palette)


p <- ggplot(df.plot, aes(x=x, y=PctExp, fill=group)) +
  geom_bar(stat="identity", position="stack", color="white") +
  geom_tile(aes(x=x, y=-3, height=3, fill=family)) +
  scale_fill_manual(values=color_palette) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0)) +
  theme(text=element_text(size=8),
        legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Variance explained") + xlab("") + ggtitle("Immune frequencies")
p
pdf("varExpl-pop.pdf", height=7, width=7)
print(p)
dev.off()

#Extract all associations ####
df.lm.all.cytof <- data.frame()
for (n in 1:length(lm.cytof.all)){
  print(paste0(n, "/", length(lm.cytof.all)))
  df.lm.all.cytof <- rbind(df.lm.all.cytof,
                           lm.cytof.all[[n]] %>% mutate(adj.pval=p.adjust(`Pr(>|t|)`, method="BH")) %>% mutate(gene=names(lm.cytof.all)[n]) )
}
df.lm.all.cytof %>% filter(adj.pval<0.05) %>% pull(var) %>% table() %>% sort(decreasing=T) %>% head(n=10)
df.lm.all.cytof %>% filter(adj.pval<0.05) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% head(n=10)

# Variance explained for major immune frequencies
GWAS.cytof.group.LD <- GWAS.cytof.LD
for (g in names(pop.group)){
  GWAS.cytof.group.LD$group[GWAS.cytof.group.LD$group %in% pop.group[[g]]] <- g
}

common.samples.lifestyle <- intersect(rownames(cytof), rownames(clinical)) %>% intersect(rownames(lifestyle))
X <- cbind(clinical[common.samples.lifestyle,], lifestyle[common.samples.lifestyle,]) %>% scale() %>% as.data.frame()
Y <- cytof.group[common.samples.lifestyle,] %>% scale() %>% as.data.frame()
lm.cytof.out <- lm_fun(X=X, 
                 Y=Y, 
                 GWAS.summary=GWAS.cytof.group.LD)

out.varExpl <- plot_varExpl(lm.cytof.out)
out.varExpl$p

pdf("varExpl-pop.group.pdf", height=3, width=4)
print(out.varExpl$p)
dev.off()

#saveRDS(lm.cytof.out, "lm.cytof.RDS")
lm.cytof.out <- readRDS("lm.cytof.RDS")

# Variance explained for gene modules
common.samples.lifestyle <- intersect(rownames(clinical), rownames(lifestyle))
df.plot <- data.frame()
pdf(paste0("VarExpl-all modules", ".pdf"), height=4, width=7)
for (m in names(gene.module)){
  print(m)
  X <- cbind(clinical[common.samples.lifestyle,], lifestyle[common.samples.lifestyle,]) %>% scale() %>% as.data.frame()

  Y <- rna.log[common.samples.lifestyle, gene.module[[m]]] %>% scale() %>% as.data.frame()
  lm.out <- lm_fun(X=X, 
                   Y=Y, 
                   GWAS.summary=GWAS.rna.LD)
  
  out.varExpl <- plot_varExpl(lm.out)
  p <- out.varExpl$p + ggtitle(m)
  
  df.plot <- rbind(df.plot, data.frame(out.varExpl$df2, module=m))
  
  print(p)
}
dev.off()

df <- df.plot %>% group_by(by) %>% summarise(totPct=sum(sumPct)) %>% as.data.frame()
df.plot$by <- factor(df.plot$by, levels=df$by[order(df$totPct, decreasing=F)])
p <- ggplot(df.plot, aes(x=module, y=meanPct, fill=by)) +
  geom_bar(stat="identity", position="stack", width=0.8, color="gray28") +
  scale_fill_manual(values=color_palette) + 
  ylim(c(0,100)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0)) +
  theme(text=element_text(size=10, family="Arial"),
        legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Variance explained") + xlab("")
p
pdf("varExpl-tot.pdf", height=5, width=4)
print(p)
dev.off()


# Genes explained variance ####
common.samples.lifestyle <- intersect(rownames(clinical), rownames(lifestyle))
X <- cbind(clinical[common.samples.lifestyle,], lifestyle[common.samples.lifestyle,]) %>% scale() %>% as.data.frame()
Y <- rna.log[common.samples.lifestyle, ] %>% scale() %>% as.data.frame()
lm.out.rna <- lm_fun(X=X, 
                 Y=Y, 
                 GWAS.summary=GWAS.rna.LD)

#saveRDS(lm.out.rna, "lm.rna.RDS")
lm.rna.out <- readRDS("lm.rna.RDS")


# Number of significant associations
df.lm.all <- data.frame()
for (n in 1:length(lm.out.rna)){
  print(paste0(n, "/", length(lm.out.rna)))
  df.lm.all <- rbind(df.lm.all,
                   lm.out.rna[[n]] %>% mutate(adj.pval=p.adjust(`Pr(>|t|)`, method="BH")) %>% mutate(gene=names(lm.out.rna)[n]) )
}
#saveRDS(df.lm.all, "lm.rna.all.RDS")
df.lm.all <- readRDS("lm.rna.all.RDS")

df.lm.all %>% filter(adj.pval<0.05) %>% pull(var) %>% table() %>% sort(decreasing=T) %>% head(n=10)
df.lm.all %>% filter(adj.pval<0.05) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% head(n=10)

df.lm.all %>% group_by(gene) %>% summarise(nass=sum(adj.pval<0.05)) %>% as.data.frame() %>% pull(nass) %>% table()

group.uniq <-  unique(df.lm.all$group)
df.fisher <- data.frame()
for (g in group.uniq){
  g.g <- df.lm.all %>% filter(adj.pval<0.05) %>% filter(group==g) %>% pull(gene)
  f.vec <- lapply(gene.pop, function(x){my_fisher_test(x, colnames(rna.log), g.g)}) 
  for (n in 1:length(f.vec)){
      df.fisher <- rbind(df.fisher, 
                         data.frame(f.vec[[n]], pop=names(f.vec)[n], group=g, 
                                    n.pop=length(gene.pop[[names(f.vec)[n]]])))
    
  }
}
df.fisher$adj.pval <- p.adjust(df.fisher$p_value, method="BH")
#df.fisher <- df.fisher %>% filter(adj.pval<0.05)

df.fisher$prc <- 100*df.fisher$n_genes/df.fisher$n.pop
p <- ggplot(df.fisher, aes(x=pop, y=prc)) +
  facet_wrap(~group,ncol=1)+
  geom_bar(stat="identity", position="dodge") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab("Percentage of genes") + xlab("")
p

pdf("prc_gene_var.pdf", height=10, width=4)
print(p)
dev.off()

# Prepare Suppl Tables of linear model ####
df.lm.all.cytof <- data.frame()
for (n in 1:length(lm.cytof.out)){
  print(paste0(n, "/", length(lm.cytof.out)))
  df.lm.all.cytof <- rbind(df.lm.all.cytof,
                           lm.cytof.out[[n]] %>% mutate(adj.pval=p.adjust(`Pr(>|t|)`, method="BH")) %>% mutate(gene=names(lm.cytof.out)[n]) )
}
#writexl::write_xlsx(list(df.lm.all.cytof %>% filter(adj.pval<0.05),df.lm.all %>% filter(adj.pval<0.05)), "Supplementary Table 6_lm.xlsx")
df.lm.all.cytof %>% filter(adj.pval<0.05) %>% pull(var) %>% table() %>% sort(decreasing=T) %>% head(n=10)
df.lm.all.cytof %>% filter(adj.pval<0.05) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% head(n=10)

# N. significant associations per variable ####

df.plot <- df.lm.all %>% filter(adj.pval<0.05) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% as.data.frame()
colnames(df.plot) <- c("var", "n")
df.plot$var <- as.character(df.plot$var)
df.plot$var <- factor(df.plot$var, levels=df.plot$var %>% rev())

p1 <- ggplot(df.plot, aes(y=var, x=n, fill=var)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_manual(values=color_palette) +
  theme_classic() + xlab("Number of genes") +
  theme(legend.position="bottom")


df.plot <- df.lm.all.cytof %>% filter(adj.pval<0.05) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% as.data.frame()
colnames(df.plot) <- c("var", "n")
df.plot$var <- as.character(df.plot$var)
df.plot$var <- factor(df.plot$var, levels=df.plot$var %>% rev())

p2 <- ggplot(df.plot, aes(y=var, x=n, fill=var)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_manual(values=color_palette) +
  theme_classic() + xlab("Number of immune populations")  +
  theme(legend.position="bottom")

pdf("lm_nass.pdf", width=8, height=4)
ggarrange(p1,p2, nrow=1)
dev.off()


# PBMC changes from variables ####
group.uniq <- unique(df.fisher$group)
for (group.n in group.uniq) {
  gene.group.n <- df.lm.all %>% filter(adj.pval<0.05) %>% filter(group==group.n) %>% pull(gene)

  pop.n <- df.fisher$pop[df.fisher$adj.pval<0.05 & df.fisher$group==group.n]
  
  #Overlap of genes associated with each pop
  if (length(pop.n)>=2){
    pop.comb <- combn(pop.n,2)
    gene.overlap <- list()
    k <- 1
    for (n in 1:ncol(pop.comb)){
      pop1 <- pop.comb[1,n]
      pop2 <- pop.comb[2,n]
      pop1.g <- intersect(gene.group.n, cell.gene.sub$gene[cell.gene.sub$pop==pop1])
      pop2.g <- intersect(gene.group.n, cell.gene.sub$gene[cell.gene.sub$pop==pop2])
      gene.overlap[[k]] <- intersect(pop1.g, pop2.g)
      names(gene.overlap)[k] <- paste0(pop1, ";", pop2)
      k <- k + 1
    }
    gene.overlap <- gene.overlap[( lapply(gene.overlap, length) %>% unlist() )>0]
  } else {gene.overlap <- list()}
  
  #Remaining unique genes
  gene.uniq <- list()
  for (pop in pop.n){
    gene.uniq[[pop]] <- intersect(gene.group.n, cell.gene.sub$gene[cell.gene.sub$pop==pop])
    for (k in grep(pop, names(gene.overlap))){
      gene.uniq[[pop]] <- setdiff(gene.uniq[[pop]], gene.overlap[[k]])
    }
  }
  
  #IGRAPH
  node_df <- data.frame(node=c(names(gene.uniq),
                               paste0(names(gene.uniq), "_gene")),
                        label=c(names(gene.uniq),
                                lapply(gene.uniq, length) %>% unlist()))
  
  if (length(gene.overlap)>0){
    node_df <- rbind(node_df, data.frame(node=paste(names(gene.overlap), "edge", sep="_"), label=lapply(gene.overlap, length) %>% unlist()))
  }
  
  
  #size
  node_df$size <- 6
  for (n in grep("gene|edge", node_df$node)){
    n.genes <- as.numeric(node_df$label[n])
    if (n.genes>499){
      node_df$size[n] <- 12
    } else if (n.genes>99){
      node_df$size[n] <- 10
    } else if (n.genes>9){
      node_df$size[n] <- 8
    }
  }
  node_df$size[node_df$node %in% names(gene.pop)] <- 10
  
  #color
  pop.palette <- c("Monocytes_classical"="sienna2", "TEMRA_CD8"="#E31A1C", "NK_Cells"="turquoise", 
                   "Naive_CD8"="#FB9A99", "Naive_CD4"="#B2DF8A",             
                   "Basophils"="goldenrod2", "Monocytes_nonclassical"="sienna2", "Central_Memory_CD8"="#E31A1C", 
                   "Monocytes_intermediate"="sienna2", "B_Cells"="dodgerblue4",               
                   "Myeloid_DC"="lemonchiffon1", "Plasmacytoid_DC"="lemonchiffon1", "Naive_B_Cells"="dodgerblue2", 
                   "Central_Memory_CD4"="#33A02C", "Effector_Memory_CD4"="#33A02C",   
                   "TEMRA_CD4"="#33A02C", "Effector_Memory_CD8"="#E31A1C", "ILCs"="gray78")
  pop.palette <- pop.palette[names(pop.palette) %in% node_df$node]
  node_df$color <- "lightsteelblue"
  node_df$color[match(names(pop.palette), node_df$node)] <- pop.palette
  
  
  #edge df, for interactions
  edge_df <- data.frame()
  if (length(gene.overlap)>0){
    for (n in 1:length(gene.overlap)){
      pop1 <- strsplit(names(gene.overlap)[n], ";")[[1]][1]
      pop2 <- strsplit(names(gene.overlap)[n], ";")[[1]][2]
      edge_node <- paste0(pop1, ";", pop2, "_edge")
      edge_df <- rbind(edge_df,
                       data.frame(source=c(pop1,edge_node),
                                  target=c(edge_node, pop2),
                                  weight=node_df$label[node_df$node==edge_node]))
    }
  }
  
  #edge df, pop-specific genes
  for (pop in names(gene.uniq)){
    edge_df <- rbind(edge_df, 
                     data.frame(source=pop, 
                                target=paste0(pop, "_gene"),
                                weight=node_df$label[node_df$node==paste0(pop, "_gene")]))
  }
  edge_df$weight <- as.numeric(edge_df$weight)
  
  net <- graph_from_data_frame(d=edge_df, vertices = node_df, directed=F)
  layout <- layout_nicely(net)
  
  
  #Adjust size of gene nodes using ratio
  for (n in grep("gene", V(net)$name)){
    node_name <- V(net)$name[n]
    generatio <- 0
    for (pop in pop.n){
      if (length(grep(pop, node_name))>0){
        g.pop <- cell.gene.sub$gene[cell.gene.sub$pop==pop]
        Ntot <- length(g.pop)
        Nnode <- length(intersect(g.pop, gene.group.n))
      }
    }
    V(net)$size[n] <- 6 + 6*Nnode/Ntot
  }
  
  # for (n in grep("gene", V(net)$name)){
  #   
  #   pop <- gsub("_gene", "", V(net)$name[n])
  #   g.node <- gene.uniq[[pop]]
  #   
  #   f <- my_fisher_test(g.node, colnames(rna), gene.group.n)
  #   
  #   node_name <- V(net)$name[n]
  #   generatio <- 0
  #   for (pop in pop.n){
  #     if (length(grep(pop, node_name))>0){
  #       g.pop <- cell.gene.sub$gene[cell.gene.sub$pop==pop] 
  #       Ntot <- length(g.pop)
  #       Nnode <- length(intersect(g.pop, gene.group.n))
  #     }
  #   }
  #   V(net)$size[n] <- 6 + 6*Nnode/Ntot
  # }
  
  #Adjust size of pop nodes using varExpl
  for (pop in pop.n){
    ord <- match(pop, V(net)$name)
    df.pop <- lm.cytof.out[[pop]]
    idx.group <- which(df.pop$group==group.n)
    prc <- df.pop$PctExp[idx.group] %>% sum()

    V(net)$size[ord] <- 10 + 10*prc/100
  }
  
  plot(net, layout=layout, vertex.label.cex=.6, 
       vertex.label.color="black", edge.label.color="black", edge.label.cex=.8,
       main=group.n)

}


# Enrichment of associations in each immune module ####
group.uniq <- unique(df.lm.all$group)
df.fisher <- data.frame()
for (pop in names(gene.module)){
  for (g in group.uniq){
    
    x <- df.lm.all %>% filter(adj.pval<0.05) %>% filter(group==g) %>% pull(gene)
    df <- my_fisher_test(gene.module[[pop]], colnames(rna.log), x)
    
    df.fisher <- rbind(df.fisher, data.frame(df, ntot=length(gene.module[[pop]]), group=g, module=pop))
  }
}
df.fisher$adj.pval <- p.adjust(df.fisher$p_value, method="BH")
group.ord <- df.lm.all %>% filter(adj.pval<1) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% names()
df.fisher$group <- factor(df.fisher$group, group.uniq %>% rev())
cairo_pdf("ass-module-enrich.pdf", height=4, width=5)
ggplot(df.fisher %>% filter(adj.pval<1), aes(x=group, y=module)) +
  geom_tile(aes(fill=adj.pval)) +
  scale_fill_gradientn(colors=c("darkred", "tomato", "gray", "gray88"), values=c(0,0.05, 0.1,1), breaks=0.05, labels="0.05") +
  theme_classic() +
  xlab("") + ylab("")
dev.off()

# Enrichment of associations in each population ####
group.uniq <- unique(df.lm.all$group)
df.fisher <- data.frame()
for (pop in unique(cell.gene.sub$pop)){
  for (g in group.uniq){
    
    x <- df.lm.all %>% filter(adj.pval<0.05) %>% filter(group==g) %>% pull(gene)
    df <- my_fisher_test(cell.gene.sub$gene[cell.gene.sub$pop==pop], colnames(rna.log), x)
    
    df.fisher <- rbind(df.fisher, data.frame(df, ntot=sum(cell.gene.sub$pop==pop), group=g, pop=pop))
  }
}
df.fisher$adj.pval <- p.adjust(df.fisher$p_value, method="BH")
group.ord <- df.lm.all %>% filter(adj.pval<1) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% names()
df.fisher$group <- factor(df.fisher$group, group.uniq %>% rev())
cairo_pdf("ass-pop-enrich.pdf", height=4, width=5)
ggplot(df.fisher %>% filter(adj.pval<1), aes(x=group, y=pop)) +
  geom_tile(aes(fill=adj.pval)) +
  scale_fill_gradientn(colors=c("darkred", "tomato", "gray", "gray88"), values=c(0,0.05, 0.1,1), breaks=0.05, labels="0.05") +
  theme_classic() +
  xlab("") + ylab("")
dev.off()


# Modify PBMC network to reflect effects ####
PBMC.network <- readRDS("PBMC-network.RDS")
net.ori <- PBMC.network$net
layout.ori <- PBMC.network$layout

group.vec <- df.lm.all %>% filter(adj.pval<0.05) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% head(n=5)
for (group.n in names(group.vec) ) {
  #group.n <- "Body composition"
  gene.group.n <- df.lm.all %>% filter(adj.pval<0.05) %>% filter(group==group.n) %>% pull(gene)
  
  net.group <- net.ori
  layout.group <- net.ori
  
  #for pop, use varExpl as size
  ord <- match(names(lm.cytof.out), V(net.ori)$name)
  for (n in 1:length(lm.cytof.out)){
    df.pop <- lm.cytof.out[[n]]
    idx.group <- which(df.pop$group==group.n)
    prc <- df.pop$PctExp[idx.group] %>% sum()
    
    V(net.group)$size[ord[n]] <- max(4,V(net.group)$size[ord[n]]*prc/100)
  }
  
  #for genes, use ratio of total genes
  scaling_function <- function(x){
    if (x>0.05){
      y <- 4
    } else {
      y <- min(8, -4*log10(x))
    }
  }
  pop.uniq <- unique(cell.gene.sub$pop)
  g.back <- df.lm.all$gene %>% unique()
  for (n in grep("gene", V(net.ori)$name)){
    node_name <- V(net.ori)$name[n]
    generatio <- 0
    for (pop in pop.uniq){
      if (length(grep(pop, node_name))>0){
        g.pop <- cell.gene.sub$gene[cell.gene.sub$pop==pop]
        Ntot <- length(g.pop)
        Nnode <- length(intersect(g.pop, gene.group.n))
        if (Nnode/Ntot>generatio){
          scaling_factor <- Nnode/Ntot
        }
        
        pval <- df.fisher$adj.pval[df.fisher$pop==pop & df.fisher$group==group.n]
        scaling_factor <- scaling_function(pval)
        if (pval>0.05){
          V(net.group)$color[n] <- "gray77"
        }
      }
    }
    V(net.group)$size[n] <- scaling_factor
  }
  
  plot(net.group, layout=layout.ori, 
       vertex.label.cex=.6, vertex.label=NA, vertex.label.color="black", 
       edge.label.color="black", edge.label.cex=.8,
       label.family="Arial")
  
  pdf(paste0("PBMC-network-", group.n, ".pdf"), height=7, width=7)
  plot(net.group, layout=layout.ori, 
       vertex.label.cex=.6, vertex.label=NA, vertex.label.color="black", 
       edge.label.color="black", edge.label.cex=.8,
       main=group.n)
  dev.off()
}

pdf("PBMC-net.pdf", height=7, width=7)
plot(net.ori, layout=layout.ori, 
     vertex.label.cex=.6, vertex.label.color="black", 
     edge.label.color="black", edge.label.cex=.8,
     label.family="Arial")
dev.off()

# Map changes to the PBMC network ####
group.uniq <- df.lm.all$group %>% unique()
cell.gene <- readRDS("ILR-cytof.group-310125.RDS")$cell.gene.sub[,1:5]

group.n <- "Genetics"
gene.group.n <- df.lm.all %>% filter(adj.pval<0.05) %>% filter(group==group.n) %>% pull(gene)
net.group <- cell.gene %>% filter(gene %in% gene.group.n) 
net.group <- net.group[!net.group$pop %in% c("Effector_Memory_CD8", "ILCs"),] 

#select pairs with significant overlap
source('C:/Users/albze08/Desktop/postDoc/functions/my_fisher_test.R')
pop.uniq <- unique(net.group$pop)
k <- 1
df <- data.frame()
for (pop1 in pop.uniq){
  for (pop2 in pop.uniq){
    if (pop1 != pop2){
      if (!any(df$source==pop1 & df$target==pop2)){
        if (!any(df$target==pop1 & df$source==pop2)){
          
          pop1.g <- net.group$gene[net.group$pop==pop1]
          pop2.g <- net.group$gene[net.group$pop==pop2]
          f <- my_fisher_test(pop1.g, colnames(rna), pop2.g)
          df <- rbind(df, 
                      data.frame(source=pop1, target=pop2, pval=f$p_value, OR=f$odds_ratio, 
                                 n1=length(pop1.g), n2=length(pop2.g), n=f$n_genes))
          k <- k + 1
        }
      }
    }
  }
}
df$adj.pval <- p.adjust(df$pval, method="BH")
df$weights <- df$n
df_overlap_pop <- df
df <- df[df$adj.pval<0.05,]

#genes in edges
gene.overlap <- vector(mode="list", nrow(df))
for (n in 1:length(gene.overlap)){
  pop1 <- df$source[n]
  pop2 <- df$target[n]
  g  <- intersect(net.group$gene[net.group$pop==pop1],
                  net.group$gene[net.group$pop==pop2])
  gene.overlap[[n]] <- g
}

#genes that are unique per cell pop
gene.pop <- vector(mode="list", length(unique(net.group$pop)))
names(gene.pop) <- unique(net.group$pop)
for (pop in names(gene.pop)){
  idx <- which(df$source == pop | df$target == pop)
  gene.pop[[pop]] <- setdiff(net.group$gene[net.group$pop==pop], gene.overlap[idx] %>% unlist())
}
gene.pop.uniq <- gene.pop
for (pop in names(gene.pop)){
  for (pop2 in names(gene.pop)){
    if (pop!=pop2){
      gene.pop.uniq[[pop]] <- setdiff(gene.pop.uniq[[pop]], gene.pop[[pop2]])
    }
  }
}
gene.pop <- gene.pop.uniq

# node df
MIN_GENES <- 0 #if less, plot as individual genes
node_df <- data.frame(node=names(gene.pop), label=names(gene.pop))
for (n in 1:length(gene.pop)){
  g <- gene.pop[[n]]
  if (length(g)>0){
    if (length(g)<MIN_GENES){
      node_df <- rbind(node_df, 
                       data.frame(node=g, label=""))
    } else {
      node_df <- rbind(node_df, 
                       data.frame(node=paste0(names(gene.pop)[n], "_gene"),
                                  label=length(g)))
    }
  }
}
for (n in 1:nrow(df)){
  pop1 <- df$source[n]
  pop2 <- df$target[n]
  g  <- intersect(net.group$gene[net.group$pop==pop1],
                  net.group$gene[net.group$pop==pop2])
  if (length(g)<MIN_GENES){
    node_df <- rbind(node_df, 
                     data.frame(node=g, label=""))
  } else {
    node_df <- rbind(node_df, 
                     data.frame(node=paste0(pop1, "_", pop2, "_edge_gene"),
                                label=length(g)))
  }
}

node_df <- node_df[!duplicated(node_df),]

#size of nodes
node_df$size <- 1 #for individual genes
pop.uniq <- unique(cell.gene$pop)
for (n in grep("gene", node_df$node)){
  node_name <- node_df$node[n]
  generatio <- 0
  for (pop in pop.uniq){
    if (length(grep(pop, node_name))>0){
      Ntot <- sum(pop == cell.gene.sub$pop)
      Nnode <- as.numeric(node_df$label[n])
      if (Nnode/Ntot>generatio){
        generatio <- Nnode/Ntot
        node_df$size[n] <- 4*10*generatio
      }
    }
  }
}
for (pop in pop.uniq){
  idx <- which(node_df$node==pop)
  df.pop <- lm.cytof.out[[pop]]
  idx.group <- which(df.pop$group==group.n)
  prc <- df.pop$PctExp[idx.group] %>% sum()
  
  node_df$size[idx] <- max(2,1.5*prc/10)
}

#color of nodes
pop.palette <- c("Monocytes_classical"="sienna2", "TEMRA_CD8"="#E31A1C", "NK_Cells"="turquoise", 
                 "Naive_CD8"="#FB9A99", "Naive_CD4"="#B2DF8A",             
                 "Basophils"="goldenrod2", "Monocytes_nonclassical"="sienna2", "Central_Memory_CD8"="#E31A1C", 
                 "Monocytes_intermediate"="sienna2", "B_Cells"="dodgerblue4",               
                 "Myeloid_DC"="lemonchiffon1", "Plasmacytoid_DC"="lemonchiffon1", "Naive_B_Cells"="dodgerblue2", 
                 "Central_Memory_CD4"="#33A02C", "Effector_Memory_CD4"="#33A02C",   
                 "TEMRA_CD4"="#33A02C")
pop.palette <- pop.palette[names(pop.palette) %in% node_df$node]
node_df$color <- "lightsteelblue"
node_df$color[match(names(pop.palette), node_df$node)] <- pop.palette


#edge df, for pop-specific genes
edge_df <- data.frame()
for (n in 1:length(gene.pop)){
  g <- gene.pop[[n]]
  if (length(g)>0){
    if (length(g)<MIN_GENES){
      edge_df <- rbind(edge_df,
                       data.frame(source=names(gene.pop)[n], target=g, weight=20))
    } else {
      edge_df <- rbind(edge_df,
                       data.frame(source=names(gene.pop)[n], 
                                  target=paste0(names(gene.pop)[n], "_gene"), 
                                  weight=20))
    }
  }
}

#edge df, for shared pop-pop genes
for (n in 1:nrow(df)){
  pop1 <- df$source[n]
  pop2 <- df$target[n]
  g <- intersect(net.group$gene[net.group$pop==pop1], net.group$gene[net.group$pop==pop2])
  if (length(g)<MIN_GENES){
    edge_df <- rbind(edge_df,
                     data.frame(source=pop1, target=g, weight=df$weights[n]),
                     data.frame(source=pop2, target=g, weight=df$weights[n]))
  } else {
    edge_node <- paste0(pop1, "_", pop2, "_edge_gene")
    edge_df <- rbind(edge_df,
                     data.frame(source=pop1, target=edge_node, weight=df$weights[n]),
                     data.frame(source=pop2, target=edge_node, weight=df$weights[n]))
    
  }
}            

net <- graph_from_data_frame(d=edge_df, vertices = node_df, directed=F)
layout <- layout_nicely(net)
plot(net, layout=layout, vertex.label.cex=.6, 
     vertex.label.color="black", edge.label.color="black", edge.label.cex=.8,
     label.family="Arial")


pdf("tmp.pdf", height=7, width=8)
plot(net, layout=layout, vertex.label.cex=.6, vertex.label=NA, 
     vertex.label.color="black", edge.label.color="black", edge.label.cex=.8,
     label.family="Arial")
dev.off()


# Map changes to the PBMC network ####
extract_beta <- function(g, variable_name, lm.df){
  
  df.beta <- data.frame()
  for (n in 1:length(g)){
    df.beta <- rbind(df.beta, 
                     lm.df[[ g[n] ]] %>% mutate(adj.pval=p.adjust(`Pr(>|t|)`, method="BH")) %>% 
                       filter(var==variable_name) %>% subset(select=c(Estimate, adj.pval)) %>% mutate(var=g[n]))
  }
  return(df.beta[, c("var", "Estimate", "adj.pval")])
}

THS_COR <- 0.7
MAX_EDGES <- 1e2
N_TOP_GENES <- 10
for (pop in names(gene.module)[c(3,7)]){
  g <- gene.module[[pop]]
  g.cor <- rcorr(x=rna.log[,g] %>% as.matrix(), type="spearman")
  pval.mat <- g.cor$P
  pval.mat[upper.tri(pval.mat)] <- NA
  g.cor.df <- reshape2::melt(pval.mat) %>% na.omit() %>% filter(Var1 != Var2) %>% 
    mutate(bonferroni=value*n()) %>% filter(bonferroni<1e3) %>% subset(select=-bonferroni)
  
  if (nrow(g.cor.df)>0){
    # select genes by centrality
    net <- graph_from_data_frame(d=g.cor.df, directed=F)
    net_btw <- betweenness(net, directed=F)
    #net_btw <- degree(net)
    #net_btw <- harmonic_centrality(net)

    # keep in network only selected genes and corresponing edges
    g.top <- names(net_btw %>% sort(decreasing=T))[1:min(length(net_btw), N_TOP_GENES)]
    g.cor.df <- g.cor.df %>% filter(Var1 %in% g.top | Var2 %in% g.top)
    g.cor.df <- g.cor.df[order(g.cor.df$value, decreasing=F)[1:MAX_EDGES],]
    
    #g.cor <- cor(rna.log[,g], method="spearman")
    #g.cor.df <- reshape2::melt(g.cor) %>% filter(abs(value)>THS_COR) %>% filter(Var1 != Var2)
    
    g.ass <- df.lm.all %>% filter(adj.pval<0.05) %>% filter(gene %in% g) %>% 
      group_by(gene) %>% summarise(var=group[which.min(adj.pval)]) %>% as.data.frame() 
    
    g.ass$color <- color_palette[g.ass$var]
    g.ass <- rbind(g.ass, data.frame(gene=setdiff(g, g.ass$gene), var="none", color="gray77"))
    
    g.ass <- g.ass %>% filter(gene %in% union(g.cor.df$Var1, g.cor.df$Var2))
    
    g.ass$btw <- net_btw[g.ass$gene] 
    g.ass$size <- 10*g.ass$btw / max(g.ass$btw) 
    g.ass$label <- ifelse(g.ass$gene %in% g.top, g.ass$gene, "")
    
    net <- graph_from_data_frame(d=g.cor.df, vertices=g.ass, directed=F)
    #layout <- layout.forceatlas2(net, plotstep = Inf)
    layout <- layout_with_gem(net)
    
    #cairo_pdf(paste0("net-", pop, ".pdf"), height=6, width=6)
    plot(net, layout=layout, 
         vertex.label.cex=.6, vertex.label.color="black", 
         edge.label.color="black", edge.label.cex=.8, edge.lty=2, edge.curved=F,
         main=pop) %>% print()
    #dev.off()
  }
}


# Stable vs unstable genes ####
g.stable <- rownames(df.CV)[df.CV$inter>df.CV$intra]
g.unstable <- rownames(df.CV)[df.CV$inter<df.CV$intra]
lm.out.stable <- lm.out.rna[g.stable]
lm.out.unstable <- lm.out.rna[g.unstable]

out.varExpl.stable <- plot_varExpl(lm.out.stable, Nmax=length(lm.out.stable))
out.varExpl.unstable <- plot_varExpl(lm.out.unstable, Nmax=length(lm.out.unstable))


df.varExpl.summary.stable <- out.varExpl.stable$df2
df.varExpl.summary.stable <- df.varExpl.summary.stable[order(df.varExpl.summary.stable$sumPct, decreasing=T),]

df.varExpl.summary.unstable <- out.varExpl.unstable$df2
df.varExpl.summary.unstable <- df.varExpl.summary.unstable[order(df.varExpl.summary.unstable$sumPct, decreasing=T),]


# Variance explained of stable and not stable pop ####
df.CV <- readRDS("df-CV-pop.RDS")

g.stable <- rownames(df.CV)[df.CV$inter>df.CV$intra]
g.unstable <- rownames(df.CV)[df.CV$inter<df.CV$intra]

common.samples.lifestyle <- intersect(rownames(clinical), rownames(lifestyle))

Y=rna.log[common.samples.lifestyle, g.stable] %>% scale() %>% as.data.frame()
lm.out.stable <- lm_fun(X=X, 
                        Y=Y, 
                        GWAS.summary=GWAS.rna.LD)

Y=rna.log[common.samples.lifestyle, g.unstable] %>% scale() %>% as.data.frame()
lm.out.unstable <- lm_fun(X=X, 
                          Y=Y, 
                          GWAS.summary=GWAS.rna.LD)

out.varExpl.stable <- plot_varExpl(lm.out.stable)
out.varExpl.unstable <- plot_varExpl(lm.out.unstable)




#
# select genes by hclust
dist_matrix <- as.dist(1 - g.cor) # Convert correlation to distance
hc <- hclust(dist_matrix, method = "ward.D")

num_clusters <- 50
clusters <- cutree(hc, k = num_clusters)

g <- sapply(unique(clusters), function(cluster) {
  cluster_genes <- names(clusters[clusters == cluster])
  cluster_genes[1] # Choose the first gene or any preferred criteria
})



plot_volcano <- function(df.out, variable_name){
  
  df.plot <- data.frame()
  for (n in 1:length(df.out)){
    print(paste0(n, "/", length(df.out)))
    df.plot <- rbind(df.plot,
                     df.out[[n]] %>% mutate(adj.pval=p.adjust(`Pr(>|t|)`, method="BH")) %>% filter(var==variable_name) %>% mutate(gene=names(df.out)[n]) )
  }
  
  df.plot$DE <- "no"
  df.plot$DE[df.plot$adj.pval<0.05 & df.plot$Estimate>0] <- "Upreg"
  df.plot$DE[df.plot$adj.pval<0.05 & df.plot$Estimate<0] <- "Downreg"
  
  df.plot$label <- df.plot$gene
  df.plot$label[df.plot$DE=="no"] <- NA
  color.DEG <- c(no="gray77", Upreg="tomato", Downreg="skyblue")
  p <- ggplot(df.plot %>% na.omit(), aes(x=Estimate, y=-log10(adj.pval), color=DE)) + 
    geom_hline(yintercept=-log10(0.05), linetype=2) + 
    geom_vline(xintercept=0, linetype=2) +
    geom_point(size=1,alpha=0.7) + 
    geom_text_repel(aes(label=label), size=3, color="black") + 
    scale_color_manual(values=color.DEG) +
    theme_classic2() + 
    xlab("beta") + ylab("-Log10(P)") +
    theme(text=element_text(size=10),
          legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    ggtitle(variable_name)
  
  return(list(p=p, df=df.plot))
  
}

volcano.BMI <- plot_volcano(lm.out.rna, "BMI")
volcano.CRP <- plot_volcano(lm.out.rna, "CRP")
volcano.Sex <- plot_volcano(lm.out.rna, "Gender")
volcano.fys <- plot_volcano(lm.out.rna, "PhysicalActivity")

p <- ggarrange(volcano.Sex$p, volcano.BMI$p, volcano.CRP$p, volcano.fys$p, nrow=2, ncol=2)
cairo_pdf("volcano-gene.pdf", width=6, height=6)
print(p)
dev.off()


# Variance explained of stable and unstable genes ####
df.CV <- readRDS("df-CV-gene.RDS")

#all genes
common.samples.lifestyle <- intersect(rownames(clinical), rownames(lifestyle)) %>% intersect(rownames(rna.log))
Y <- rna.log[common.samples.lifestyle, ] %>% scale() %>% as.data.frame()
lm.out.rna <- lm_fun(X=X, 
                     Y=Y, 
                     GWAS.summary=GWAS.rna.LD) 
df.varExpl.rna <- plot_varExpl(lm.out.rna, Nmax=ncol(rna.log))
df.varExpl.summary <- df.varExpl.rna$df2
df.varExpl.summary <- df.varExpl.summary[order(df.varExpl.summary$sumPct, decreasing=T),]

df.varExpl.rna <- df.varExpl.rna$df


# Map changes to the PBMC network ####
group.uniq <- df.lm.all$group %>% unique()
cell.gene <- readRDS("ILR-cytof.group-310125.RDS")$cell.gene.sub[,1:5]

group.n <- "Genetics"
gene.group.n <- df.lm.all %>% filter(adj.pval<0.05) %>% filter(group==group.n) %>% pull(gene)
net.group <- cell.gene %>% filter(gene %in% gene.group.n) 
net.group <- net.group[!net.group$pop %in% c("Effector_Memory_CD8", "ILCs"),] 

#select pairs with significant overlap
source('C:/Users/albze08/Desktop/postDoc/functions/my_fisher_test.R')
pop.uniq <- unique(net.group$pop)
k <- 1
df <- data.frame()
for (pop1 in pop.uniq){
  for (pop2 in pop.uniq){
    if (pop1 != pop2){
      if (!any(df$source==pop1 & df$target==pop2)){
        if (!any(df$target==pop1 & df$source==pop2)){
          
          pop1.g <- net.group$gene[net.group$pop==pop1]
          pop2.g <- net.group$gene[net.group$pop==pop2]
          f <- my_fisher_test(pop1.g, colnames(rna), pop2.g)
          df <- rbind(df, 
                      data.frame(source=pop1, target=pop2, pval=f$p_value, OR=f$odds_ratio, 
                                 n1=length(pop1.g), n2=length(pop2.g), n=f$n_genes))
          k <- k + 1
        }
      }
    }
  }
}
df$adj.pval <- p.adjust(df$pval, method="BH")
df$weights <- df$n
df_overlap_pop <- df
df <- df[df$adj.pval<0.05,]

#edge df
edge_df <- data.frame(source=net.group$pop, target=net.group$gene) %>% rbind(
  data.frame(source=df$source, target=df$target) )

#node df
node_df <- data.frame(node=unique(c(edge_df$source, edge_df$target)))
node_df$label <- ifelse(node_df$node %in% colnames(cytof.group), node_df$node, "")
node_df$size <- ifelse(node_df$node %in% colnames(cytof.group), 4, 2)

#color
pop.palette <- c("Monocytes_classical"="sienna2", "TEMRA_CD8"="#E31A1C", "NK_Cells"="turquoise", 
                 "Naive_CD8"="#FB9A99", "Naive_CD4"="#B2DF8A",             
                 "Basophils"="goldenrod2", "Monocytes_nonclassical"="sienna2", "Central_Memory_CD8"="#E31A1C", 
                 "Monocytes_intermediate"="sienna2", "B_Cells"="dodgerblue4",               
                 "Myeloid_DC"="lemonchiffon1", "Plasmacytoid_DC"="lemonchiffon1", "Naive_B_Cells"="dodgerblue2", 
                 "Central_Memory_CD4"="#33A02C", "Effector_Memory_CD4"="#33A02C",   
                 "TEMRA_CD4"="#33A02C")
pop.palette <- pop.palette[names(pop.palette) %in% node_df$node]
node_df$color <- "lightsteelblue"
node_df$color[match(names(pop.palette), node_df$node)] <- pop.palette

net <- graph_from_data_frame(d=edge_df, vertices = node_df, directed=F)
layout <- layout_with_gem(net)
plot(net, layout=layout, 
     vertex.label.cex=.6, vertex.label.color="black", 
     edge.label.color="black", edge.label.cex=.8)