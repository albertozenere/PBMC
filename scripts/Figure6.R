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
metadata$subject_id <- gsub("1-", "", metadata$subject)
metadata$id <- paste0(metadata$subject_id, ":", metadata$visit)
rownames(metadata) <- metadata$id

# Filter metadata
clinical <- metadata
rownames(clinical) <- metadata$id
clinical$Gender <- ifelse(clinical$Gender=="f", 0, 1)
clinical <- subset(clinical, select=c(Gender, Height, Cap_Gluc, SBP, DBP, Weight, Waist, Hip, Bioimp_fat, Bioimp_muscle, Bioimp_bone, PhysicalActivity, Age_at_Visit, BMI, 
                                      ProBNP, TNT, CRP, ALAT, GGT, HDL, Chol, LDL, TG, HbA1c, Urate, Gluc, CystC, Crea, ApoB.apoA1, ApoB, ApoA1))

clinical$PhysicalActivity[clinical$PhysicalActivity=="aldrig"] <- 0
clinical$PhysicalActivity[clinical$PhysicalActivity=="da_och_da"] <- 1
clinical$PhysicalActivity[clinical$PhysicalActivity=="1_2_ggr_vecka"] <- 2
clinical$PhysicalActivity[clinical$PhysicalActivity=="2_3_ggr_vecka"] <- 3
clinical$PhysicalActivity[clinical$PhysicalActivity=="mer_an_3_ggr_vecka"] <- 4
clinical$PhysicalActivity <- as.numeric(clinical$PhysicalActivity)


#RNA-seq S3WP
rna_s3wp <- read.table("data/wellness_PBMC_v16_norm.txt", sep="\t", header = T)


#Cytof S3WP
cytof <- read.table("data/original.cytof.txt", sep="\t", header = T)
rownames(cytof) <- cytof$SampleID
cytof <- subset(cytof, select=-SampleID)
rownames(cytof) <- gsub("_", ":", rownames(cytof))


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
  
#calculate number independent pop, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6693916/
# #pval.cytof <- 5e-8/7.788615 #CLR
# pval.cytof <- 5e-8/12.55736 #orig
# pval.cytof <- 5e-8/4.09191 #mean
# #pval.rna <- 5e-8/21.64967
# pval.rna <- 5e-8/15.96234 #mean

pval.rna <- 5e-8/ncol(rna.log)
pval.cytof <- 5e-8/ncol(cytof)

# gather cQTLs
gather_assoc <- function(folder_QTL, pval_loose=5e-8){
  
  QTL_files <- grep("linear", list.files(folder_QTL), value=T)
  df_summary <- data.frame()
  for (n in 1:length(QTL_files)){
    f <- read.table(paste0(folder_QTL, QTL_files[n]), header=F) 
    colnames(f) <- c("CHR", "BP", "SNP", paste0("V", 1:8), "beta", "V9", "V10", "P", "V11")
    f <- f[f$P<pval_loose,]
    if (nrow(f)>0){
      f$group <- gsub("out\\.", "", QTL_files[n])
      f$group <- gsub("\\.glm\\.linear", "", f$group)
      
      df_summary <- rbind(df_summary, f)
    }
  }
  df_summary$id <- paste0(df_summary$CHR, ":", df_summary$BP)
  df_summary$BP <- as.numeric(df_summary$BP)
  
  return(df_summary)
}

gather_conditional_assoc <- function(folder_QTL, pval_loose=5e-8){
  var <- grep("linear", list.files(paste0(folder_QTL, "conditional")), value=T)
  var <-  gsub("\\..*", "", var)
  QTL_files <- grep(paste(var, collapse="|"), list.files(folder_QTL), value=T)
  QTL_files <- QTL_files[grep("linear", QTL_files)]
  df_summary <- data.frame()
  for (n in 1:length(QTL_files)){
    f <- read.table(paste0(folder_QTL, QTL_files[n]), header=F) 
    colnames(f) <- c("CHR", "BP", "SNP", paste0("V", 1:8), "beta", "V9", "V10", "P", "V11")
    f <- f[f$P<pval_loose,]
    if (nrow(f)>0){
      f$group <- gsub("out\\.", "", QTL_files[n])
      f$group <- gsub("\\.glm\\.linear", "", f$group)
      
      df_summary <- rbind(df_summary, f)
    }
  }
  df_summary$family <- macro.anno.my[match(df_summary$group, names(macro.anno.my))]
  GWAS.cytof <- df_summary
  GWAS.cytof$id <- paste0(GWAS.cytof$CHR, ":", GWAS.cytof$BP)
  GWAS.cytof$BP <- as.numeric(GWAS.cytof$BP)
  
  return(GWAS.cytof)
}

#LD
gather_assoc_LD <- function(folder_QTL, pval.ths){
  QTL_files <- grep("clumps", list.files(paste0(folder_QTL, "LD")), value=T)
  df_summary.QTL.clump <- data.frame()
  for (n in 1:length(QTL_files)){
    var <- gsub("LD\\.out\\.", "", QTL_files[n])
    var <- gsub("\\.glm\\.linear\\.clumps", "", var)
    
    f <- read.table(paste0(folder_QTL, "LD/", QTL_files[n]), header=F)
    colnames(f) <- c("CHR", "BP", "SNP", "P", paste0("V", 1:6))
    f <- f[f$P<pval.ths,]
    if (nrow(f)>0){
      f$group <- var
      df_summary.QTL.clump <- rbind(df_summary.QTL.clump, f[which.min(f$P),c("SNP", "P", "group")])
    }
    if (nrow(f)>1){
      conditional_file <- grep(paste0(var, "\\."), list.files(paste0(folder_QTL, "conditional")), value=T)
      conditional_file <- grep("linear", conditional_file, value=T)
      if (length(conditional_file)==1){
        f_conditional <- read.table(paste0(folder_QTL, "conditional/", conditional_file), header=F)
        colnames(f_conditional) <- c("CHR", "BP", "SNP", paste0("V", 1:11), "P", "ERRCODE")
        f_conditional$group <- var
        
        idx <- which(f_conditional$P<0.1) #this is the conditional p-value, i.e. when doing linear model that includes the dominant variant
        f_conditional <- f_conditional[idx, ]
        f_conditional$P <- f$P[match(f_conditional$SNP, f$SNP)] #store the actual p-value from GWAS, not the conditional
        df_summary.QTL.clump <- rbind(df_summary.QTL.clump, f_conditional[,c("SNP", "P", "group")])
      } else {
        print(paste0("[WARNING] Number of files found for ", var, " : ", length(conditional_file)))
      }
    }
  }
  GWAS.cytof.LD <- df_summary.QTL.clump
  
  str <- strsplit(GWAS.cytof.LD$SNP, "\\:") %>% unlist()
  GWAS.cytof.LD$CHR <- str[seq(1,length(str),4)] %>% as.numeric()
  GWAS.cytof.LD$BP <- str[seq(2,length(str),4)] %>% as.numeric()
  GWAS.cytof.LD$id <- paste0(GWAS.cytof.LD$CHR, ":", GWAS.cytof.LD$BP)
  
  GWAS.cytof.LD <- GWAS.cytof.LD %>% arrange(CHR,BP,P)
  return(GWAS.cytof.LD)
}

closest_snp_fun <- function(df1, df2){
  df.out <- data.frame()
  for (n in 1:nrow(df1)){
    idx <- which(df2$CHR == df1$CHR[n])
    d <- abs(df1$BP[n] - df2$BP[idx])
    
    dist.min <- min(d) #to account for multiple matches with the same distance
    idx.min.vec <- idx[which(d==dist.min)]
    
    for (idx.min in idx.min.vec){
      df.out <- rbind(df.out, 
                      data.frame(df1[n,],
                                 idx=idx.min,
                                 distance=dist.min,
                                 CHR_other=df2$CHR[idx.min],
                                 BP_other=df2$BP[idx.min]) )
    }
  }
  return(df.out %>% arrange(distance))
}

GWAS.lmm <- gather_assoc("cQTL-mean-53/", pval_loose = 1e-6)
GWAS.lmm$family <- macro.anno.my[match(GWAS.lmm$group, names(macro.anno.my))]

GWAS.lmm.LD <- gather_assoc_LD("cQTL-mean-53/", pval.cytof)
GWAS.lmm.LD$family <- macro.anno.my[match(GWAS.lmm.LD$group, names(macro.anno.my))]
GWAS.lmm.LD$beta <- GWAS.lmm$beta[match(GWAS.lmm.LD$SNP, GWAS.lmm$SNP)]

SNP_LD <- GWAS.lmm$SNP %>% unique()

GWAS.lmm.eqtl <- gather_assoc("eQTL-mean-8K/", pval_loose = 5e-8)
GWAS.lmm.LD.eqtl <- gather_assoc_LD("eQTL-mean-8K/", pval.rna)
GWAS.lmm.LD.eqtl$beta <- GWAS.lmm.eqtl$beta[match(GWAS.lmm.LD.eqtl$SNP, GWAS.lmm.eqtl$SNP)]

GWAS.lmm.eqtl.filt <- GWAS.lmm.eqtl %>% filter(P<5e-8)
GWAS.shared <- closest_snp_fun(GWAS.lmm.LD %>% filter(beta>0), GWAS.lmm.eqtl.filt) 
GWAS.shared$gene <- GWAS.lmm.eqtl.filt$group[GWAS.shared$idx]
GWAS.shared$gene[GWAS.shared$distance==0] %>% unique()

# add gene column
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=115)
genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol",
                 "chromosome_name", "start_position", "end_position", "strand"),
  mart = mart
)

gene.anno <- genes %>% filter(chromosome_name %in% 1:22) %>%
  filter(hgnc_symbol!="") %>%
  mutate(CHR=chromosome_name, BP=(start_position+end_position)/2)

snp_gene_anno <- closest_snp_fun(GWAS.lmm.LD, gene.anno)
snp_gene_anno$gene <- gene.anno$hgnc_symbol[snp_gene_anno$idx]

GWAS.lmm.LD$gene <- snp_gene_anno$gene[match(GWAS.lmm.LD$SNP, snp_gene_anno$SNP)]


GWAS.lmm$gene <- snp_gene_anno$gene[match(GWAS.lmm$SNP, snp_gene_anno$SNP)]
GWAS.lmm$gene[GWAS.lmm$P > pval.cytof] <- NA

GWAS.lmm.eqtl$gene <- ifelse(GWAS.lmm.eqtl$P<pval.rna, GWAS.lmm.eqtl$group, NA)

#debug
GWAS.1 <- gather_assoc("eQTL-mean/", pval_loose = 5e-8)
GWAS.2 <- gather_assoc("eQTL-mean/", pval_loose = 1e-7)

GWAS.shared.1 <- closest_snp_fun(GWAS.lmm.LD, GWAS.1) 
GWAS.shared.1$gene <- GWAS.1$group[GWAS.shared.1$idx]

GWAS.shared.2 <- closest_snp_fun(GWAS.lmm.LD, GWAS.2) 
GWAS.shared.2$gene <- GWAS.2$group[GWAS.shared.2$idx]


# Supplementary Table 4 ####
df <- GWAS.lmm.LD[, c("SNP", "beta", "P", "group", "family", "gene")]
colnames(df) <- c("SNP", "Coefficient", "P-value", "Population", "Major population", "Closest gene")

df2 <- GWAS.lmm.LD.eqtl[, c("SNP", "beta", "P", "group")]
colnames(df2) <- c("SNP", "Coefficient", "P-value", "Gene")

wb <- createWorkbook()
addWorksheet(wb, "Immune frequencies GWAS")
addWorksheet(wb, "RNA-seq GWAS")
writeData(wb, "Immune frequencies GWAS", df)
writeData(wb, "RNA-seq GWAS", df2)

saveWorkbook(wb, "Supplementary Table 4.xlsx", overwrite = TRUE)


#Manhattan plot function ####
manhattan_plot <- function(df, pval.ths=5e-8, MAX_SNP=5e4){
  
  #Prepare GWAS CyTOF
  df_plot <- df[order(df$P, decreasing=F), c("CHR", "BP", "P", "group", "family", "SNP", "gene")]
  df_plot <- df_plot[!df_plot$CHR %in% c("PAR1", "PAR2"),]
  df_plot$CHR <- as.numeric(df_plot$CHR)
  
  #merge
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

  p <- ggplot(df_plot, aes(x=BP_cum, y=-log10(P), fill=color)) +
    geom_point(shape=21) + 
    geom_hline(yintercept=-log10(pval.ths), linetype=2) +
    theme_bw() +
    geom_text(aes(label=gene), cex=2, color="black") +
    scale_x_continuous(label = chr.df$CHR, breaks = chr.df$BP_mean) +
    scale_fill_manual(values=c(even="gray77", odd="gray23")) + 
    scale_size_continuous(range = c(0.5,3)) + 
    labs(x = NULL, y = "-log<sub>10</sub>(p)") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=0, size=8, vjust=0.5),
          legend.position="none",
          text=element_text(size=12)) +
    ylab("-Log10 P-value")
  #p  
  
  return(list(df=df_plot, p=p, chr.df=chr.df))
}

# Manhattan of specific pop ####
manhattan.Bcells <- manhattan_plot(GWAS.lmm[grep("B_Cells", GWAS.lmm$family),], pval.ths = pval.cytof)$p + ggtitle("B cells") + theme(plot.title = element_text(hjust = 0.5))
manhattan.CD8 <- manhattan_plot(GWAS.lmm[grep("CD8", GWAS.lmm$family),], pval.ths = pval.cytof)$p + ggtitle("CD8+ T cells") + theme(plot.title = element_text(hjust = 0.5))

manh

pdf("manhattan_specific_221025.pdf", width=7, height=8)
ggarrange(manhattan.Bcells, manhattan.CD8, ncol=1)
dev.off()


# Combined frequency-expression Manhattan ####
GWAS.both <- rbind(GWAS.lmm %>% mutate(omic="CyTOF"),
                   GWAS.lmm.eqtl %>% mutate(family="", omic="RNA")) %>% arrange(P)
idx.cytof <- which(GWAS.both$omic=="CyTOF")
idx.rna <- which(GWAS.both$omic=="RNA")

manhattan.both <- manhattan_plot(GWAS.both, MAX_SNP=Inf)
manhattan.cytof <- manhattan.both$df[idx.cytof,] %>% arrange(P)

manhattan.rna <- manhattan.both$df[idx.rna,] %>% arrange(P)
manhattan.rna <- manhattan.rna[1:3e4,]

chr.df <- manhattan.both$chr.df

# CyTOF
manhattan.cytof$gene <- ifelse(manhattan.cytof$P<pval.cytof, manhattan.cytof$gene, NA)

pop.palette <- c("Monocytes_classical"="sienna2", "TEMRA_CD8"="#E31A1C", "NK_Cells"="turquoise", 
                 "Naive_CD8"="#FB9A99", "Naive_CD4"="#B2DF8A",             
                 "Basophils"="goldenrod2", "Monocytes_nonclassical"="sienna2", "Central_Memory_CD8"="#E31A1C", 
                 "Monocytes_intermediate"="sienna2", "B_Cells"="dodgerblue4",               
                 "Myeloid_DC"="lemonchiffon1", "Plasmacytoid_DC"="lemonchiffon1", "Naive_B_Cells"="dodgerblue2", 
                 "Central_Memory_CD4"="#33A02C", "Effector_Memory_CD4"="#33A02C",   
                 "TEMRA_CD4"="#33A02C", "Effector_Memory_CD8"="#E31A1C", "ILCs"="gray78")
manhattan.cytof$color <- ifelse(manhattan.cytof$P<pval.cytof, manhattan.cytof$family, manhattan.cytof$color)

p1 <- ggplot(manhattan.cytof, aes(x=BP_cum, y=-log10(P), fill=color)) +
  geom_point(shape=21) + 
  geom_hline(yintercept=-log10(pval.cytof), linetype=2) +
  theme_bw() +
  geom_text(aes(label=gene), cex=2, color="black") +
  scale_x_continuous(label = chr.df$CHR, breaks = chr.df$BP_mean) +
  scale_fill_manual(values=c(pop.palette, even="gray77", odd="gray23")) + 
  scale_size_continuous(range = c(0.5,3)) + 
  labs(x = NULL, y = "-log<sub>10</sub>(p)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0, size=8, vjust=0.5),
        legend.position="none",
        text=element_text(size=12)) +
  ylab("-Log10 P-value")

#RNA-seq
df_min_p <- manhattan.rna %>%
  group_by(gene) %>%
  slice_min(P, n = 1, with_ties = FALSE) %>%
  ungroup()

p2 <- ggplot(manhattan.rna, aes(x=BP_cum, y=-log10(P), fill=color)) +
  geom_point(shape=21) + 
  geom_hline(yintercept=-log10(pval.rna), linetype=2) +
  theme_bw() +
  geom_text(data=df_min_p, aes(label=gene), cex=2, color="black") +
  scale_x_continuous(label = chr.df$CHR, breaks = chr.df$BP_mean) +
  scale_fill_manual(values=c(even="gray77", odd="gray23")) + 
  scale_size_continuous(range = c(0.5,3)) + 
  labs(x = NULL, y = "-log<sub>10</sub>(p)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0, size=8, vjust=0.5),
        legend.position="none",
        text=element_text(size=12)) +
  ylab("-Log10 P-value")

ggarrange(p1, p2, ncol=1)


#Plot
pdf("manhattan_CyTOF.pdf", width=7, height=4)
p1
dev.off()



# Load Genotype ####
library("vcfR")
vcf <- read.vcfR("C:/Users/albze08/Desktop/postDoc/wellness/Wellness/Data/WGS/HPA_101_qc_maf0.05.vcf")
snp <- vcf@fix %>% as.data.frame()
geno <- vcf@gt %>% as.data.frame()
rm(vcf)
geno <- geno[,2:ncol(geno)]
colnames(geno) <- anno$Subject[match(colnames(geno), paste0("P", anno$Barcode))]

stopifnot(all(SNP_LD %in% snp$ID))

idx <- match(SNP_LD, snp$ID) 
snp <- snp[idx,]
geno <- geno[idx,]


# Transform geno in a matrix ####
geno.X <- matrix(NA, nrow(geno), ncol(geno))
rownames(geno.X) <- rownames(geno)
colnames(geno.X) <- colnames(geno)
for (n in 1:nrow(geno.X)){
  vec <- geno[n,]
  vec[vec=="0/0"] <- 0
  vec[vec=="0/1"] <- 1
  vec[vec=="1/1"] <- 2
  vec <- as.numeric(vec)
  geno.X[n,] <- vec
}


# Contribution of genetic and non genetic variables to immune frequencies ####
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



var.group.SNP <- data.frame(var=unique(SNP_LD),
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
    geom_bar(stat="identity", position="stack", color="black") +
    scale_fill_manual(values=color_palette) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0)) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Variance explained") + xlab("")
  
  
  #Plot total variance explained
  df$meanPct <- df$sumPct/length(unique(df.plot$x))
  
  
  return(list(p=p, df=df.plot, df2=df))  
}


color_palette <- c("darkseagreen",	"#A30059",	"#9e5a28ff",		
                   "tomato",
                   "#4a6fe3",	"lightblue",	"steelblue" ,
                   "bisque",	"#1CE6FF",
                   "goldenrod")
names(color_palette) <- c("Genetics", "Sex", "Age", 
                          "Blood pressure", 
                          "Body composition", "Lipid profile","Glucose homeostasis",
                          "Organ biomarker", "CRP",
                          "Lifestyle")


clinical <- metadata %>%
  subset(select=-c(SCAPISorIGT_id, subject, Study, Visitdate, id, visit, Birthdate, subject_id,
                   Hb, WBC, Plt, RBC, Hct, MCV, MCH, MCHC, Neut, Lymph, Mono, Eos, Baso,
                   Housing_change, Housing_current, Address_change, MaritalStatus_change, MaritalStatus_current, ShareHousehold_change, 
                   ShareHousehold_current, Employment_change, Employment_current, Tobacco_change, Tobacco_current, PerceivedHealth,       
                   Stress, PhysicalActivity, SedentaryTime_hours, SedentaryTime_minutes, SedentaryTime_unknown, TravelAbroad,          
                   TravelAbroad_country, Animals, Animals_type,
                   Calcitriol, Calcidiol, Health_Other, Common_cold_influenzae, NSAID_painmed, Bp_med, Lipid_med,             
                   Antibiotics_med, Diab_med, Med_details, Smoking, Smoking_hours)) %>%
  mutate(Gender=ifelse(Gender=="m",1,0))


#major pop
GWAS.cytof.group.LD <- GWAS.lmm.LD
for (g in names(pop.group)){
  GWAS.cytof.group.LD$group[GWAS.lmm.LD$group %in% pop.group[[g]]] <- g
}

X <- clinical[rownames(cytof.group),] %>% scale() %>% as.data.frame()
Y <- cytof.group %>% scale() %>% as.data.frame()
lm.cytof.out <- lm_fun(X=X, 
                       Y=Y, 
                       GWAS.summary=GWAS.cytof.group.LD)

out.varExpl <- plot_varExpl(lm.cytof.out)
out.varExpl$p

pdf("anova.pdf", width=4, height=3)
out.varExpl$p
dev.off()

df <- out.varExpl$df
df <- df %>% arrange(group, desc(PctExp))

#specific pop
X <- clinical[rownames(cytof),] %>% scale() %>% as.data.frame()
Y <- cytof %>% scale() %>% as.data.frame()
lm.cytof.all <- lm_fun(X=X, 
                       Y=Y, 
                       GWAS.summary=GWAS.lmm.LD)
#saveRDS(lm.cytof.all, "lm.cytof.all.RDS")

out.varExpl <- plot_varExpl(lm.cytof.all)

df.plot <- data.frame(out.varExpl$df)
df.plot <- merge(df.plot, macro.anno.my, by.x="x", by.y="row.names")
colnames(df.plot)[ncol(df.plot)] <- "family"

pop.palette <- c("NK_cells"="turquoise",  "monocytes"="sienna2", "Dendritic_cells"="lemonchiffon1", "basophils"="goldenrod2", "Other"="gray78", 
                 "Naive_B_cells"="dodgerblue2", "B_cells"="dodgerblue4", "CD4pos_naive" ="#B2DF8A",   
                 "CD8pos_naive"="#FB9A99", "CD4pos"="#33A02C", "CD8pos"="#E31A1C")
color_palette <- c(color_palette, pop.palette)


p <- ggplot(df.plot, aes(x=x, y=PctExp, fill=group)) +
  geom_bar(stat="identity", position="stack", color="black") +
  geom_tile(aes(x=x, y=-3, height=3, fill=family)) +
  scale_fill_manual(values=color_palette) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0)) +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Variance explained") + xlab("") + ggtitle("Immune frequencies")
p




# Plot examples of genotype---frequency associations  ####

plot_SNP_frequency <- function(snp_id, pop, major.pop=NULL){
  idx.snp <- which(snp$ID==snp_id)
  if (is.null(major.pop)){
    df.plot <- data.frame(sample=rownames(cytof.nonnegative), value=dfcytof.nonnegative[,pop], omic="CyTOF")
    pop.name <- pop
  } else {
    df.plot <- data.frame(sample=rownames(cytof.nonnegative), value=cytof.group[,major.pop], omic="CyTOF")
    pop.name <- major.pop
  }
  df.plot$rank <- 100*(df.plot$value %>% rank())/nrow(df.plot)
  df.plot$ind <- gsub("\\:.*", "", df.plot$sample)
  df.plot$visit <- gsub(".*\\:", "", df.plot$sample)
  df.plot$geno <- geno[idx.snp, df.plot$ind] %>% as.character()
  df.plot <- na.omit(df.plot)
  df.plot$geno <- ifelse(df.plot$geno =="0/0", 0, 1) %>% as.character()
  
  p1 <- ggplot(df.plot, aes(x=geno, y=value)) +
    geom_quasirandom(dodge.width=0.8, shape=21, size=2, color="black") +  
    geom_boxplot(fill=NA, outlier.shape=NA) + 
    xlab(snp_id) + ylab("Percentile in the cohort") +
    theme_classic() + theme(legend.position="none", plot.title = element_text(size=10))
  
  p2 <- ggplot(df.plot, aes(x=visit, y=value, color=geno, fill=geno)) + 
    geom_line(aes(group=ind), alpha=0.5) +
    geom_point(data=df.plot %>% na.omit() %>% filter(geno=="0/0"), size=3, shape=21, color="black", alpha=0.8) +
    geom_point(data=df.plot %>% na.omit() %>% filter(geno=="0/1"), size=3, shape=21, color="black", alpha=0.8) +
    geom_point(data=df.plot %>% na.omit() %>% filter(geno=="1/1"), size=3, shape=21, color="black", alpha=0.8) +
    theme_classic() + theme(legend.position = "bottom") +
    ylab("Percentile in the cohort") + ggtitle(pop.name) + theme(plot.title = element_text(hjust = 0.5))
  
  df.plot.summ <- df.plot %>% group_by(geno, visit) %>% 
    summarise(min_value=min(value), max_value=max(value), median_value=median(value), low_05=quantile(value, 0.025), high_05=quantile(value, 0.975)) %>% 
    as.data.frame()
  
  p3 <- ggplot(df.plot.summ, aes(x=visit, fill=geno)) + 
    geom_line(data=df.plot, aes(x=visit, y=value, group=ind, color=geno), alpha=0.5) +
    geom_ribbon(aes(ymin=min_value, ymax=max_value)) +
    geom_point(size=3, aes(y=median_value), shape=21) +
    theme_classic() + theme(legend.position="bottom") +
    ylab("Percentile in the cohort")
  p3
  
  p4 <- ggplot(df.plot.summ, aes(x = visit, fill = geno, group = geno)) +
    geom_ribbon(aes(ymin = low_05, ymax = high_05), alpha = 0.3)  +
    geom_line(aes(group=geno, y=median_value, color=geno)) +
    geom_point(size=4, aes(y=median_value), shape=21) +
    theme_classic() +
    theme(legend.position = "bottom") +
    ylab("Percentile in the cohort") + xlab(snp_id) +
    ggtitle(pop.name)
  p4
  
  return(list(df=df.plot, p1=p1, p2=p2, p3=p3, p4=p4))
}

GWAS.lmm.LD <- GWAS.lmm.LD %>% arrange(desc(beta))

GWAS.lmm.LD.B <- GWAS.lmm.LD %>% filter(family=="B_Cells") %>% arrange(beta)
pdf(paste0("pdfs/2610-Bcells-major", ".pdf"), height=3, width=4)
for (n in 1:nrow(GWAS.lmm.LD.B)){
  
  snp_n <- GWAS.lmm.LD.B$SNP[n]
  pop_n <- GWAS.lmm.LD.B$group[n]
  
  out <- plot_SNP_frequency(snp_n, pop_n, major.pop=GWAS.lmm.LD.B$family[n])
  #out <- plot_SNP_frequency(snp_n, pop_n)
  p1 <- out$p1
  p2 <- out$p3
  p4 <- out$p4
  print(p4)
}
dev.off()

pdf(paste0("pdfs/2610-CD8", ".pdf"), height=3, width=4)
GWAS.lmm.pop <- GWAS.lmm.LD[grep("CD8", GWAS.lmm.LD$group), ]
for (n in 1:nrow(GWAS.lmm.pop)){
  
  snp_n <- GWAS.lmm.pop$SNP[n]
  pop_n <- GWAS.lmm.pop$group[n]
  
  out <- plot_SNP_frequency(snp_n, pop_n, major.pop=GWAS.lmm.pop$family[n])
  p1 <- out$p1
  p2 <- out$p3
  p4 <- out$p4
  print(p4)
}
dev.off()


# Genetic associations with clustering ####


# GWAS on the best match ####
#Frequency
snp_ENTPD1 <- GWAS.lmm.LD %>% filter(gene=="ENTPD1-AS1") %>% pull(SNP)
idx.snp <- which(snp$ID == snp_ENTPD1[1])
pop_ENTPD1 <- GWAS.lmm.LD %>% filter(gene=="ENTPD1-AS1") %>% pull(group)
family_ENTPD1 <- GWAS.lmm.LD %>% filter(gene=="ENTPD1-AS1") %>% pull(family)
p.freq <- plot_SNP_frequency(snp_ENTPD1[1], pop_ENTPD1[1], major.pop = family_ENTPD1[1])$p4

pdf(paste0("ENTPD1", ".pdf"), height=3, width=4)
p.freq
dev.off()

#Expression
df.plot <- data.frame(sample=rownames(rna.log), value=rna.log[,"ENTPD6"])
df.plot$rank <- 100*(df.plot$value %>% rank())/nrow(df.plot)
df.plot$ind <- gsub("\\:.*", "", df.plot$sample)
df.plot$visit <- gsub(".*\\:", "", df.plot$sample)
df.plot$geno <- geno[idx.snp, df.plot$ind] %>% as.character()
df.plot <- na.omit(df.plot)
df.plot$geno <- ifelse(df.plot$geno =="0/0", 0, 1) %>% as.character()

df.plot.summ <- df.plot %>% group_by(geno, visit) %>% 
  summarise(min_value=min(value), max_value=max(value), median_value=median(value), low_05=quantile(value, 0.025), high_05=quantile(value, 0.975)) %>% 
  as.data.frame()

p.expr <- ggplot(df.plot.summ, aes(x = visit, fill = geno, group = geno)) +
  geom_ribbon(aes(ymin = low_05, ymax = high_05), alpha = 0.3)  +
  geom_line(aes(group=geno, y=median_value, color=geno)) +
  geom_point(size=4, aes(y=median_value), shape=21) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ylab("Percentile in the cohort") + xlab(snp_ENTPD1[1]) 
p.expr




# MAF vs p-value ####
GWAS.lmm.LD$MAF <- NA
for (n in 1:nrow(GWAS.lmm.LD)){
  idx <- which(snp$ID == GWAS.lmm.LD$SNP[n])[1]
  str <- strsplit(geno[idx,] %>% as.character() %>% na.omit(), "/") %>% unlist()
  allel.freq <- str %>% table()
  GWAS.lmm.LD$MAF[n] <- sum ((allel.freq %>% names() %>% as.numeric()) * (allel.freq %>% as.numeric()))/length(str)
}

df.plot <- GWAS.lmm.LD %>% filter(family %in% c("B_Cells", "Central_Memory_CD8"))
family.ord <- GWAS.lmm.LD %>% group_by(family) %>% summarise(mean_MAF=median(MAF, na.rm=T)) %>% 
  arrange(desc(mean_MAF)) %>%
  pull(family)
df.plot$family <- factor(df.plot$family, levels=family.ord)

label.n.SNP <- df.plot %>% group_by(family) %>% summarise(n=n())
df.plot <- df.plot %>% arrange(family, MAF)
p <- ggplot(df.plot, aes(x=family, y=MAF)) +
  geom_boxplot(width=0.3, outlier.shape=NA) +
  geom_quasirandom(alpha=0.5, shape=21, fill="gray77", color="black") +
  geom_text(data=label.n.SNP, aes(x=family, label=n, y=0.2)) +
  theme_classic()  +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
pdf("MAF-pop.pdf", height=4, width=3)
print(p)
dev.off()


# B cell freq based on PGS ####
major.pop <- "B_Cells"
df.snp <- GWAS.lmm.LD %>% filter(family == major.pop)
df.snp$beta_weighted <- NA

pop.spec <- unique(df.snp$group)
for (pop in pop.spec){
  idx <- which(df.snp$group == pop)
  relative.freq <- mean(cytof.nonnegative[,pop])/mean(cytof.group[,major.pop])
  
  df.snp$beta_weighted[idx] <- df.snp$beta[idx]*relative.freq
}

idx.geno <- match(df.snp$SNP, snp$ID)
PRS <- apply(geno.X[idx.geno,], 2, function(x){sum(x*df.snp$beta_weighted, na.rm=T)}) # scale by relative freq of that specific pop?

pheno <- data.frame(y=cytof.group[,major.pop], sample=rownames(cytof.group))
pheno$ind <- gsub("\\:.*", "", pheno$sample)
pheno$PRS <- PRS[pheno$ind]
pheno$PRS_rank <- rank(pheno$PRS)
pheno$PRS_norm <- scale(pheno$PRS) %>% as.numeric()

ind.uniq <- unique(pheno$ind)
color.UMAP <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] %>% sample(length(ind.uniq))
color.UMAP <- setNames(color.UMAP, ind.uniq)

p <- ggplot(pheno, aes(x=PRS, y=y)) +
  geom_point(aes(fill=ind), shape=21, size=2) +
  geom_smooth(method="lm", color="red", size=1) +
  geom_text(label=paste0("Spearman's r= ", cor(pheno$PRS, pheno$y, method="spearman") %>% round(digits=2)), x=5, y=0.03) +
  scale_fill_manual(values=color.UMAP) +
  theme_classic() + theme(legend.position = "none") +
  xlab("Polygenic score") + ylab("Memory B cells frequency")
p

pdf("pdfs/PGS-memory B-mean.pdf", width=3, height=3)
print(p)
dev.off()

cor(pheno$y, pheno$PRS_norm, method="spearman")

#PGS based on cluster
pheno$cluster <- primary_cluster$primary_cluster[match(pheno$sample, primary_cluster$sample)]

comparisons <- list(c(1,2), c(1,3), c(2,3))
p <- ggplot(pheno %>% na.omit(), aes(x=cluster, y=PRS_rank)) +
  geom_boxplot() +
  geom_quasirandom(shape=21, aes(fill=cluster), size=2) +
  scale_fill_manual(values=c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73")) + 
  stat_compare_means(method="wilcox.test", comparisons=comparisons, size=3) +
  theme_classic()
p
pdf("PGS-cluster.pdf", width=4, height=3)
print(p)
dev.off()
  
#Intersect with B unit genes
idx.B <- match(GWAS.lmm.LD$SNP[GWAS.lmm.LD$family=="B_Cells"], snp$ID)
geno.X.B <- geno.X[idx.B,]


# CD4:CD8 freq based on PGS ####


calculate_PGS <- function(df.snp, tot.freq){
  df.snp$beta_weighted <- NA
  
  pop.spec <- unique(df.snp$group)
  for (pop in pop.spec){
    pop.freq <- mean(cytof.nonnegative[,pop])
    tot.mean <- tot.freq %>% mean()
    relative.freq <- pop.freq/tot.mean
    
    idx <- which(df.snp$group == pop)
    df.snp$beta_weighted[idx] <- df.snp$beta[idx]*relative.freq
  }
  
  idx.geno <- match(df.snp$SNP, snp$ID)
  PRS <- apply(geno.X[idx.geno,], 2, function(x){sum(x*df.snp$beta_weighted, na.rm=T)}) # scale by relative freq of that specific pop?
  
  pheno <- data.frame(y=tot.freq, sample=rownames(cytof.group))
  pheno$ind <- gsub("\\:.*", "", pheno$sample)
  pheno$PRS <- PRS[pheno$ind]
  pheno$PRS_rank <- rank(pheno$PRS)
  pheno$PRS_norm <- scale(pheno$PRS) %>% as.numeric()
  
  ind.uniq <- unique(pheno$ind)
  color.UMAP <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] %>% sample(length(ind.uniq))
  color.UMAP <- setNames(color.UMAP, ind.uniq)
  
  p <- ggplot(pheno, aes(x=PRS, y=y)) +
    geom_point(aes(fill=ind), shape=21, size=2) +
    geom_smooth(method="lm", color="red", size=1) +
    geom_text(label=paste0("Spearman's r= ", cor(pheno$PRS, pheno$y, method="spearman") %>% round(digits=2)), x=5, y=0.03) +
    scale_fill_manual(values=color.UMAP) +
    theme_classic() + theme(legend.position = "none") +
    xlab("Polygenic score") + ylab("CD8+ cells frequency")
  p
  
  return(list(p=p, pheno=pheno))
}

df <- GWAS.lmm.LD[grep("Central_Memory_CD8", GWAS.lmm.LD$family), ]
tot.CD8 <- cytof.group[,grep("Central_Memory_CD8", colnames(cytof.group))]
calculate_PGS(df, tot.CD8)$p
pheno <- calculate_PGS(df, tot.CD8)$pheno

df <- GWAS.lmm.LD[grep("CD4", GWAS.lmm.LD$family), ]
tot.CD4 <- cytof.group[,grep("CD4", colnames(cytof.group))] %>% rowSums()
calculate_PGS(df, tot.CD4)$p



pheno$cluster <- primary_cluster$primary_cluster[match(pheno$sample, primary_cluster$sample)]

comparisons <- list(c(1,2), c(1,3), c(2,3))
ggplot(pheno %>% na.omit(), aes(x=cluster, y=PRS_rank)) +
  geom_boxplot() +
  stat_compare_means(method="wilcox.test", comparisons=comparisons, size=3) +
  theme_classic()


pdf("pdfs/PGS-CD8-mean.pdf", width=3, height=3)
print(p)
dev.off()

cor(pheno$y, pheno$PRS_norm, method="spearman")

#PGS based on cluster
pheno$cluster <- primary_cluster$primary_cluster[match(pheno$sample, primary_cluster$sample)]

comparisons <- list(c(1,2), c(1,3), c(2,3))
ggplot(pheno %>% na.omit(), aes(x=cluster, y=PRS_rank)) +
  geom_boxplot() +
  stat_compare_means(method="wilcox.test", comparisons=comparisons, size=3) +
  theme_classic()


# CD8 freq based on PGS ####
major.pop <- "Central_Memory_CD8"
df.snp <- GWAS.lmm.LD %>% filter(family == major.pop)
df.snp$beta_weighted <- NA

pop.spec <- unique(df.snp$group)
for (pop in pop.spec){
  idx <- which(df.snp$group == pop)
  relative.freq <- mean(cytof.nonnegative[,pop])/mean(cytof.group[,major.pop])
  
  df.snp$beta_weighted[idx] <- df.snp$beta[idx]*relative.freq
}

idx.geno <- match(df.snp$SNP, snp$ID)
PRS <- apply(geno.X[idx.geno,], 2, function(x){sum(x*df.snp$beta_weighted, na.rm=T)}) # scale by relative freq of that specific pop?

pheno <- data.frame(y=cytof.group[,major.pop], sample=rownames(cytof.group))
pheno$ind <- gsub("\\:.*", "", pheno$sample)
pheno$PRS <- PRS[pheno$ind]
pheno$PRS_rank <- rank(pheno$PRS)
pheno$PRS_norm <- scale(pheno$PRS) %>% as.numeric()

ind.uniq <- unique(pheno$ind)
color.UMAP <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] %>% sample(length(ind.uniq))
color.UMAP <- setNames(color.UMAP, ind.uniq)

p <- ggplot(pheno, aes(x=PRS, y=y)) +
  geom_point(aes(fill=ind), shape=21, size=2) +
  geom_smooth(method="lm", color="red", size=1) +
  geom_text(label=paste0("Spearman's r= ", cor(pheno$PRS, pheno$y, method="spearman") %>% round(digits=2)), x=5, y=0.03) +
  scale_fill_manual(values=color.UMAP) +
  theme_classic() + theme(legend.position = "none") +
  xlab("Polygenic score") + ylab("Central Memory CD8")
p

pdf("pdfs/PGS-Central CD8-mean.pdf", width=3, height=3)
print(p)
dev.off()

cor(pheno$y, pheno$PRS_norm, method="spearman")

#PGS based on cluster
pheno$cluster <- primary_cluster$primary_cluster[match(pheno$sample, primary_cluster$sample)]

comparisons <- list(c(1,2), c(1,3), c(2,3))
ggplot(pheno %>% na.omit(), aes(x=cluster, y=PRS_rank)) +
  geom_boxplot() +
  stat_compare_means(method="wilcox.test", comparisons=comparisons, size=3) +
  theme_classic()

# Load network ####
cell.gene <- readRDS("Ensemble_ILR-306025.RDS") %>% bind_rows() %>% 
  filter(estimate>0) %>% 
  #mutate(pval.bon=pval*n()) %>% filter(pval.bon<0.05) %>% 
  mutate(pval.BH=p.adjust(pval, method="BH")) %>% filter(pval.BH<0.05) %>% 
  filter(gene %in% colnames(rna)) 
colnames(cell.gene)[1] <- "pop"
table(cell.gene$pop) %>% sort()

my_fisher_test(GWAS.lmm.LD.eqtl$group %>% unique(), colnames(rna.log), cell.gene %>% filter(pop=="B_Cells") %>% pull(gene) %>% unique())


# PRS and gene expression ####
Bcell.g <- cell.gene %>% filter(pop %in% c("Naive_B_Cells", "B_Cells")) %>% pull(gene) %>% unique()
df.DEG <- data.frame()
for (g in Bcell.g){
  
  df <- pheno
  df$gene <- rna.log[match(df$sample, rownames(rna.log)),g] 
  
  lmFit <- lm(gene~PRS_norm, df) %>% summary() %>% coefficients() %>% as.data.frame() %>% mutate(gene=g)
  
  df.DEG <- rbind(df.DEG,
                  lmFit["PRS",])
}
df.DEG$adj.pval <- p.adjust(df.DEG$`Pr(>|t|)`, method="BH")
df.DEG$fill <- ifelse(df.DEG$adj.pval<0.05, "notDEG", "DEG")

df.DEG <- df.DEG %>% arrange(adj.pval)

df.DEG$label <- NA
df.DEG$label[1:50] <- df.DEG$gene
p <- ggplot(df.DEG, aes(x=Estimate, y=-log10(adj.pval), fill=fill)) +
  geom_point(shape=21, size=2) +
  geom_text(aes(label=label), size=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_hline(yintercept=0, linetype=2) +
  geom_hline(yintercept=-log10(0.05), linetype=2) +
  scale_fill_manual(values=c(notDEG="gray77", DEG="gray35")) +
  theme_classic() + theme(legend.position = "none")
p

pdf("pdfs/PGS-memory B-genes-mean.pdf", width=3, height=3)
print(p)
dev.off()
    

# PRS and gene expression of Central_Memory_CD8 ####
Bcell.g <- cell.gene %>% filter(pop %in% c("Central_Memory_CD8")) %>% pull(gene) %>% unique()
df.DEG <- data.frame()
for (g in Bcell.g){
  
  df <- pheno
  df$gene <- rna.log[match(df$sample, rownames(rna.log)),g] 
  
  lmFit <- lm(gene~PRS_norm, df) %>% summary() %>% coefficients() %>% as.data.frame() %>% mutate(gene=g)
  
  df.DEG <- rbind(df.DEG,
                  lmFit["PRS",])
}
df.DEG$adj.pval <- p.adjust(df.DEG$`Pr(>|t|)`, method="BH")
df.DEG$fill <- ifelse(df.DEG$adj.pval<0.05, "notDEG", "DEG")

df.DEG <- df.DEG %>% arrange(adj.pval)

df.DEG$label <- NA
df.DEG$label[1:50] <- df.DEG$gene
p <- ggplot(df.DEG, aes(x=Estimate, y=-log10(adj.pval), fill=fill)) +
  geom_point(shape=21, size=2) +
  geom_text(aes(label=label), size=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_hline(yintercept=0, linetype=2) +
  geom_hline(yintercept=-log10(0.05), linetype=2) +
  scale_fill_manual(values=c(notDEG="gray77", DEG="gray35")) +
  theme_classic() + theme(legend.position = "none")
p

pdf("pdfs/PGS-memory CD8cm-genes-mean.pdf", width=3, height=3)
print(p)
dev.off()


# Examples of B cell co-regulation ####
snp_shared <- intersect(GWAS.lmm.LD$SNP, GWAS.lmm.eqtl$SNP)
k <- 0
for (snp_n in snp_shared){
  idx.snp <- which(snp$ID==snp_n)
  stopifnot(length(idx.snp)==1)
  
  idx.cqtl <- which(GWAS.lmm.LD$SNP == snp_n & GWAS.lmm.LD$family=="B_Cells")
  idx.eqtl <- which(GWAS.lmm.eqtl$SNP == snp_n)
  for (idx1 in idx.cqtl){
    pop <- GWAS.lmm.LD$group[idx1]
    for (idx2 in idx.eqtl){
      gene <- GWAS.lmm.eqtl$group[idx2]
      
      df.plot.pop <- data.frame(sample=rownames(cytof.nonnegative), variable=cytof.nonnegative[,pop], omic="CyTOF")
      df.plot.pop$rank <- 100*(df.plot.pop$variable %>% rank())/nrow(df.plot.pop)
      df.plot.pop$ind <- gsub("\\:.*", "", df.plot.pop$sample)
      df.plot.pop$geno <- geno[idx.snp, df.plot.pop$ind] %>% as.character()
      
      df.plot.gene <- data.frame(sample=rownames(rna.log), variable=rna.log[,gene], omic="RNA-seq")
      df.plot.gene$rank <- 100*(df.plot.gene$variable %>% rank())/nrow(df.plot.gene)
      df.plot.gene$ind <- gsub("\\:.*", "", df.plot.gene$sample)
      df.plot.gene$geno <- geno[idx.snp, df.plot.gene$ind] %>% as.character()
      
      df.plot.n <- rbind(df.plot.pop, df.plot.gene) %>% na.omit()
      
      print(paste0(pop, "; ", gene, "; ", snp_n, ";"))
      print(table(df.plot.pop$geno))
      print(table(df.plot.gene$geno))
      
      p <- ggplot(df.plot.n, aes(x=geno, y=rank, fill=omic)) +
        geom_quasirandom(dodge.width=0.8, shape=21, size=2) +  
        geom_boxplot() + 
        scale_fill_manual(values=c(CyTOF="#B2E6D4",`RNA-seq`="#CABAC8")) +
        xlab(snp_n) + ylab("Percentile in the cohort") +
        ggtitle(paste0(pop, "; ", gene)) +
        theme_classic() + theme(legend.position="none", plot.title = element_text(size=10))
      
      cairo_pdf(paste0("pdfs/", pop, "; ", gene, "; ", spn_n, ".pdf"), height=3, width=4)
      print(p)
      dev.off()
      
      k <- k + 1
    }
  }
}
  


# Logistic regression to sample cluster ####
logistic_fun <- function(X,Y, GWAS.summary){
  
  snp_n <- GWAS.summary$SNP
  
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
  df <- data.frame(y = Y %>% as.numeric(), cbind(X, X.geno)) 
  
  
  #logistic
  lmFit <- glm(formula(df), df, family="binomial") 
  df <- lmFit %>% summary() %>% coef()
  df <- df[-1, ]
  OR <- exp(coef(lmFit))
  CI <- exp(confint(lmFit))
  
  df.out <- data.frame(df, CI[rownames(df),], OR=OR[rownames(df)]) %>% rownames_to_column("var")
  
  df.out$group <- var.group$group[match(df.out$var, var.group$var)]
  
  return(df.out)
}


clinical <- metadata %>%
  subset(select=-c(SCAPISorIGT_id, subject, Study, Visitdate, id, visit, Birthdate, subject_id,
                   Hb, WBC, Plt, RBC, Hct, MCV, MCH, MCHC, Neut, Lymph, Mono, Eos, Baso,
                   Housing_change, Housing_current, Address_change, MaritalStatus_change, MaritalStatus_current, ShareHousehold_change, 
                   ShareHousehold_current, Employment_change, Employment_current, Tobacco_change, Tobacco_current, PerceivedHealth,       
                   Stress, PhysicalActivity, SedentaryTime_hours, SedentaryTime_minutes, SedentaryTime_unknown, TravelAbroad,          
                   TravelAbroad_country, Animals, Animals_type,
                   Calcitriol, Calcidiol, Health_Other, Common_cold_influenzae, NSAID_painmed, Bp_med, Lipid_med,             
                   Antibiotics_med, Diab_med, Med_details, Smoking, Smoking_hours)) %>%
  mutate(Gender=ifelse(Gender=="m",1,0))


#major pop
GWAS.cytof.group.LD <- GWAS.lmm.LD
for (g in names(pop.group)){
  GWAS.cytof.group.LD$group[GWAS.lmm.LD$group %in% pop.group[[g]]] <- g
}

X <- clinical[rownames(primary_cluster),] %>% scale() %>% as.data.frame()
Y <- ifelse(primary_cluster$primary_cluster==1, 1, 0)
logistic.A <- logistic_fun(X=X, 
                       Y=Y, 
                       GWAS.summary=GWAS.cytof.group.LD) 

Y <- ifelse(primary_cluster$primary_cluster==2, 1, 0)
logistic.B <- logistic_fun(X=X, 
                           Y=Y, 
                           GWAS.summary=GWAS.cytof.group.LD) 

Y <- ifelse(primary_cluster$primary_cluster==3, 1, 0)
logistic.C <- logistic_fun(X=X, 
                           Y=Y, 
                           GWAS.summary=GWAS.cytof.group.LD) 

color_palette <- c("darkseagreen",	"#A30059",	"#9e5a28ff",		
                   "tomato",
                   "#4a6fe3",	"lightblue",	"steelblue" ,
                   "bisque",	"#1CE6FF",
                   "goldenrod")
names(color_palette) <- c("Genetics", "Sex", "Age", 
                          "Blood pressure", 
                          "Body composition", "Lipid profile","Glucose homeostasis",
                          "Organ biomarker", "CRP",
                          "Lifestyle")


ggplot(logistic.A %>% filter(Pr...z..<0.05) %>% filter(X2.5..>1), aes(x=var, y=OR)) +
  geom_point() +
  geom_hline(yintercept=1, linetype=2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ggtitle("Cluster A")
ggplot(logistic.B %>% filter(Pr...z..<0.05) %>% filter(X2.5..>1), aes(x=var, y=OR)) +
  geom_point() +
  geom_hline(yintercept=1, linetype=2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ggtitle("Cluster B")
ggplot(logistic.C %>% filter(Pr...z..<0.05) %>% filter(X2.5..>1), aes(x=var, y=OR)) +
  geom_point() +
  geom_hline(yintercept=1, linetype=2) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ggtitle("Cluster C")


#prepare genotype
snp_n <- GWAS.lmm.LD$SNP

idx <- match(snp_n, snp$ID)
geno_n <- geno[idx,] %>% t() %>% as.data.frame()
colnames(geno_n) <- paste0("SNP:", snp_n)
colnames(geno_n) <- gsub(":", "\\.", colnames(geno_n) ) #formula transforms them into . anyway

subject <- gsub("\\:.*", "", rownames(primary_cluster))
geno.cytof <- geno_n[subject,] %>% as.data.frame()
geno.cytof[geno.cytof=="0/0"] <- 0
geno.cytof[geno.cytof=="0/1"] <- 1
geno.cytof[geno.cytof=="1/1"] <- 2

for (k in 1:ncol(geno.cytof)){
  geno.cytof[,k] <- as.numeric(geno.cytof[,k])
}


df.A <- data.frame()
for (n in 1:ncol(geno.cytof)){
  df <- data.frame(y=ifelse(primary_cluster$primary_cluster==1, 1, 0), x=geno.cytof[, n])
  
  #Fisher's exact test
  ind.alt <- which(df$x %in% c(1,2))
  ind.clust <- which(df$y== 1)
  
  f <- my_fisher_test(ind.alt %>% as.character(), 
                      1:nrow(df) %>% as.character(), 
                      ind.clust %>% as.character())
  
  #Logistic
  glmFit <- glm(y~x, df, family="binomial")
  df <- glmFit %>% summary() %>% coef() %>% as.data.frame()
  OR <- exp(coef(glmFit))
  CI <- exp(confint(glmFit))
  
  
  #Store
  df.A <- rbind(df.A,
                data.frame(snp_id=colnames(geno.cytof)[n],
                           family=GWAS.lmm.LD$family[n],
                           p=df["x","Pr(>|z|)"],
                           b=df["x", "Estimate"],
                           OR=OR["x"],
                           CI_low=CI["x", 1],
                           CI_high=CI["x",2],
                           p_fisher=f$p_value,
                           OR_fisher=f$odds_ratio))
}

df.B <- data.frame()
for (n in 1:ncol(geno.cytof)){
  df <- data.frame(y=ifelse(primary_cluster$primary_cluster==2, 1, 0), x=geno.cytof[, n])
  
  #Fisher's exact test
  ind.alt <- which(df$x %in% c(1,2))
  ind.clust <- which(df$y== 1)
  
  f <- my_fisher_test(ind.alt %>% as.character(), 
                      1:nrow(df) %>% as.character(), 
                      ind.clust %>% as.character())
  
  
  #Logistic
  glmFit <- glm(y~x, df, family="binomial")
  df <- glmFit %>% summary() %>% coef() %>% as.data.frame()
  OR <- exp(coef(glmFit))
  CI <- exp(confint(glmFit))
  
  
  #Store
  df.B <- rbind(df.B,
                data.frame(snp_id=colnames(geno.cytof)[n],
                           family=GWAS.lmm.LD$family[n],
                           p=df["x","Pr(>|z|)"],
                           b=df["x", "Estimate"],
                           OR=OR["x"],
                           CI_low=CI["x", 1],
                           CI_high=CI["x",2],
                           p_fisher=f$p_value,
                           OR_fisher=f$odds_ratio))
}

df.C <- data.frame()
for (n in 1:ncol(geno.cytof)){
  df <- data.frame(y=ifelse(primary_cluster$primary_cluster==3, 1, 0), x=geno.cytof[, n])
  
  #Fisher's exact test
  ind.alt <- which(df$x %in% c(1,2))
  ind.clust <- which(df$y== 1)
  
  f <- my_fisher_test(ind.alt %>% as.character(), 
                 1:nrow(df) %>% as.character(), 
                 ind.clust %>% as.character())
  
  
  #Logistic
  glmFit <- glm(y~x, df, family="binomial")
  df <- glmFit %>% summary() %>% coef() %>% as.data.frame()
  OR <- exp(coef(glmFit))
  CI <- exp(confint(glmFit))
  
  
  #Store
  df.C <- rbind(df.C,
                data.frame(snp_id=colnames(geno.cytof)[n],
                           family=GWAS.lmm.LD$family[n],
                           p=df["x","Pr(>|z|)"],
                           b=df["x", "Estimate"],
                           OR=OR["x"],
                           CI_low=CI["x", 1],
                           CI_high=CI["x",2],
                           p_fisher=f$p_value,
                           OR_fisher=f$odds_ratio))
}








##################### OLD ############################################
# Polygenic risk score per cluster ####
pop.uniq <- unique(GWAS.lmm.LD$group)
#pop.uniq <- "B_Cells"
for (pop in pop.uniq){
  idx <- which(GWAS.lmm.LD$group == pop)
  if (length(idx)>=5){
    beta <- GWAS.lmm.LD$beta[idx]
    stopifnot(!any(duplicated(idx)))
    
    idx.geno <- match(GWAS.lmm.LD$SNP[idx], snp$ID)
    stopifnot(!any(duplicated(idx.geno)))
    
    PRS <- apply(geno.X[idx.geno,], 2, function(x){sum(x*beta, na.rm=T)})
    
    pheno <- data.frame(y=cytof.nonnegative[,pop], sample=rownames(cytof.nonnegative))
    pheno$ind <- gsub("\\:.*", "", pheno$sample)
    pheno <- pheno %>% group_by(ind) %>% summarise(y=mean(y))
    pheno$PRS <- PRS[pheno$ind]
    pheno$PRS_rank <- rank(pheno$PRS)
    pheno$PRS_norm <- scale(pheno$PRS) %>% as.numeric()
    pheno$cluster <- cluster_ind$cluster[match(pheno$ind, cluster_ind$ind)]
    
    #p1: PRS vs cluster
    p1 <- ggplot(pheno %>% na.omit(), aes(x=cluster, y=PRS_rank)) + 
      geom_boxplot(aes(fill=cluster)) +
      geom_quasirandom() +
      scale_fill_manual(values=c("#e69f0099", "#56b4e999", "#009e7399")) +
      stat_compare_means(label="p.signif",
        method="wilcox.test", size=4,
        comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3"))) +
      theme_classic() + 
      xlab("") +
      theme(legend.position = "none")

    
    #p3: number variants vs cluster
    n.snp.ind <- colSums(geno.X[idx.geno,], na.rm=T)
    df.plot <- data.frame(n=n.snp.ind, ind=names(n.snp.ind))
    df.plot$cluster <- cluster_ind$cluster[match(df.plot$ind, cluster_ind$ind)]
    df.plot$rank <- rank(df.plot$n)
    
    p3 <- ggplot(df.plot %>% na.omit(), aes(x=cluster, y=rank, fill=cluster)) +
      geom_boxplot() +
      geom_quasirandom() +
      scale_fill_manual(values=c("#e69f0099", "#56b4e999", "#009e7399")) +
      stat_compare_means(label="p.signif",
        method="wilcox.test", size=4,
        comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3"))) +
      theme_classic() + theme(legend.position = "none") +
      xlab("") + ylab("Number of variants")
    
    
    #p2: PRS vs phenotype
    #fit exponential
    lm_start <- lm(log(y) ~ PRS_norm, data = pheno)
    b_start <- coef(lm_start)[2]
    a_start <- exp(coef(lm_start)[1])
    
    fit <- nls(y ~ a * exp(b * PRS_norm),
               data = pheno,
               start = list(a = a_start, b = b_start))
    
    newdf <- data.frame(PRS_norm = seq(min(df$PRS_norm),
                                       max(df$PRS_norm),
                                       length.out = 200))
    newdf$y_pred <- predict(fit, newdata = newdf)
    
    p2 <- ggplot(pheno, aes(x=PRS_norm, y=y)) +
      geom_point() +
      geom_line(data=newdf, aes(x=PRS_norm, y=y_pred), color="red", size=1) +
      theme_classic() +
      ylab("Frequency") + ggtitle(pop) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
    
    p <- ggarrange(p2,p3,p1,ncol=1)
    print(p)
    # pdf(paste0("pdfs/PRS-", pop, ".pdf"), width=3, height=5)
    # print(p)
    # dev.off()
  }
}


# Fisher's test for each SNP ####
pheno.B <- data.frame(y=cytof.group[,"B_Cells"], sample=rownames(cytof.group))
pheno.B$ind <- gsub("\\:.*", "", pheno.B$sample)
pheno.B <- pheno.B %>% group_by(ind) %>% summarise(y=mean(y))

pheno.CD8 <- data.frame(y=cytof.group[,"Central_Memory_CD8"], sample=rownames(cytof.group))
pheno.CD8$ind <- gsub("\\:.*", "", pheno.CD8$sample)
pheno.CD8 <- pheno.CD8 %>% group_by(ind) %>% summarise(y=mean(y))


Bcell_SNP <- GWAS.lmm.LD %>% filter(family=="B_Cells") %>% pull(SNP) %>% unique()
CD8_SNP <- GWAS.lmm.LD %>% filter(family=="Central_Memory_CD8") %>% pull(SNP) %>% unique()

idx.B <- which(snp$ID %in% Bcell_SNP)
idx.CD8 <- which(snp$ID %in% CD8_SNP)

geno.B <- colSums(geno.X[idx.B,], na.rm=T)
geno.CD8 <- colSums(geno.X[idx.CD8,], na.rm=T)

stopifnot(all(names(geno.B)==names(geno.CD8)))

df.geno <- data.frame(ind=names(geno.B), geno.B=geno.B, geno.CD8=geno.CD8)
df.geno$cluster <- cluster_ind$cluster[match(df.geno$ind, cluster_ind$ind)]

df.geno$B <- pheno.B$y[match(df.geno$ind, pheno.B$ind)]
df.geno$CD8 <- pheno.CD8$y[match(df.geno$ind, pheno.CD8$ind)]



ggplot(df.geno %>% na.omit(), aes(x=cluster, y=B, fill=cluster)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#e69f0099", "#56b4e999", "#009e7399")) +
  stat_compare_means(label="p.signif",
                     method="t.test", size=4,
                     comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3"))) +
  theme_classic() + theme(legend.position = "none") +
  xlab("") + ylab("Memory B cells")

ggplot(df.geno %>% na.omit(), aes(x=cluster, y=CD8, fill=cluster)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#e69f0099", "#56b4e999", "#009e7399")) +
  stat_compare_means(label="p.signif",
                     method="t.test", size=4,
                     comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3"))) +
  theme_classic() + theme(legend.position = "none") +
  xlab("") + ylab("Central Memory CD8+")




df.snp <- data.frame()
for (n in 1:nrow(GWAS.lmm.LD)){
  snp_n <- GWAS.lmm.LD$SNP[n]
  
  idx <- which(snp$ID == snp_n)
  if (length(idx)==1){
    geno_n <- data.frame(geno=geno[idx, ] %>% as.character(), ind=colnames(geno))
    
    geno_n$cluster <- cluster_ind$cluster[match(geno_n$ind, cluster_ind$ind)]
    geno_n <- na.omit(geno_n)
    geno_n$variant <- NA
    geno_n$variant[geno_n$geno == "0/0"] <- 0
    geno_n$variant[geno_n$geno == "0/1"] <- 1
    geno_n$variant[geno_n$geno == "1/1"] <- 1
    
    chi2_test <- chisq.test(table(geno_n$cluster, geno_n$variant))
    
    df.snp <- rbind(df.snp,
                    data.frame(snp=snp_n, 
                               family=GWAS.lmm.LD$family[n],
                               pval=chi2_test$p.value,
                               n.variant=sum(geno_n$variant),
                               n.variant1=sum(geno_n$variant[geno_n$cluster=="cluster1"]),
                               n.variant2=sum(geno_n$variant[geno_n$cluster=="cluster2"]),
                               n.variant3=sum(geno_n$variant[geno_n$cluster=="cluster3"])))
  }
}


# Fisher's test for cluster enrichment of SNPs ####
idx.B <- which(snp$ID %in% unique(GWAS.lmm.LD$SNP[GWAS.lmm.LD$family=="B_Cells"]))
n.snp.ind.B <- colSums(geno.X[idx.B,], na.rm=T)

df.plot <- data.frame(n=n.snp.ind.B, ind=names(n.snp.ind.B))
df.plot$cluster <- cluster_ind$cluster[match(df.plot$ind, cluster_ind$ind)]
df.plot$rank <- rank(df.plot$n)

ggplot(df.plot %>% na.omit(), aes(x=cluster, y=rank, fill=cluster)) +
  geom_boxplot() +
  geom_quasirandom() +
  scale_fill_manual(values=c("#e69f0099", "#56b4e999", "#009e7399")) +
  stat_compare_means(
    method="wilcox.test",
    comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3"))) +
  theme_classic() + theme(legend.position = "none") +
  xlab("") + ylab("Number of variants")


# primary clusters ####
gene.main <- cell.gene  %>% pull(gene) %>% unique()
da_forcluster <- rna.log[, gene.main] %>% as.data.frame()
da_forcluster$ind <- gsub("\\:.*", "", rownames(da_forcluster))

da_forcluster <- da_forcluster %>% group_by(ind) %>% summarise(across(everything(), mean, .names = "{.col}")) %>% as.data.frame()
rownames(da_forcluster) <- da_forcluster$ind
da_forcluster <- da_forcluster %>% select(-ind)
da_forcluster <- scale(da_forcluster) %>% as.data.frame()

da_forcluster <- readRDS("rna.lmm.RDS")

# KNN
kn <- 15
mat <-  t(da_forcluster)  
knn.info <- RANN::nn2( t(mat), k=kn )

# get adjacency matrix
knn <- knn.info$nn.idx
adj_527 <- matrix(0, ncol(mat), ncol(mat))
rownames(adj_527) <- colnames(adj_527) <- colnames(mat)
for(i in seq_len(ncol(mat))) {
  adj_527[i,colnames(mat)[knn[i,]]] <- 1
}

# Create graphs from adjacency matrices of kNN
g_knn_527 <- igraph::graph.adjacency(adj_527, mode="undirected")
g_knn_527 <- igraph::simplify(g_knn_527) # remove self loops

# community detection
#km_527 <- igraph::cluster_louvain(g_knn_527, resolution=.5) #kn=20 
km_527 <- igraph::cluster_fast_greedy(g_knn_527) #kn=20 

community_527 <- km_527$membership
names(community_527) <- km_527$names

cluster_ind <- data.frame(community_527)
cluster_ind$ind <- rownames(cluster_ind)
colnames(cluster_ind)[1] <- "primary_cluster"

cluster_ind$cluster <- paste0("cluster", cluster_ind$primary_cluster)

table(cluster_ind$cluster)


# Sanity checks ####
pheno.df <- data.frame(B=cytof.group[,"B_Cells"], 
                       CD8=cytof.group[,"Central_Memory_CD8"], 
                       sample=rownames(cytof.group))
pheno.df$cluster <- primary_cluster$cluster[match(pheno.df$sample, primary_cluster$sample)]
pheno.df <- na.omit(pheno.df)

ggplot(pheno.df, aes(x=cluster, y=B, fill=cluster)) +
  geom_quasirandom(color="black", shape=21) +
  geom_boxplot(color="black", fill=NA, outlier.shape=NA) +
  scale_fill_manual(values=c(cluster1="#e69f0099", cluster2="#56b4e999", cluster3="#009e7399")) +
  stat_compare_means(label="p.signif",
    method="wilcox.test",
    comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3"))) +
  theme_classic() + theme(legend.position = "none") +
  xlab("") + ylab("Memory B cells")




most_frequent <- function(x){
  tab <- table(x)
  if (any(tab>3)){
    freq <- names(tab)[tab>3]
  } else {
    freq <- NA
  }
  
  return(freq)
}


pheno.df$ind <- gsub("\\:.*", "", pheno.df$sample)
pheno.df.ind <- pheno.df %>% group_by(ind) %>% summarise(B_mean=mean(B), CD8_mean=mean(CD8))
pheno.df.ind$cluster <- cluster_ind$cluster[match(pheno.df.ind$ind, cluster_ind$ind)]

ggplot(pheno.df.ind %>% na.omit(), aes(x=cluster, y=B_mean, fill=cluster)) +
  geom_quasirandom(color="black", shape=21) +
  geom_boxplot(color="black", fill=NA, outlier.shape=NA) +
  scale_fill_manual(values=c(cluster1="#e69f0099", cluster2="#56b4e999", cluster3="#009e7399")) +
  stat_compare_means(
                     method="wilcox.test",
                     comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3"))) +
  theme_classic() + theme(legend.position = "none") +
  xlab("") + ylab("Memory B cells")

ggplot(pheno.df.ind %>% na.omit(), aes(x=cluster, y=CD8_mean, fill=cluster)) +
  geom_quasirandom(color="black", shape=21) +
  geom_boxplot(color="black", fill=NA, outlier.shape=NA) +
  scale_fill_manual(values=c(cluster1="#e69f0099", cluster2="#56b4e999", cluster3="#009e7399")) +
  stat_compare_means(
                     method="wilcox.test",
                     comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3"))) +
  theme_classic() + theme(legend.position = "none") +
  xlab("") + ylab("Central_Memory_CD8")




# Clinical factors and clusters ####
clinical.variables <- c("Gender", "Age_at_Visit", "SBP", "DBP", "ProBNP", "TNT", "CRP", "ALAT", "GGT", "HDL", "LDL", "TG")

df.plot <- cbind(primary_cluster, clinical[rownames(primary_cluster),] %>% apply(2, rank)) %>% 
  reshape2::melt(id.vars=c("sample", "primary_cluster", "ind", "visit", "cluster")) %>% filter(variable %in% clinical.variables)

p <- ggplot(df.plot, 
       aes(x=cluster, y=value)) +
  facet_wrap(~variable) +
  geom_boxplot(aes(fill=cluster)) +
  geom_quasirandom(alpha=0.6) +
  scale_fill_manual(values=c("#e69f0099", "#56b4e999", "#009e7399")) +
  stat_compare_means(label="p.signif",
    method="wilcox.test",
    comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3"))) +
  theme_classic() +
  theme(strip.background = element_blank())
p

pdf("boxplot-clusters.pdf", height=9, width=9)
print(p)
dev.off()

df.lm <- data.frame()
for (var in clinical.variables){
  df <- df.plot %>% filter(variable == var)
  df$sex <- metadata$Gender[match(df$sample, metadata$id)]
  df$age <- metadata$Age_at_Visit[match(df$sample, metadata$id)]
  
  lmFit <- lm(value ~ cluster, df)
  
  emm <- emmeans(lmFit, ~ cluster)
  contr <- pairs(emm)
  df.lm <- rbind(df.lm, 
                 summary(contr) %>% as.data.frame() %>% mutate(var=var))
}


# Plot example SNP-frequency associations ####
df.plot <- data.frame(snp=GWAS.lmm.LD$SNP, pop=GWAS.lmm.LD$group) 

for (n in 1:nrow(df.plot)){
  idx <- match(df.plot$snp[n], snp$ID)
  
  df.plot.pop <- data.frame(sample=rownames(cytof), variable=cytof.nonnegative[,df.plot$pop[n]], omic="CyTOF")
  df.plot.pop$ind <- gsub("\\:.*", "", df.plot.pop$sample)
  df.plot.pop$geno <- geno[idx, df.plot.pop$ind] %>% as.character()

  df.plot.pop <- df.plot.pop %>% group_by(ind,omic,geno) %>% summarise(variable_mean=mean(variable))
  df.plot.pop$rank <- 100*(df.plot.pop$variable_mean %>% rank())/nrow(df.plot.pop)
  
  # df.plot.gene <- data.frame(sample=rownames(rna.log), variable=rna.log[,df.plot$gene[n]], omic="RNA-seq")
  # df.plot.gene$rank <- 100*(df.plot.gene$variable %>% rank())/nrow(df.plot.gene)
  # df.plot.gene$ind <- gsub("\\:.*", "", df.plot.gene$sample)
  # df.plot.gene$geno <- geno[idx, df.plot.gene$ind] %>% as.character()
  # 
  # df.plot.n <- rbind(df.plot.pop, df.plot.gene) %>% na.omit()
  
  # print(paste0(df.plot$pop[n], "; ", df.plot$gene[n],"; ", df.plot$snp[n], ";"))
  # print(table(df.plot.pop$geno))
  # print(table(df.plot.gene$geno))
   
  p <- ggplot(df.plot.pop %>% na.omit(), aes(x=geno, y=rank)) + 
    geom_boxplot() + 
    geom_quasirandom(alpha=0.6) + 
    #geom_quasirandom(dodge.width=0.8,  alpha=0.6) + 
    theme_classic() + xlab("Genotype") + ylab("Percentile in the cohort") +
    ggtitle(paste0(df.plot$pop[n], "; ", df.plot$gene[n],"; ", df.plot$snp[n])) +
    theme(plot.title = element_text(size = 20))
  
  
  cairo_pdf(paste0("pdfs/", df.plot$pop[n], "-", n, ".pdf"), height=3, width=4)
  print(p)
  dev.off()
}









######################## OLD ###########################################

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


var.group.SNP <- data.frame(var=unique(c(GWAS.lmm$SNP)),
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
X <- clinical[rownames(cytof),] %>% scale() %>% as.data.frame()
Y <- cytof %>% scale() %>% as.data.frame()
lm.cytof.all <- lm_fun(X=X, 
                       Y=Y, 
                       GWAS.summary=GWAS.lmm.LD)
#saveRDS(lm.cytof.all, "lm.cytof.all.RDS")

out.varExpl <- plot_varExpl(lm.cytof.all)

df.plot <- data.frame(out.varExpl$df)
df.plot <- merge(df.plot, macro.anno.my, by.x="x", by.y="row.names")
colnames(df.plot)[ncol(df.plot)] <- "family"

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

# Extract all associations ####
df.lm.all.cytof <- data.frame()
for (n in 1:length(lm.cytof.all)){
  print(paste0(n, "/", length(lm.cytof.all)))
  df.lm.all.cytof <- rbind(df.lm.all.cytof,
                           lm.cytof.all[[n]] %>% mutate(adj.pval=p.adjust(`Pr(>|t|)`, method="BH")) %>% mutate(gene=names(lm.cytof.all)[n]) )
}
df.lm.all.cytof %>% filter(adj.pval<0.05) %>% pull(var) %>% table() %>% sort(decreasing=T) %>% head(n=10)
df.lm.all.cytof %>% filter(adj.pval<0.05) %>% pull(group) %>% table() %>% sort(decreasing=T) %>% head(n=10)


# Variance explained for major immune frequencies ####
GWAS.cytof.group.LD <- GWAS.lmm.LD
for (g in names(pop.group)){
  GWAS.cytof.group.LD$group[GWAS.lmm.LD$group %in% pop.group[[g]]] <- g
}

X <- clinical[rownames(cytof.group),] %>% scale() %>% as.data.frame()
Y <- cytof.group %>% scale() %>% as.data.frame()
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

#####################


# Manhattan of GWAS ####
MAX_SNP <- 5e4

#Prepare GWAS CyTOF
df_plot_cytof <- GWAS.lmm[order(GWAS.lmm$P, decreasing=F), c("CHR", "BP", "P", "group", "family", "SNP")]
df_plot_cytof <- df_plot_cytof[!df_plot_cytof$CHR %in% c("PAR1", "PAR2"),]
df_plot_cytof$CHR <- as.numeric(df_plot_cytof$CHR)

df_plot_cytof$pop <- df_plot_cytof$group
df_plot_cytof <- df_plot_cytof[, c("CHR", "BP", "P", "pop", "family", "group", "SNP")]
N_cqtl <- nrow(df_plot_cytof)

#Prepare GWAS RNA-seq
df_plot_rna <- GWAS.lmm.eqtl[order(GWAS.lmm.eqtl$P, decreasing=F)[1:MAX_SNP], c("CHR", "BP", "P", "group", "SNP")]
df_plot_rna <- df_plot_rna[!df_plot_rna$CHR %in% c("PAR1", "PAR2"),]
df_plot_rna$CHR <- as.numeric(df_plot_rna$CHR)

df_plot_rna$pop <- ""
df_plot_rna$family <- ""
df_plot_rna <- df_plot_rna[, c("CHR", "BP", "P", "pop", "family", "group", "SNP")]
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
df_plot_cytof <- df_plot[idx.cytof, c("CHR", "BP", "P", "pop", "family", "group", "BP_cum", "color")]
df_plot_rna <- df_plot[idx.rna, c("CHR", "BP", "P", "group", "SNP", "BP_cum", "color")]


# Plot CyTOF, RNA-seq and cytokine GWAS ####
pop.palette <- c("Monocytes_classical"="sienna2", "TEMRA_CD8"="#E31A1C", "NK_Cells"="turquoise", 
                 "Naive_CD8"="#FB9A99", "Naive_CD4"="#B2DF8A",             
                 "Basophils"="goldenrod2", "Monocytes_nonclassical"="sienna2", "Central_Memory_CD8"="#E31A1C", 
                 "Monocytes_intermediate"="sienna2", "B_Cells"="dodgerblue4",               
                 "Myeloid_DC"="lemonchiffon1", "Plasmacytoid_DC"="lemonchiffon1", "Naive_B_Cells"="dodgerblue2", 
                 "Central_Memory_CD4"="#33A02C", "Effector_Memory_CD4"="#33A02C",   
                 "TEMRA_CD4"="#33A02C", "Effector_Memory_CD8"="#E31A1C", "ILCs"="gray78")
#df_plot_cytof$symbol <- ifelse(df_plot_cytof$P<pval.cytof, df_plot_cytof$group, "")

df_plot_cytof$color <- ifelse(df_plot_cytof$P<pval.cytof, df_plot_cytof$family, df_plot_cytof$color)
p1 <- ggplot(df_plot_cytof, aes(x=BP_cum, y=-log10(P), fill=color)) +
  geom_point() + 
  geom_hline(yintercept=-log10(pval.cytof), linetype=2) +
  theme_bw() +
  #geom_text_repel(aes(label=symbol), cex=2, color="black") +
  scale_x_continuous(label = chr.df$CHR, breaks = chr.df$BP_mean) +
  scale_fill_manual(values=c(pop.palette, even="gray77", odd="gray23")) + 
  scale_size_continuous(range = c(0.5,3)) + 
  labs(x = NULL, y = "-log<sub>10</sub>(p)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0, size=8, vjust=0.5),
        legend.position="none",
        text=element_text(size=12)) +
  ylab("-Log10 P-value")
#p1

df_plot_rna <- df_plot_rna %>% arrange(P)

#Add label for top genes
N_GENES_LABEL <- 40
GENES_LABEL <- character(0)
n <- 1
counter <- 0
df_plot_rna$label <- ""
while (counter<N_GENES_LABEL){
  g <- df_plot_rna$group[n]
  if (!g %in% GENES_LABEL){
    GENES_LABEL <- c(GENES_LABEL, g)
    df_plot_rna$label[n] <- g
    counter <- counter + 1
  }
  n <- n + 1
}

p2 <- ggplot(df_plot_rna, aes(x=BP_cum, y=-log10(P), fill=color)) +
  geom_point() +
  geom_hline(yintercept=-log10(pval.rna), linetype=2) +
  scale_fill_manual(values=c(even="gray77", odd="gray23")) + 
  geom_text(aes(label=label), cex=2, color="black") +
  theme_bw() +
  scale_x_continuous(label = chr.df$CHR, breaks = chr.df$BP_mean) +
  scale_size_continuous(range = c(0.5,3)) + 
  labs(x = NULL, y = "-log<sub>10</sub>(p)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0, size=8, vjust=0.5),
        legend.position="none",
        text=element_text(size=12)) +
  ylab("-Log10 P-value")
#p2

p <- ggarrange(p1,p2, ncol=1)
p
cairo_pdf("pdfs/GWAS-cytof-rna-221025.pdf", width=5, height=6.6)
print(p)
dev.off()


#################################