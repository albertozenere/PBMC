
pack_R <- c("dplyr","ggplot2","ggrepel","umap", "edgeR",
            "RColorBrewer", "pheatmap", "tidyverse",
            "igraph", "ForceAtlas2", "biomaRt", "extrafont",
            "org.Hs.eg.db", "ggpubr",
            "clusterProfiler", "RANN", "dbscan", "cluster",
            "corrplot", "igraph",
            "ggbeeswarm", "glmnet")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(1)


# Load ####

# Annotation of individuals
anno <- read.csv("C:/Users/albze08/Desktop/postDoc/genome-protein/data/WELLNESS/Wellness/Data/Wellness_barcodes.txt", sep="\t", header=T)
anno <- anno[anno$Sample.type=="Helblod",]

# metadata
metadata <- read.csv("C:/Users/albze08/Desktop/postDoc/genome-protein/data/WELLNESS/Wellness/Data/Metadata/complete.clinical.data.wellness.t2d.txt", sep="\t", header=T)
metadata <- metadata[- which(metadata$Study=="T2D"), ]
metadata$visit <- as.character(metadata$visit)
metadata$subject_id <- gsub("1-", "", metadata$subject_id)
metadata$id <- paste0(metadata$subject_id, ":", metadata$visit)
rownames(metadata) <- metadata$id

# Filter metadata
clinical <- metadata
rownames(clinical) <- metadata$id
clinical$Gender <- ifelse(clinical$Gender=="f", 0, 1)
clinical <- subset(clinical, select=-c(visit, Number, subject_id, Study, Visitdate, subject_short, id, 
                                       WBC, Neut, Lymph, Mono, Eos, Baso,
                                       Hct, Hb, MCH, MCHC, RBC, MCV, Plt)) %>%
  na.omit()

#RNA-seq S3WP
rna_s3wp <- read.table("data/wellness_PBMC_v16_norm.txt", sep="\t", header = T)

#Cytof S3WP
cytof <- read.table("data/original.cytof.txt", sep="\t", header = T)
rownames(cytof) <- cytof$SampleID
cytof <- subset(cytof, select=-SampleID)
rownames(cytof) <- gsub("_", ":", rownames(cytof))


#HPA elevated genes
Bmemory <- read.table("HPA/memory b cell elevated.tsv", sep="\t", header=T)
intermediate <- read.table("HPA/intermediate monocyte elevated.tsv", sep="\t", header=T)
classical <- read.table("HPA/classical monocyte elevated.tsv", sep="\t", header=T)
myeloid <- read.table("HPA/myeloid elevated.tsv", sep="\t", header=T)
Bnaive <- read.table("HPA/naive b cell elevated.tsv", sep="\t", header=T)
plasmacytoid <- read.table("HPA/plasmacytoid elevated.tsv", sep="\t", header=T)
nonclassical <- read.table("HPA/nonclassical monocyte elevated.tsv", sep="\t", header=T)
basophil <- read.table("HPA/basophil elevated.tsv", sep="\t", header=T)
cd8_naive <- read.table("HPA/naive cd8 elevated.tsv", sep="\t", header=T)
cd4_naive <- read.table("HPA/naive cd4 elevated.tsv", sep="\t", header=T)
cd8_memory <- read.table("HPA/memory cd8 elevated.tsv", sep="\t", header=T)
cd4_memory <- read.table("HPA/memory cd4 elevated.tsv", sep="\t", header=T)

monocytes <- rbind(intermediate,classical,nonclassical) %>% distinct()
DC <- rbind(plasmacytoid, myeloid) %>% distinct()

HPA.genes <- c(Bmemory$Ensembl, monocytes$Ensembl,DC$Ensembl,Bnaive$Ensembl,
               cd8_naive$Ensembl,cd4_naive$Ensembl,cd8_memory$Ensembl,cd4_memory$Ensembl) %>% unique()

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


# ICC of frequencies, expression, cytokines ####
library("lme4")
library("lmerTest")

icc_fun <- function(df){  
  lme_model <- lmer(y ~ 1 + sex + age + (1 | ind), data = df)
  lme.var <- VarCorr(lme_model) %>% as.data.frame() %>% column_to_rownames("grp")
  s.btw <- lme.var["ind", "vcov"]
  s.wtn <- lme.var["Residual", "vcov"]
  icc <- s.btw/(s.btw+s.wtn)
  return(data.frame(icc=icc, 
                    sex_b=summary(lme_model)$coefficient["sex", "Estimate"],
                    age_b=summary(lme_model)$coefficient["age", "Estimate"], 
                    sex_p=summary(lme_model)$coefficient["sex", "Pr(>|t|)"],
                    age_p=summary(lme_model)$coefficient["age", "Pr(>|t|)"]))
}

icc.df <- data.frame()
print.prc <- 0
for (n in 1:ncol(rna.log)){
  if (n/ncol(rna.log)>=print.prc){
    print(paste0(100*print.prc, "% done"))
    print.prc <- print.prc + 0.1
  }
  df <- data.frame(y=rna.log[,n] %>% scale(), 
                   ind= gsub("\\:.*", "", rownames(rna.log)) %>% as.character())
  df$sex <- ifelse(metadata$Gender[match(rownames(rna.log), rownames(metadata))]=="m",1,0)
  df$age <- metadata$Age_at_Visit[match(rownames(rna.log), rownames(metadata))] %>% scale()
  icc <- icc_fun(df)
  icc.df <- rbind(icc.df, data.frame(icc=icc, var=colnames(rna.log)[n], omic="RNA"))
}

for (n in 1:ncol(cytof)){
  df <- data.frame(y=cytof[,n] %>% scale(), ind= gsub("\\:.*", "", rownames(cytof)) %>% as.character())
  df$sex <- ifelse(metadata$Gender[match(rownames(cytof), rownames(metadata))]=="m",1,0)
  df$age <- metadata$Age_at_Visit[match(rownames(cytof), rownames(metadata))] %>% scale()
  icc <- icc_fun(df)
  icc.df <- rbind(icc.df, data.frame(icc=icc, var=colnames(cytof)[n], omic="CyTOF"))
}

for (n in 1:ncol(protein)){
  df <- data.frame(y=protein[,n] %>% scale(), ind= gsub("\\:.*", "", rownames(protein)) %>% as.character())
  df$sex <- ifelse(metadata$Gender[match(rownames(protein), rownames(metadata))]=="m",1,0)
  df$age <- metadata$Age_at_Visit[match(rownames(protein), rownames(metadata))] %>% scale()
  icc <- icc_fun(df)
  icc.df <- rbind(icc.df, data.frame(icc=icc, var=colnames(protein)[n], omic="Protein"))
}


#saveRDS(icc.df, "icc.df.RDS")
icc.df <- readRDS("icc.df.RDS")

top.rna <- icc.df %>% filter(omic=="RNA") %>% arrange(desc(icc.icc)) %>% head(n=50) %>% pull(var)
top.imm <- icc.df %>% filter(omic=="CyTOF") %>% arrange(desc(icc.icc)) %>% head(n=10) %>% pull(var)
top.prot <- icc.df %>% filter(omic=="Protein") %>% arrange(desc(icc.icc)) %>% head(n=50) %>% pull(var)

icc.df$label <- ifelse(icc.df$var %in% c(top.rna, top.prot), icc.df$var, "")
label.df <- icc.df %>% filter(label!="")

icc.df$omic <- factor(icc.df$omic, levels=c("Protein", "RNA", "CyTOF") %>% rev())
p1 <- ggplot(icc.df, aes(x=omic, y=icc.icc)) +
  geom_violin() +
  geom_boxplot(width=0.3, outlier.shape=NA) +
  theme_classic() +
  geom_point(data=label.df, aes(x=omic, y=icc.icc)) +
  geom_text_repel(data=label.df, aes(label=label), 
                  size=2, 
                  max.overlaps=Inf, 
                  force=2, 
                  box.padding=0.2, 
                  point.padding=0.1) +
  xlab("") + ylab("ICC on normalized data")
p1
pdf("ICC-omics-1012.pdf", width=6, height=4)
print(p1)
dev.off()

# age
icc.df$icc.age_p.BH[icc.df$omic=="RNA"] <- p.adjust(icc.df$icc.age_p[icc.df$omic=="RNA"], method="BH")
p2 <- ggplot(icc.df %>% filter(omic=="RNA"), aes(x=icc.age_b, y=-log10(icc.age_p.BH))) + 
  geom_point(color="gray") +
  geom_hline(yintercept=-log10(0.05), linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Estimate") + ylab("-Log10 (adj-pval)") + ggtitle("Age")
p2

#sex
df.plot.rna <- icc.df %>% filter(omic=="RNA")
gene.sex <- (df.plot.rna %>% arrange(icc.sex_p) %>% pull(var) %>% head(20))
df.plot.rna$label <- ifelse(df.plot.rna$var %in% gene.sex, df.plot.rna$var, "")
df.plot.rna$icc.sex_p.BH <- p.adjust(df.plot.rna$icc.sex_p, method="BH")
df.plot.rna$group <- ifelse(df.plot.rna$icc.sex_p.BH>0.05, "not", ifelse(df.plot.rna$icc.sex_b>0, "up", "down"))
p3 <- ggplot(df.plot.rna, aes(x=icc.sex_b, y=-log10(icc.sex_p.BH), color=group)) + 
  geom_point() +
  geom_hline(yintercept=-log10(0.05), linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_text_repel(data=df.plot.rna, aes(label=label), 
                  size=2, 
                  max.overlaps=Inf, 
                  force=2, 
                  box.padding=0.2, 
                  point.padding=0.1) +
  scale_color_manual(values=c(up="skyblue", down="tomato", not="gray")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
  xlab("Estimate") + ylab("-Log10 (adj-pval)") + ggtitle("Sex")
p3

p <- ggarrange(p1,p2,p3, nrow=1, widths=c(4,3,3))

pdf("ICC.pdf", height=4, width=10)
print(p)
dev.off()

#some stats
df.plot.rna <- icc.df %>% filter(omic=="RNA")

p.sex <- df.plot.rna$icc.sex_p
p.sex.adj <- p.adjust(p.sex, method="BH")
table(p.sex.adj<0.05)

df.plot.rna$var[p.sex.adj<0.05 & df.plot.rna$icc.icc>0.8]

p.age <- df.plot.rna$icc.age_p
p.age.adj <- p.adjust(p.age, method="BH")
table(p.age.adj<0.05)



# Frequencies ICC ####

#Major pop
icc.df.major <- data.frame()
for (n in 1:ncol(cytof.group)){
  df <- data.frame(y=cytof.group[,n] %>% scale(), ind= gsub("\\:.*", "", rownames(cytof.group)) %>% as.character())
  df$sex <- ifelse(metadata$Gender[match(rownames(cytof.group), rownames(metadata))]=="m",1,0)
  df$age <- metadata$Age_at_Visit[match(rownames(cytof.group), rownames(metadata))] %>% scale()
  icc <- icc_fun(df)
  icc.df.major <- rbind(icc.df.major, data.frame(icc=icc, family=colnames(cytof.group)[n]))
}
family.ord <- icc.df.major %>% arrange(desc(icc.icc)) %>% pull(family)
icc.df.major$family <- factor(icc.df.major$family, levels=family.ord)

#Specific pop
df.cytof <- icc.df %>% filter(omic=="CyTOF") 
df.cytof$family <- macro.anno.my[df.cytof$var]
df.cytof$family <- factor(df.cytof$family, levels=family.ord)

p <- ggplot() + 
  geom_point(data=df.cytof, aes(x=family, y=icc.icc, color=family), shape=8, position = position_nudge(x = 0.2)) +
  geom_point(data=icc.df.major, aes(x=family, y=icc.icc, fill=family), shape=21, size=3, position = position_nudge(x = -0.2)) +
  geom_hline(yintercept=median(df.cytof$icc.icc), linetype=2) +
  scale_color_manual(values=pop.palette) +
  scale_fill_manual(values=pop.palette) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, vjust=0, hjust=0.1), legend.position="none") +
  xlab("") + ylab("ICC on normalized data") +
  ylim(0,1)
p

pdf("ICC-CyTOF.pdf", width=6, height=4)
print(p)
dev.off()

df.cytof.family <- df.cytof %>% group_by(family) %>% 
  summarise(median_icc=median(icc.icc), n=n()) %>% as.data.frame()
df.cytof$family <- factor(df.cytof$family, levels=df.cytof.family$family[order(df.cytof.family$median_icc, decreasing=T)])

#Plot

pop.palette <- c(
  "Monocytes_classical" = "#C24100", "Monocytes_intermediate"  = "#E66101", "Monocytes_nonclassical"  = "#FDB863",  
  "Basophils" = "#FFD92F", "Myeloid_DC" = "#E6AB02", "Plasmacytoid_DC" = "#B8A200",   
  "B_Cells" = "dodgerblue4", "Naive_B_Cells" = "dodgerblue2",  
  "Naive_CD4" = "#33A02C", "Central_Memory_CD4" = "#66C2A5","Effector_Memory_CD4" = "#1B9E77", "TEMRA_CD4" = "#41AE76",  
  "Naive_CD8" = "#FB9A99", "Central_Memory_CD8" = "#E31A1C", "Effector_Memory_CD8" = "#B2182B", "TEMRA_CD8" = "#67001F",  
  "NK_Cells" = "#CC78BC",  
  "ILCs" = "#B3B3B3"   
)

p1 <- ggplot(df.cytof, aes(x=family, y=icc.icc, color=family)) +
  geom_boxplot(outlier.shape=NA) +
  scale_color_manual(values=pop.palette) +
  geom_text(data=df.cytof.family, aes(label=paste0("n=", n), y=1), color="black", size=3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, vjust=0, hjust=0.1), legend.position="none") +
  xlab("") + ylab("ICC on normalized data") +
  ylim(0,1)

df.cytof$var <- factor(df.cytof$var, levels=df.cytof$var[order(df.cytof$icc.icc, decreasing=T)])
p1 <- ggplot(df.cytof, aes(x=var, y=icc.icc, color=family)) +
  geom_point(size=2) +
  geom_hline(yintercept=median(icc.df$icc.icc[icc.df$omic=="CyTOF"]), linetype=2) +
  geom_text_repel(aes(label=label), size=3) +
  scale_color_manual(values=pop.palette) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, vjust=0, hjust=0.1), legend.position="none") +
  xlab("") + ylab("ICC on normalized data") +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    legend.position = "none"
  ) +
  ylim(0,1)
p1

p12 <- ggplot(df.cytof, aes(x=var, y=family)) +
  geom_tile() +
  theme_classic() +
  xlab("") + ylab("") +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    legend.position = "none"
  )

ggarrange(p1,p12, nrow=2)

pdf("ICC-CyTOF.pdf", width=8, height=4)
print(p1)
dev.off()


# ICC of major immune frequencies ####
icc.df.major <- data.frame()
for (n in 1:ncol(cytof.group)){
  df <- data.frame(y=cytof.group[,n] %>% scale(), ind= gsub("\\:.*", "", rownames(cytof.group)) %>% as.character())
  df$sex <- ifelse(metadata$Gender[match(rownames(cytof.group), rownames(metadata))]=="m",1,0)
  df$age <- metadata$Age_at_Visit[match(rownames(cytof.group), rownames(metadata))] %>% scale()
  icc <- icc_fun(df)
  icc.df.major <- rbind(icc.df.major, data.frame(icc=icc, var=colnames(cytof.group)[n]))
}

pop.palette <- c(
  "Monocytes_classical" = "#C24100", "Monocytes_intermediate"  = "#E66101", "Monocytes_nonclassical"  = "#FDB863",  
  "Basophils" = "#FFD92F", "Myeloid_DC" = "#E6AB02", "Plasmacytoid_DC" = "#B8A200",   
  "B_Cells" = "dodgerblue4", "Naive_B_Cells" = "dodgerblue2",  
  "Naive_CD4" = "#33A02C", "Central_Memory_CD4" = "#66C2A5","Effector_Memory_CD4" = "#1B9E77", "TEMRA_CD4" = "#41AE76",  
  "Naive_CD8" = "#FB9A99", "Central_Memory_CD8" = "#E31A1C", "Effector_Memory_CD8" = "#B2182B", "TEMRA_CD8" = "#67001F",  
  "NK_Cells" = "#CC78BC",  
  "ILCs" = "#B3B3B3"   
)

df.cytof$family <- factor(df.cytof$family, levels=df.cytof.family$family[order(df.cytof.family$median_icc, decreasing=T)])
icc.df.major$var <- factor(icc.df.major$var, levels=icc.df.major$var[order(icc.df.major$icc.icc, decreasing = T)])
p <- ggplot(icc.df.major, aes(x=var, y=icc.icc, color=var)) +
  geom_point(size=3) +
  geom_hline(yintercept=median(icc.df.major$icc.icc), linetype=2) +
  scale_color_manual(values=pop.palette) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, vjust=0, hjust=0.1), legend.position="none") +
  xlab("") + ylab("ICC on normalized data") +
  ylim(0,1)
p

pdf("ICC-CyTOF-major.pdf", width=5, height=4)
print(p)
dev.off()

# Examples
gene.sex.all <- icc.df %>% filter(omic=="RNA" & icc.sex_p.BH<0.05) %>% pull(var)
g.vec <- icc.df %>% filter(var %in% gene.sex.all) %>% arrange(desc(icc.icc)) %>% head(n=20) %>% pull(var)
for (g in g.vec){
  df.plot <- data.frame(g=rna.log[sample.cytof.rna,g] %>% scale(), sex=metadata$Gender[match(sample.cytof.rna, metadata$id)], 
                        ind=metadata$subject_id[match(sample.cytof.rna, metadata$id)],
                        visit=gsub(".*\\:", "", sample.cytof.rna) %>% as.numeric())
 
  p <- ggplot(df.plot, aes(x=visit, y=g, color=sex)) +
    geom_point() +
    geom_line(aes(group=ind)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
    scale_color_manual(values=c(f="tomato", m="skyblue")) +
    ggtitle(g) + ylab("")
   print(p)
   pdf(paste0("sex-examples/", g, ".pdf"), width=4, height=3)
   print(p)
   dev.off()
}

#GSEA on ICC ####
g <- icc.df %>% filter(omic=="RNA") %>% pull(var)

gene_entrez <- bitr(g, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_entrez$icc <- icc.df$icc.icc[match(gene_entrez$SYMBOL, icc.df$var)]
gene_entrez <- arrange(gene_entrez, desc(icc))

g.ranked <- setNames(gene_entrez$icc, gene_entrez$ENTREZID)

SnowParam(1)
gsea_result <- gseGO(
  geneList = g.ranked,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  scoreType = "pos",
  verbose = TRUE
)

pdf("GSEA-ICC.pdf", width=6, height=5)
dotplot(gsea_result, color="NES")
dev.off()

# Genes: Receptors, downstream, end-products ####
immune_genes <- list(
  Receptors=c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", # HLA (MHC)
              "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", # TLRs
              "RIGI", "IFIH1", "DHX58", # RIG-I-like receptors
              "CLEC7A", "CLEC4E", "CLEC9A", # CLRs
              "NOD1", "NOD2", # NOD-like receptors
              "IL2RA", "IL7R", "IL12RB1", "IL12RB2", "IL4R", "IL6R", # Interleukin receptors
              "CCR1", "CCR2", "CCR5", "CCR7", "CXCR3", "CXCR4",
              "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2",
              "KIR2DL1", "KIR2DL2", "KIR3DL1", "KIR3DL2", "KIR2DS1", "KIR2DS4", # KIR family
              "KLRC1", "KLRK1", # NKG2 family
              "TNFRSF1A", "TNFRSF1B", "TNFRSF4", "TNFRSF9",
              "FCGR1A", "FCGR2A", "FCGR3A", "FCER1A", "FCER1G",
              "TRA", "TRB", "TRG", "TRD", # TCR chains
              "CD3D", "CD3E", "CD3G", "CD247", # CD3 complex
              "IGLL1", "CD79A", "CD79B", # BCR complex
              "CD19", "CR2", "CD81", "CD21", "CD28", "CTLA4", "PDCD1", "ICOS", "CD40", "CD80", "CD86", "PDCD1LG2" # Co-receptors
  ),
  Downstream=c("IFIT1", "IFIT3", "MX1", "OAS1", "ISG15", # Representative ISGs
               "TRAF2", "TRAF6", #TNF response
               "CD28", "CTLA4", "PDCD1", "ICOS", "CD40", "CD80", "CD86",
               "MYD88", "TICAM1", "MAVS", "TYROBP",
               "NFKB1", "RELA", "RELB",
               "IRF3", "IRF7", "IRF9",
               "STAT1", "STAT3", "STAT4", "STAT5A",
               "CD4", "CD8A", "CD8B", # Co-receptors
               "LCK", "ZAP70", "LAT", "LCP2", "ITK", "PLCG1", "GRB2", "VAV1", "NFATC1", # Downstream signaling
               "CSK", "CBL", "PTPN6",
               "SYK", "BTK", "BLNK", "PLCG2", "PIK3CD", "VAV1", "GRB2", "AKT1", # Downstream signaling
               "NFKB1", "RELA", "NFATC1", "IRF4", "PRDM1", # Transcription factors
               "LCK", "ZAP70", "LAT", "LCP2", "ITK", "PLCG1", "GRB2", "VAV1", "NFATC1", "CD28", "CSK", # T cell signaling
               "SYK", "BTK", "BLNK", "PLCG2", "PIK3CD", "VAV1", "GRB2", "AKT1", "CD19", "FCGR2B", # B cell signaling
               "PTPN6", "FCGR2B" # Negative regulators# Regulatory proteins
  ),
  `End Product` = c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", 
                    "IL32", "IL23A", "IL1B", "IL18", "IL15", "IL16", "TGFB1", "TGFBI", "CCL5", "CCL4", "CCL3",
                    "TNFAIP1", "TNFAIP2", "TNFAIP3", "TNFAIP8", "TNFAIP8L1", "TNFAIP8L2",
                    "LYZ", "S100A8", "S100A9",
                    "FTH1", "FTL", "SERPINA1", "FTH1")
)
# Define the categories
categories <- list(
  `HLA class I` = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
                    "HHLA1", "HHLA2", "HHLA3"),
  `HLA class II` = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5",
                     "HLA-DQA1", "HLA-DQA1.1", "HLA-DQB1", "HLA-DQB2",
                     "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB",
                     "HLA-DOA", "HLA-DOB"),
  `TLRs` = grep("TLR", colnames(rna), value=T),
  `RIG-I-like receptors` = c("RIGI", "IFIH1", "DHX58", "LGP2", "DDX58"),
  `CLRs` = grep("CLEC", colnames(rna), value=T),
  cGAS_STING = c("MB21D1", "TMEM173", "TBK1", "IRF3", "IFNB1"),
  `NOD-like receptors` = c("NOD1", "NOD2", grep("NLRP|NLRC|NLRX", colnames(rna), value=T)),
  `Interleukin receptors` = grep("IL[0-9]+R", colnames(rna), value = TRUE),
  `Chemokine receptors` = c("CCR1", "CCR2", "CCR5", "CCR7", "CXCR3", "CXCR4", "CX3CR1", "CCR4"),
  `IFN receptors` = c("IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IFNLR1"),
  `KIR family` = grep("KIR[0-9]+D", colnames(rna), value = TRUE),
  `NKG2 family` = c("KLRC1", "KLRK1", "KLRC2", "KLRC3", "KLRK1"),
  `TNF receptor family` = grep("TNFR", colnames(rna), value = TRUE),
  `Fc receptors` = c("FCGR1A", "FCGR2A", "FCGR3A", "FCER1A", "FCER1G", "FCGR2B", "FCGR4"),
  `TCR receptor` = c("TRA", "TRB", "TRG", "TRD", "TRBC1", "TRBC2",
                     "CD3D", "CD3E", "CD3G", "CD247", "CD3Z",
                     "CD4", "CD8A", "CD8B"),
  `BCR complex` = c("IGLL1", "CD79A", "CD79B", "IGHD", "IGHG1",
                    "CD19", "CR2", "CD81"),
  `Generic co-receptors` = c("CD21", "CD28", "CTLA4", "PDCD1", "ICOS", "CD40", "CD80", "CD86", "PDCD1LG2"),
  `ISGs` = c("IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
             "ISG15", "RSAD2", "IFI27", "IFI44", "GBP1", "TRIM25", "EIF2AK2", 
             "CXCL10", "STAT1", "STAT2"),
  `TNF response` = c("TRAF2", "TRAF6", "TRADD", "CASP8", "FADD", "RIPK1", "RIPK3", "NFKB1", "RELA",
                     "IKBKB", "CHUK", "TNFAIP3", "BIRC2", "BIRC3", "MAP3K7", "MAPK8", "MAPK14", "CYLD"),
  `IRF family` = c("IRF3", "IRF7", "IRF9", "IRF1", "IRF5", "IRF4"),
  `JAK STAT signaling` = c("STAT1", "STAT3", "STAT4", "STAT5A", "STAT5B", "JAK1", "JAK2", "SOCS1", "SOCS3", "PTPN6"),
  NF_kB_signaling = c("NFKB1", "NFKB2", "RELA", "RELB", "CHUK", "IKBKB", "IKBKG", "TRAF6", "TAB1", "TAB2"),
  MAPK_pathway = c("MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK14", "MAP2K1", "MAP2K2", "MAP3K7", "TAOK1"),
  PI3K_AKT_mTOR = c("PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "AKT1", "AKT2", "MTOR", "RPTOR", "MLST8"),
  cAMP_PKA_signaling = c("ADCY3", "ADCY9", "PRKACA", "PRKACB", "PRKACG", "CREB1"),
  `T cell signaling` = c("LCK", "ZAP70", "LAT", "LCP2", "ITK", "PLCG1", "GRB2", "VAV1", "NFATC1", "CD28", "CSK"),
  `B cell signaling` = c("SYK", "BTK", "BLNK", "PLCG2", "PIK3CD", "VAV1", "GRB2", "AKT1", "CD19", "FCGR2B"),
  Pyroptosis = c("CASP1", "CASP3", "CASP4", "CASP6", "CASP8", "CASP9", "GPX4", "GSDMB", "GSDMD", "IL18", "IL1B"),
  Inflammasome = c("PLCG1", "PRKACA", "PYCARD", "SCAF11", "TIRAP"),
  `Cytotoxic mediators` = c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1"),
  `Cytokines Chemokines` = c("IL32", "IL23A", "IL1B", "IL18", "IL15", "IL16", "TGFB1", "TGFBI", "CCL5", "CCL4", "CCL3"),
  `TNF end products` = c("TNFAIP1", "TNFAIP2", "TNFAIP3", "TNFAIP8", "TNFAIP8L1", "TNFAIP8L2"),
  `Antimicrobial proteins` = c("LYZ", "S100A8", "S100A9"),
  `Acute phase` = c("FTH1", "FTL", "SERPINA1", "FTH1")
)
sum((unlist(categories) %>% unique()) %in% colnames(rna))

# Flatten the categories into a data frame
immune_df <- data.frame()
for (cat in names(categories)){
  if (length(categories[[cat]])>0){
    if ( any(categories[[cat]] %in% immune_genes$Receptors)){
      immune_df <- rbind(immune_df,
                         data.frame(var=categories[[cat]],
                                    PreciseCategory=cat,
                                    MacroCategory="Receptors"))
    } else if ( any(categories[[cat]] %in% immune_genes$Downstream)) {
      immune_df <- rbind(immune_df,
                         data.frame(var=categories[[cat]],
                                    PreciseCategory=cat,
                                    MacroCategory="Downstream"))
    } else {
      immune_df <- rbind(immune_df,
                         data.frame(var=categories[[cat]],
                                    PreciseCategory=cat,
                                    MacroCategory="EndProduct"))
    }
  }
}
df.plot <- left_join(immune_df, icc.df[icc.df$omic=="RNA",], by="var") %>% na.omit()
df.plot.summ <- df.plot %>% group_by(PreciseCategory) %>% 
  summarise(median_icc=median(icc.icc), n=n()) %>% as.data.frame()
df.plot$PreciseCategory <- factor(df.plot$PreciseCategory, levels=df.plot.summ$PreciseCategory[order(df.plot.summ$median_icc, decreasing=T)])

p1 <- ggplot(df.plot, aes(x=PreciseCategory, y=icc.icc, color=MacroCategory)) +
  geom_boxplot(outlier.shape=NA) +
  geom_text(data=df.plot.summ, aes(label=paste0("n=", n), y=1.2), color="black", size=3) +
  geom_hline(yintercept=median(icc.df$icc.icc[icc.df$omic=="RNA"]), linetype=2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, vjust=0.5, hjust=0), 
        plot.title = element_text(hjust = 0.5), legend.position="none") +
  xlab("") + ylab("ICC") + 
  ylim(0,1.3) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2 <- ggplot(df.plot, aes(x=PreciseCategory, y=icc.sex_b, color=MacroCategory)) +
  geom_boxplot(outlier.shape=NA) +
  #geom_text(data=df.plot.summ, aes(label=paste0("n=", n), y=1.2), color="black", size=3) +
  geom_hline(yintercept=median(icc.df$icc.sex_b[icc.df$omic=="RNA"]), linetype=2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, vjust=0.5, hjust=0), 
        plot.title = element_text(hjust = 0.5), legend.position="none") +
  xlab("") + ylab("Sex effect") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p3 <- ggplot(df.plot, aes(x=PreciseCategory, y=icc.age_b, color=MacroCategory)) +
  geom_boxplot(outlier.shape=NA) +
  #geom_text(data=df.plot.summ, aes(label=paste0("n=", n), y=1.2), color="black", size=3) +
  geom_hline(yintercept=median(icc.df$icc.age_b[icc.df$omic=="RNA"]), linetype=2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, vjust=0.5, hjust=0), 
        plot.title = element_text(hjust = 1), legend.position="bottom") +
  xlab("") + ylab("Age effect") 


ggarrange(p1,p2,p3, ncol=1)

pdf("custom-genes-ICC.pdf", width=15, height=8)
ggarrange(p1,p2,p3, ncol=1)
dev.off()

df.plot$MacroCategory <- factor(df.plot$MacroCategory ,
                                levels=df.plot %>% group_by(MacroCategory) %>% summarise(median_icc=median(icc.icc)) %>% arrange(desc(median_icc)) %>% pull(MacroCategory))
p <- ggplot(df.plot, aes(x=MacroCategory, y=icc.icc)) + 
  geom_violin() +
  geom_boxplot(width=0.2, outlier.shape=NA) +
  theme_classic() +
  xlab("")
p

pdf("custom-genes-ICC-macrocategory.pdf", width=3, height=4)
print(p)
dev.off()


#Protein macro category based on HPA annotation ####
receptor.prot <- hpa$Gene[grep("Predicted membrane proteins", hpa$Protein.class, ignore.case = T)] %>% unique() %>% intersect(colnames(protein))
intra.prot <-  hpa$Gene[grep("Predicted intracellular proteins", hpa$Protein.class, ignore.case = T)] %>% unique() %>% intersect(colnames(protein))
secreted.prot <-  hpa$Gene[grep("Predicted secreted proteins", hpa$Protein.class, ignore.case = T)] %>% unique() %>% intersect(colnames(protein))







# Plot examples of stable genes ####

plot_stable_gene <- function(g){
  #gather data of gene g
  df.plot <- rna.log[,g, drop=F] %>% scale() %>% as.data.frame()
  colnames(df.plot) <- "gene"
  df.plot$ind <- gsub("\\:.*", "", rownames(df.plot))
  df.plot$visit <- gsub(".*\\:", "", rownames(df.plot)) %>% as.numeric()
  df.plot <- df.plot %>% filter(visit %in% c("1", "2", "3", "4"))
  
  #divide in groups based on first visit
  df.plot.1 <- df.plot %>% filter(visit=="1") %>% 
    mutate(group=as.character(cut(gene, breaks = quantile(gene, probs = seq(0, 1, 0.25), na.rm = TRUE),
                                  include.lowest = TRUE, labels=F))) %>% 
    select(c(ind, group)) 
  
  df.plot$group <- df.plot.1$group[match(df.plot$ind, df.plot.1$ind)]
  df.plot.summ <- df.plot %>% group_by(group,visit) %>% summarise(median_gene=median(gene)) %>% as.data.frame()
  
  #plot
  p <- ggplot(df.plot, aes(x=visit, y=gene, color=group)) +
    geom_line(alpha=0.2, aes(group=ind)) +
    geom_point(alpha=0.2) +
    geom_line(data=df.plot.summ, aes(group=group, x=visit, y=median_gene), linewidth=1, linetype=2) +
    geom_point(data=df.plot.summ, aes(x=visit, y=median_gene), size=3) +
    scale_color_manual(values=c("1"="#EF8354", "2"="#7FC6A4", "3"="#0077B6", "4"="#750D37")) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
    xlab("Visit") + ylab("") + ggtitle(g) +
    xlim(1,4)
  return(p)
}

for (cat in names(categories)){
  g.vec <- categories[[cat]] %>% intersect(colnames(rna.log))
  if (length(g.vec)>0){
    p.list <- list()
    for (g in g.vec){
      
      p <- plot_stable_gene(g)
      
      pdf(paste0("examples-stable/", g, ".pdf"), width=4, height=3)
      print(p)
      dev.off()
      p.list[[g]] <- p
    }
  }
  p <- ggarrange(plotlist = p.list)
  print(p)
}

stable_genes <- df.plot.rna %>% filter(icc.icc>0.9) %>% pull(var)
for (g in stable_genes){
  p <- plot_stable_gene(g)
  pdf(paste0("examples-stable/", g, ".pdf"), width=4, height=3)
  print(p)
  dev.off()
  
}

# Plot examples of stable frequencies ####
for (pop in colnames(cytof.group)){
  
  #gather data of gene g
  df.plot <- cytof.group[,pop, drop=F] %>% as.data.frame()
  colnames(df.plot) <- "pop"
  df.plot$ind <- gsub("\\:.*", "", rownames(df.plot))
  df.plot$visit <- gsub(".*\\:", "", rownames(df.plot)) %>% as.numeric()
  df.plot <- df.plot %>% filter(visit %in% c("1", "2", "3", "4"))
  
  #divide in groups based on first visit
  df.plot.1 <- df.plot %>% filter(visit=="1") %>% 
    mutate(group=as.character(cut(pop, breaks = quantile(pop, probs = seq(0, 1, 0.25), na.rm = TRUE),
                                  include.lowest = TRUE, labels=F))) %>% 
    select(c(ind, group)) 
  
  df.plot$group <- df.plot.1$group[match(df.plot$ind, df.plot.1$ind)]
  df.plot.summ <- df.plot %>% group_by(group,visit) %>% summarise(median_pop=median(pop)) %>% as.data.frame()
  
  #plot
  p <- ggplot(df.plot, aes(x=visit, y=pop, color=group)) +
    geom_line(alpha=0.2, aes(group=ind)) +
    geom_point(alpha=0.2) +
    geom_line(data=df.plot.summ, aes(group=group, x=visit, y=median_pop), linewidth=1, linetype=2) +
    geom_point(data=df.plot.summ, aes(x=visit, y=median_pop), size=3) +
    scale_color_manual(values=c("1"="#EF8354", "2"="#7FC6A4", "3"="#0077B6", "4"="#750D37")) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
    xlab("Visit") + ylab("") + ggtitle(pop) +
    xlim(1,4)
  
  pdf(paste0("examples-stable/", pop, ".pdf"), width=4, height=3)
  print(p)
  dev.off()
}


# Plot examples of stable proteins ####
for (prot in top.prot){
  
  #gather data of gene g
  df.plot <- protein[,prot, drop=F] %>% as.data.frame()
  colnames(df.plot) <- "pop"
  df.plot$ind <- gsub("\\:.*", "", rownames(df.plot))
  df.plot$visit <- gsub(".*\\:", "", rownames(df.plot)) %>% as.numeric()
  df.plot <- df.plot %>% filter(visit %in% c("1", "2", "3", "4"))
  
  #divide in groups based on first visit
  df.plot.1 <- df.plot %>% filter(visit=="1") %>% 
    mutate(group=as.character(cut(pop, breaks = quantile(pop, probs = seq(0, 1, 0.25), na.rm = TRUE),
                                  include.lowest = TRUE, labels=F))) %>% 
    select(c(ind, group)) 
  
  df.plot$group <- df.plot.1$group[match(df.plot$ind, df.plot.1$ind)]
  df.plot.summ <- df.plot %>% group_by(group,visit) %>% summarise(median_pop=median(pop)) %>% as.data.frame()
  
  #plot
  p <- ggplot(df.plot, aes(x=visit, y=pop, color=group)) +
    geom_line(alpha=0.2, aes(group=ind)) +
    geom_point(alpha=0.2) +
    geom_line(data=df.plot.summ, aes(group=group, x=visit, y=median_pop), linewidth=1, linetype=2) +
    geom_point(data=df.plot.summ, aes(x=visit, y=median_pop), size=3) +
    scale_color_manual(values=c("1"="#EF8354", "2"="#7FC6A4", "3"="#0077B6", "4"="#750D37")) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
    xlab("Visit") + ylab("") + ggtitle(prot) +
    xlim(1,4)
  
  pdf(paste0("examples-stable/", gsub("/", "-", prot), ".pdf"), width=4, height=3)
  print(p)
  dev.off()
  
}

#Plasma protein-immune composition correlation ####
common.samples <- intersect(rownames(protein), rownames(cytof.group))

cor.prot.pop <- cor(protein[common.samples,], cytof.group[common.samples,], method="spearman", use="complete.obs")

df.cor <- data.frame()
for (pop in colnames(cytof.group)){
  r <- cor.prot.pop[,pop]  
  df.cor <- rbind(df.cor, data.frame(r=r, protein=colnames(protein), pop=pop))
}
df.cor <- df.cor %>% filter(r>0.2) %>% arrange(pop, desc(r))

df.cor$pop <- gsub("B_Cells", "Memory B cells", df.cor$pop)
df.cor$pop <- gsub("Central_Memory_CD4", "T-cm CD4+", df.cor$pop)
df.cor$pop <- gsub("Central_Memory_CD8", "T-cm CD8+", df.cor$pop)
df.cor$pop <- gsub("Effector_Memory_CD8", "T-em CD8+", df.cor$pop)
df.cor$pop <- gsub("Monocytes_classical", "Classical monocytes", df.cor$pop)
df.cor$pop <- gsub("Monocytes_intermediate", "Intermediate monocytes", df.cor$pop)
df.cor$pop <- gsub("Monocytes_nonclassical", "Nonclassical monocytes", df.cor$pop)
df.cor$pop <- gsub("Myeloid_DC", "Myeloid DCs", df.cor$pop)
df.cor$pop <- gsub("NK_Cells", "NK cells", df.cor$pop)
df.cor$pop <- gsub("Naive_B_Cells", "Naive B cells", df.cor$pop)
df.cor$pop <- gsub("Naive_CD4", "Naive CD4+ T cells", df.cor$pop)
df.cor$pop <- gsub("Naive_CD8", "Naive CD8+ T cells", df.cor$pop)
df.cor$pop <- gsub("TEMRA_CD4", "TEMRA CD4+ T cells", df.cor$pop)
df.cor$pop <- gsub("TEMRA_CD8", "TEMRA CD8+ T cells", df.cor$pop)


# Supplementary Table 2 ####

#CyTOF ICC
df.cytof <- icc.df %>% filter(omic=="CyTOF") %>% arrange(desc(icc.icc))
colnames(df.cytof) <- c("ICC", "Sex coefficient", "Age coefficient", "Sex p-value", "Age p-value", "Population")

#RNA ICC
df.rna <- icc.df %>% filter(omic=="RNA") %>% arrange(desc(icc.icc))
colnames(df.rna) <- c("ICC", "Sex coefficient", "Age coefficient", "Sex p-value", "Age p-value", "Gene")

#protein ICC
df.protein <- icc.df %>% filter(omic=="Protein") %>% arrange(desc(icc.icc))
colnames(df.protein) <- c("ICC", "Sex coefficient", "Age coefficient", "Sex p-value", "Age p-value", "Protein")

#Manual annotation
df <- immune_df 
colnames(df) <- c("Gene", "Immune pathway", "Protein class")
df <- df %>% filter(!`Immune pathway` %in% c("cGAS_STING", "IFN receptors", "Generic co-receptors", "ISGs", "TNF response", "cAMP_PKA_signaling", "Cytokines Chemokines"))

wb <- createWorkbook()
addWorksheet(wb, "Immune frequency ICC")
addWorksheet(wb, "RNA-seq ICC")
addWorksheet(wb, "Proteomics ICC")
addWorksheet(wb, "Manual annotation")
addWorksheet(wb, "Corr proteome-immune comp")
writeData(wb, "Immune frequency ICC", df.cytof[,1:6])
writeData(wb, "RNA-seq ICC", df.rna[,1:6])
writeData(wb, "Proteomics ICC", df.protein[,1:6])
writeData(wb, "Manual annotation", df)
writeData(wb, "Corr proteome-immune comp", df.cor)

saveWorkbook(wb, "Supplementary Table 2.xlsx", overwrite = TRUE)
