rm(list = ls())

# Setup ####
path <- .libPaths()
newpath <- "C:/Users/albze08/Desktop/R/mQTL"
.libPaths(newpath)


pack_R <- c("dplyr","ggplot2","ggrepel","umap", "edgeR", "Rtsne",
            "RColorBrewer", "pheatmap", "tidyverse", "vegan",
            "igraph", "ForceAtlas2", "biomaRt", 
            "org.Hs.eg.db", "ggpubr", "glmnet", "Rtsne",
            "clusterProfiler", "RANN", "dbscan", "cluster", 
           "corrplot", "igraph")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

setwd("C:/Users/albze08/Desktop/postDoc/PBMC/")

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

# Fix CyTOF ####
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


# ILR transform ####
nthroot = function(x,n) {
  (abs(x)^(1/n))*sign(x)
}

ilr.transform <- function(vec){ #Hron et al. J. Applied Statistics, 2009 Eq.4
  D <- length(vec)
  vec.ilr <- rep(0, D-1)
  for (i in 1:length(vec.ilr)){
    vec.ilr[i] <- sqrt((D-i)/(D-i+1)) * log((vec[i])/nthroot(prod(vec[(i+1):D]), D-i))
  }
  return(vec.ilr)
}

ilr.matrix <- function(mat){
  mat.ilr <- matrix(0, nrow(mat), ncol(mat)-1)
  rownames(mat.ilr) <- rownames(mat)
  for (n in 1:nrow(mat)){
    mat.ilr[n,] <- mat[n,] %>% as.numeric() %>% ilr.transform()
  }
  return(mat.ilr)
}


clr.transform <- function(vec){
  D <- length(vec)
  vec.ln <- log(vec)
  vec.clr <- vec.ln - (sum(vec.ln))/D
  
  # geometric_mean <- exp(mean(log(vec)))
  # vec.clr <- log(vec / geometric_mean)
  
  return(vec.clr)
}

clr.matrix <- function(mat){
  mat.clr <- mat
  for (n in 1:nrow(mat)){
    mat.clr[n,] <- mat[n,] %>% as.numeric() %>% clr.transform()
  }
  return(mat.clr)
}


condense_PCA <- function(df, Npca=50){
  pca_result <- prcomp(df, center = F, scale. = F)
  df_pca <- pca_result$x[, 1:Npca]
  return(df_pca)
}

# UMAP and silhouettte ####
umap_fun <- function(X, color.UMAP, method="UMAP"){
  
  if (method=="UMAP"){
    umap_df <- umap(X, n_neighbors = 10, n_components = 2)
    df.plot <- umap_df$layout %>% as.data.frame()
    colnames(df.plot)[1:2] <- c("UMAP1", "UMAP2")
  } else if (method=="MDS") {
    mds_df <- metaMDS(X, k = 2, trymax = 10)
    df.plot <- mds_df$points %>% as.data.frame()
    colnames(df.plot)[1:2] <- c("UMAP1", "UMAP2")
  } else if (method=="tSNE"){
    tsne_df <- Rtsne(X, dims = 2, perplexity = 5, verbose = F, max_iter = 500)
    df.plot <- tsne_df$Y[,1:2] %>% as.data.frame()
    rownames(df.plot) <- rownames(X)
    colnames(df.plot)[1:2] <- c("UMAP1", "UMAP2")
  } else if (method=="PCA"){
    pca_df <- prcomp(X, center = F, scale. = F)
    df.plot <- pca_df$x[, 1:2] %>% as.data.frame()
    rownames(df.plot) <- rownames(X)
    colnames(df.plot)[1:2] <- c("UMAP1", "UMAP2")
  } else {
    print("Method not supported")
    stopifnot(1<0)
  }
  
  df.plot$subject <- metadata$subject_id[match(rownames(df.plot), metadata$id)] %>% as.character()
  df.plot$visit <- metadata$visit[match(rownames(df.plot), metadata$id)]
  
  #Plot
  p <- ggplot(df.plot, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(size=1, aes(color=subject)) + 
    geom_line(aes(group=subject, color=subject), linetype=2) +
    scale_color_manual(values=color.UMAP) + 
    theme_classic() + theme(text=element_text(size=12), legend.position="none") +
    ylab("")
  p
  
  return(list(p=p, df=df.plot))
}

#Silhouette
silhouette_fun <- function(df.plot, color.UMAP, method.dist="euc"){
  ind <- gsub("\\:.*", "", rownames(df.plot))
  
  sil <- silhouette(as.numeric(ind), coda.base::dist(df.plot, method=method.dist))

  df.plot <- fviz_silhouette(sil)$data
  df.plot$cluster <- as.character(df.plot$cluster)
  
  df <- df.plot %>% group_by(cluster) %>% summarise(avg.score=mean(sil_width)) %>% as.data.frame()
  ind.ord <- as.character(df$cluster)[order(df$avg.score, decreasing=T)]
  
  df.plot$sample <- paste0(df.plot$cluster, df.plot$name) 
  sample.ord <- character(0)
  for (ind in ind.ord){
    sample.ord <- c(sample.ord, df.plot$sample[df.plot$cluster==ind])
  }
  df.plot$sample <- factor(df.plot$sample, levels=sample.ord)
  
  p <- ggplot(df.plot, aes(x=sample, y=sil_width, fill=cluster, color=cluster)) + 
    geom_bar(stat="identity", position="dodge") +
    geom_hline(yintercept=0) +
    geom_hline(aes(yintercept=mean(sil_width)), linetype=2, color="red") +
    theme_classic2() +
    theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_fill_manual(values=color.UMAP) + scale_color_manual(values=color.UMAP) + 
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=12), legend.position="none") +
    xlab("Sample") + ylab("") +
    ylim(c(-1,1))
  
  return(p)
}


sample.cytof.rna <- intersect(rownames(cytof), rownames(rna)) 
ind.uniq <- gsub("\\:.*", "", sample.cytof.rna) %>% unique()
color.UMAP <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] %>% sample(length(ind.uniq))
color.UMAP <- setNames(color.UMAP, ind.uniq)

cytof.transf <- cytof.group[sample.cytof.rna,] %>% as.data.frame()
rna.pca <- rna.log[sample.cytof.rna,] %>% condense_PCA(Npca=50) %>% scale() %>% as.data.frame()

method <- "tSNE"
out.cytof <- umap_fun(cytof.transf %>% scale(), color.UMAP, method=method)
out.rna <- umap_fun(rna.pca, color.UMAP, method=method)

p1 <- silhouette_fun(cytof.transf, color.UMAP, method.dist = "ait") + ylab("Silhouette score")
p2 <- silhouette_fun(rna.pca, color.UMAP)

p <- ggarrange(out.cytof$p + ylab("UMAP2"), out.rna$p,
               p1, p2, nrow=2, ncol=2)
p

pdf("UMAP-silhouette-2.pdf", width=6, height=6)
print(p)
dev.off()


#PBMC omics combined
pbmc.data <- cbind(cytof.transf %>% scale(), rna.pca)
out.pbmc <- umap_fun(pbmc.data, color.UMAP, method=method)

pdf("UMAP-PBMC.pdf", width=3, height=3)
out.pbmc$p
dev.off()


# Proteomics t-SNE
protein.pca <- protein %>% makeX(na.impute = T) %>% condense_PCA(Npca=50)
out.protein <- umap_fun(protein.pca, color.UMAP, method=method)


pdf("UMAP-all.pdf", width=6, height=6)
ggarrange(out.cytof$p, out.rna$p, out.pbmc$p, out.protein$p, nrow=2, ncol=2)
dev.off()



# Show genes most associated with the clustering
coord.rna <-  out.rna$df

cor.rna.coord <- data.frame(cor1=cor(rna.log[sample.cytof.rna,], coord.rna[sample.cytof.rna,"UMAP1"], method="spearman"),
                            cor2=cor(rna.log[sample.cytof.rna,], coord.rna[sample.cytof.rna,"UMAP2"], method="spearman"),
                            gene=colnames(rna.log)) %>% 
  mutate(cor_tot=abs(cor1)+abs(cor2)) %>% arrange(desc(cor_tot))

cor.rna.coord %>% arrange(desc(abs(cor1))) %>% pull(gene) %>% head(n=20)
cor.rna.coord %>% arrange(desc(abs(cor2))) %>% pull(gene) %>% head(n=20)


# Compare distances between inter-ind visit 1-4, inter-ind visit 1-4 vs 5-6, inter-ind ####
rna.pca <- rna.log %>% condense_PCA(Npca=50) %>% scale() %>% as.data.frame()

sample_id <- rownames(rna.log)
visit <- gsub(".*\\:", "", sample_id)
ind <- gsub("\\:.*", "", sample_id)  
ind.uniq <- unique(ind)

df.plot <- data.frame()
for (i in ind.uniq){
  
  sample.ind.14 <- paste0(i, ":", 1:4)
  sample.ind.56 <- paste0(i, ":", 5:6)
  sample.ind <- c(sample.ind.14, sample.ind.56)
  
  idx.14 <- which(visit %in% c("1", "2", "3", "4") & ind==i)
  idx.56 <- which(visit %in% c("5", "6") & ind==i)
  
  #distance between visits 1-4
  d <- dist(rna.pca[idx.14,])
  d.vec <- as.matrix(d)[upper.tri(as.matrix(d))]

  df.plot <- rbind(df.plot, data.frame(d=d.vec, ind=i, group="intra_14"))
  
  #distance between visit 1-4 and 5-6
  d <- dist(rbind(rna.pca[idx.14,], rna.pca[idx.56,])) %>% as.matrix()
  d.vec <- d[sample.ind.14, sample.ind.56] %>% c()
  
  df.plot <- rbind(df.plot, data.frame(d=d.vec, ind=i, group="intra_56"))
  
  #distance to other individuals
  d <- dist(rna.pca) %>% as.matrix()
  
  for (i.other in setdiff(ind, i)){
    df.plot <- rbind(df.plot,
                     data.frame(d=d[sample.ind, grep(i.other, colnames(d))] %>% colMeans() %>% mean(), ind=i, group="inter"))
    
  }
}
df.plot.summ <- df.plot %>% group_by(ind, group) %>% summarise(d_med=median(d), d_min=min(d), d_max=max(d)) 

#df.plot$ind <- factor(df.plot$ind, levels=df.plot.summ %>% filter(group=="intra_14_dup") %>% arrange(desc(d_med)) %>% pull(ind))
p <- ggplot() + 
  geom_point(data=df.plot.summ, aes(y=ind,x=d_med,color=group), size=3) +
  geom_segment(data=df.plot.summ, aes(y=ind, x=d_min, xend=d_max, color=group), linewidth=1) +
  #geom_point(data=df.plot.summ, aes(y=ind,x=d_med), shape=18, size=3) + 
  scale_color_manual(values=c(intra_14="#FFCCD1", intra_56="#C9B7DC", inter="gray77", intra_14_dup="#FF828F", intra_56_dup="#AF86DB")) +
  theme_classic() +
  theme(legend.position = "none") +  
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab("Euclidean distance") + ylab("Participants")
p

pdf("dist-samples.pdf", height=10, width=6)
print(p)
dev.off()



#distance between samples at visits 6 and 1-4
sample_6 <- grep("\\:6", sample_id, value=T)
sample_14 <- grep("\\:1|\\:2|\\:3|\\:4", sample_id, value=T) 
df.dist.6 <- data.frame()
for (s1 in sample_6){
  ind1 <- gsub("\\:.*", "", s1)
  for (s2 in sample_14){
    ind2 <- gsub("\\:.*", "", s2)
    
    d <- sqrt(sum((rna.pca[s1,] - rna.pca[s2,])^2))
    if (ind1==ind2){
      df.dist.6 <- rbind(df.dist.6,
                         data.frame(ind=ind1, d=d, group="intra") )
    } else {
      df.dist.6 <- rbind(df.dist.6,
                         data.frame(ind=ind1, d=d, group="inter") )
    }
  }
}
df.dist.6$sample <- 1:nrow(df.dist.6)


df.dist.6.inter <- df.dist.6 %>% filter(group=="inter") %>% group_by(ind) %>% summarise(d_min=min(d), d_max=max(d), d_med=median(d))
df.dist.6.intra <- df.dist.6 %>% filter(group=="intra") %>% group_by(ind) %>% summarise(d_min=min(d), d_max=max(d), d_med=median(d))
df.plot <- rbind(df.dist.6.intra %>% mutate(group="intra"),
                 df.dist.6.inter %>% mutate(group="inter"))
df.plot$ind <- factor(df.plot$ind, levels=df.dist.6.intra %>% arrange(d_med) %>% pull(ind))

p <- ggplot(df.plot) + 
  geom_segment(aes(y=ind, x=d_max, xend=d_min, color=group)) + 
  geom_point(aes(y=ind, x=d_med, color=group)) +
  theme_classic() +
  theme(legend.position = "none") +  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  +
  xlab("Euclidean distance") + ylab("Participants") +
  coord_flip()
p

cairo_pdf("dist-RNA-ind.pdf", height=3.5, width=6)
print(p)
dev.off()



df.plot.summ$group <- factor(df.plot.summ$group, levels=c("intra_14", "intra_56", "inter"))
p <- ggplot(df.plot.summ, aes(x=group, y=d_med)) +
  geom_violin() +
  geom_quasirandom() +
  stat_compare_means(method="wilcox.test", comparisons=list(c("intra_14", "intra_56"), c("intra_14", "inter"), c("intra_56", "inter")), label="p.signif") +
  theme_classic()
p

cairo_pdf("dist_RNA.pdf", width=3, height=3)
print(p)
dev.off()




# CyTOF CV ####
cytof$ind <- gsub("\\:.*", "", rownames(cytof))
inter.cv <- cytof %>% group_by(ind) %>% 
  summarise(across(everything(), mean)) %>%
  subset(select=-ind) %>% apply(2,function(n){sd(n)/abs(mean(n))}) 

intra.cv <- cytof %>% group_by(ind) %>%
  summarise(across(everything(), function(n){sd(n)/abs(mean(n))})) %>%
  subset(select=-ind) %>% apply(2, mean)

stopifnot(all(names(inter.cv)==names(intra.cv)))
df.plot <- data.frame(inter=inter.cv, intra=intra.cv)
df.plot$bin <- 3
df.plot$bin[df.plot$intra>2] <- 2
df.plot$bin[df.plot$intra>10] <- 1

df.plot$label <- ifelse(df.plot$inter>0.25, rownames(df.plot), NA)
p1 <- df.plot %>%  
  ggplot(aes(x=inter, y=intra, label=rownames(df.plot))) + 
  geom_point(size=2, fill="gray77", shape=21) + 
  geom_text() +
  geom_abline(slope=1, intercept=0) + 
  xlab("Inter-individual CV") + ylab("Intra-individual CV") + ggtitle("CyTOF") + 
  theme_classic() + ylim(c(-1,50)) +
  facet_grid(bin ~ ., scale='free_y') +
  scale_y_continuous() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.y = element_blank()) 
p1

print(paste0("There are ", sum(df.plot$inter<df.plot$intra, na.rm=T), " populations with inter-CV < intra-CV (BAD) and ",
             sum(df.plot$inter>df.plot$intra, na.rm=T), " populations with inter-CV > intra-CV (GOOD)"))

colMeans(cytof[, rownames(df.plot)[df.plot$inter<df.plot$intra]]) %>% sort()


# RNA CV ####
rna.log$ind <- gsub("\\:.*", "", rownames(rna.log))
inter.cv <- rna.log %>% group_by(ind) %>% 
  summarise(across(everything(), mean)) %>%
  subset(select=-ind) %>% apply(2,function(n){sd(n)/abs(mean(n))}) 

intra.cv <- rna.log %>% group_by(ind) %>%
  summarise(across(everything(), function(n){sd(n)/abs(mean(n))})) %>%
  subset(select=-ind) %>% apply(2, mean)

stopifnot(all(names(inter.cv)==names(intra.cv)))
df.plot <- data.frame(inter=inter.cv, intra=intra.cv) -> df.plot.rna
df.plot$label <- ifelse(df.plot$inter>0.25, rownames(df.plot), NA)
p2 <- ggplot(df.plot, aes(x=inter, y=intra)) + 
  geom_point(size=2, fill="gray77", shape=21) + 
  geom_text_repel(aes(label=label), size=2) + 
  geom_abline(slope=1, intercept=0, linetype=2) + 
  xlim(c(0,max(df.plot$inter))) + ylim(c(0,max(df.plot$intra))) +
  xlab("Inter-individual CV") + ylab("Intra-individual CV") + ggtitle("RNA-seq") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))
p2

pdf("CV-rna.pdf", height=3, width=3)
print(p2)
dev.off()

p
png("CV-rna.png", width = 3, height = 3, units = "in", res = 300)
print(p2)
dev.off()


print(paste0("There are ", sum(df.plot$inter<df.plot$intra, na.rm=T), " genes with inter-CV < intra-CV (BAD) and ",
             sum(df.plot$inter>df.plot$intra, na.rm=T), " genes with inter-CV > intra-CV (GOOD)"))




g <- rownames(df.plot)[df.plot$inter>df.plot$intra] 
gene_entrez <- bitr(g, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_ids <- gene_entrez$ENTREZID

go_result <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
pdf("inter-genes.pdf", width=6, height=3)
barplot(go_result)
dev.off()



# protein CV ####
protein$ind <- gsub("\\:.*", "", rownames(protein))
inter.cv <- protein %>% group_by(ind) %>% 
  summarise(across(everything(), mean)) %>%
  subset(select=-ind) %>% apply(2,function(n){sd(n)/abs(mean(n))}) 

intra.cv <- protein %>% group_by(ind) %>%
  summarise(across(everything(), function(n){sd(n)/abs(mean(n))})) %>%
  subset(select=-ind) %>% apply(2, mean)

stopifnot(all(names(inter.cv)==names(intra.cv)))
df.plot.protein <- df.plot <- data.frame(inter=inter.cv, intra=intra.cv)
df.plot$label <- ifelse(df.plot$inter>0.15, rownames(df.plot), NA)
p3 <- ggplot(df.plot, aes(x=inter, y=intra)) + 
  geom_point(size=2, fill="gray77", shape=21) + 
  geom_text_repel(aes(label=label), size=2) + 
  geom_abline(slope=1, intercept=0, linetype=2) + 
  xlim(c(0,max(df.plot$inter))) + ylim(c(0,max(df.plot$intra))) +
  xlab("Inter-individual CV") + ylab("Intra-individual CV") + ggtitle("Olink") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))
p3
















######### OLD ############

# Predict individuals from RNA-seq of visit 5-6 ####

# for each sample, choose the individual who has the avg. shorter distance from visits 1-4
Npca <- 50
ind.14 <- gsub("\\:.*", "", sample.cytof.rna) %>% unique()

idx.56 <- which(gsub(".*\\:", "", rownames(rna.log)) %in% c("5","6"))  %>% intersect(which( gsub("\\:.*", "", rownames(rna.log)) %in% ind.14))
rna.56 <- rna.log[idx.56,]

pca.14 <- prcomp(rna.log[sample.cytof.rna,], center = F, scale. = F)  
scores.14 <- predict(pca.14, newdata = rna.log[sample.cytof.rna,])[,1:Npca] %>% scale()
scores.56 <- predict(pca.14, newdata = rna.56)[,1:Npca] %>% scale()

df.pred <- data.frame(sample=rownames(rna.56))
df.pred$ind <- gsub("\\:.*", "", df.pred$sample)
df.pred$ind_pred <- NA
df.plot <- data.frame()
for (n in 1:nrow(df.pred)){
  
  x <- scores.56[n,]
  d <- data.frame(sample=rownames(scores.14), dist=NA)
  d$ind <- gsub("\\:.*", "", d$sample)
  for (k in 1:nrow(d)){
    d$dist[k] <- sqrt(sum((scores.14[k,] - x)^2))
  }
  
  d.summ <- d %>% group_by(ind) %>% summarise(mean_dist=mean(dist)) %>% as.data.frame()
  df.pred$ind_pred[n] <- d.summ$ind[which.min(d.summ$mean_dist)]
  
  df.plot <- rbind(df.plot, data.frame(sampleA=df.pred$sample[n], sampleB=d$sample, dist=d$dist))
}

df.plot$indA <- gsub("\\:.*", "", df.plot$sampleA)
df.plot$indB <- gsub("\\:.*", "", df.plot$sampleB)
df.plot.summ <- df.plot %>% filter(indA==indB) %>% group_by(sampleA, indB) %>% summarise(avg_dist=mean(dist)) %>% as.data.frame()
df.plot.summ <- df.plot %>% group_by(sampleA, indA, indB) %>% summarise(avg_dist=mean(dist)) %>% as.data.frame()

sample.ord <- df.plot.summ %>% filter(indA==indB) %>% arrange(avg_dist) %>% pull(sampleA)
df.plot.summ$sampleA <- factor(df.plot.summ$sampleA, levels=sample.ord)
p <- ggplot(df.plot.summ %>% filter(indA!=indB), aes(y=sampleA, x=avg_dist)) +
  geom_violin(fill="gray68") + 
  geom_point(data=df.plot.summ %>% filter(indA==indB), aes(y=sampleA, x=avg_dist), color="black", size=1, shape=1) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
p
pdf("dist-56.pdf", height=5, width=3)
print(p)
dev.off()

# Silhouette
X <- rbind(scores.14, scores.56)
ind <- gsub("\\:.*", "", rownames(X))

sil <- silhouette(as.numeric(ind), dist(X))

df.plot <- fviz_silhouette(sil)$data
df.plot$cluster <- as.character(df.plot$cluster)
df.plot$id <- rownames(X)[as.numeric(as.character(df.plot$name))]
df.plot$visit <- gsub(".*\\:", "", df.plot$id)

df <- df.plot %>% group_by(cluster) %>% summarise(avg.score=mean(sil_width)) %>% as.data.frame()
ind.ord <- as.character(df$cluster)[order(df$avg.score, decreasing=T)]

df.plot$sample <- paste0(df.plot$cluster, df.plot$name) 
sample.ord <- character(0)
for (ind in ind.ord){
  sample.ord <- c(sample.ord, df.plot$sample[df.plot$cluster==ind])
}
df.plot$sample <- factor(df.plot$sample, levels=sample.ord)

p <- ggplot(df.plot %>% filter(visit %in% c("5", "6")), aes(x=sample, y=sil_width, fill=cluster, color=cluster)) + 
  geom_bar(stat="identity", position="dodge") +
  geom_hline(yintercept=0) +
  geom_hline(aes(yintercept=mean(sil_width)), linetype=2, color="red") +
  theme_classic2() +
  theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_manual(values=color.UMAP) + scale_color_manual(values=color.UMAP) + 
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=12, family="Arial"), legend.position="none") +
  xlab("Sample") + ylab("") +
  ylim(c(-1,1))
p

cairo_pdf("UMAP-visit56.pdf", width=3, height=3)
print(p)
dev.off()

