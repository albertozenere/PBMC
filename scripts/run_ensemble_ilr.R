
pack_R <- c("dplyr","ggplot2","ggrepel","umap", "edgeR",
            "RColorBrewer", "pheatmap", "tidyverse",
            "igraph", "ForceAtlas2", "biomaRt", 
            "org.Hs.eg.db", "ggpubr",
            "clusterProfiler", "RANN", "dbscan", "cluster", 
            "corrplot", "igraph",
            "ggbeeswarm", "glmnet", "metafor")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

setwd("C:/Users/albze08/Desktop/postDoc/wellness-rna/")

set.seed(1)


# Load ####

# Annotation of individuals
anno <- read.csv("C:/Users/albze08/Desktop/postDoc/genome-protein/data/WELLNESS/Wellness/Data/Wellness_barcodes.txt", sep="\t", header=T)
anno <- anno[anno$Sample.type=="Helblod",]

#Proteomics
protein_data <- read.csv("C:/Users/albze08/Desktop/postDoc/genome-protein/data/WELLNESS/Wellness/Data/Proteome/wellness_norm_final_794_wj.txt", sep="\t", header=T)
protein_data$subject_id <- anno$Subject[match(protein_data$sample, anno$Wellness.id)]
protein_data$id <- paste0(protein_data$subject_id, ":", protein_data$visit)
rownames(protein_data) <- protein_data$id
protein_data <- subset(protein_data, select=-c(iid, sample, visit, id, subject_id))

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
filter.genes <- filterByExpr(t(rna), min.count = 1, min.total.count = 10, large.n = 10, min.prop = 0.7)
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

# Compositional regression functions ####
#ILR transformation, Hron et al. J. Applied Statistics, 2009
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

circ.shift <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

lm_ilr <- function(y,X.ilr,covar=NULL){
  
  D <- length(X.ilr)
  
  #init output
  lm.results <- data.frame(estimate=rep(0,D), pvalue=rep(0,D))
  
  #do lm
  for (n in 1:D){  
    
    #nth ilr
    X.ilr.n <- X.ilr[[n]]
    var.n <- colnames(X.ilr.n)[1]
    rownames(lm.results)[n] <- var.n
    
    #lm
    df <- cbind(y,X.ilr.n, covar)
    lmFit <- lm(formula(df), df) %>% summary()
    
    #store
    lm.results[var.n,] <- lmFit$coefficients[var.n, c("Estimate", "Pr(>|t|)")]
    
  }
  
  return(lm.results)
}


# Estimate gene-pop associations using ILR+ensemble ####
common.samples <- intersect(rownames(cytof), rownames(rna.log))

visit.uniq <- 1:4 %>% as.character()
df.net <- vector(mode="list", ncol(rna.log))
names(df.net) <- colnames(rna.log)
for (n.gene in 1:ncol(rna.log)){
  print(n.gene)
  #Gather data
  g <- colnames(rna.log)[n.gene]
  y <- rna.log[common.samples,g] %>% scale(scale = F)
  X.freq <- cytof.group[common.samples,]
  covar <- data.frame(visit=gsub(".*\\:", "", common.samples),
                      sex=ifelse(metadata$Gender[match(common.samples, metadata$id)]=="m",1,0),
                      age=metadata$Age_at_Visit[match(common.samples, metadata$id)])
  
  
  #Divide by visit
  df.gene <- data.frame()
  for (visit in visit.uniq){
    idx <- which(covar$visit == visit)
    
    y.visit <- y[idx]
    X.freq.visit <- X.freq[idx,]
    covar.visit <- covar[idx,]
    
    #One lm per transform
    D <- ncol(X.freq.visit)
    
    
    X.ilr <- matrix(0, nrow(X.freq.visit), D-1)
    
    for (n in 1:ncol(X.freq.visit)){
      
      idx.n <- circ.shift(1:D, n-1) #D
      vec.names <- colnames(X.freq.visit)[idx.n[1:(D-1)]] #D-1
      colnames(X.ilr) <- vec.names
      
      
      for (m in 1:nrow(X.ilr)){ #Transform
        vec <- X.freq.visit[m,] %>% as.numeric()
        vec.shift <- vec[idx.n] #D
        X.ilr[m,] <- ilr.transform(vec.shift) #D-1
      }
      df <- data.frame(y=y.visit, X.ilr, covar.visit[,c("sex", "age")])
      
      lm.visit <- lm(formula(df), df) %>% summary() %>% coef() %>% as.data.frame()
      
      df.gene <- rbind(df.gene,
                       data.frame(lm.visit[vec.names[1],], visit=visit, pop=vec.names[1]))
      
    } 
  }
  
  # Ensemble by using meta-analysis
  out.ensemble <- df.gene %>%
    group_by(pop) %>%
    summarise(meta = list(
      tryCatch(
        rma(yi = Estimate, sei = Std..Error, data = pick(Estimate, Std..Error)),
        error = function(e) NULL  
      ))) %>%
    mutate(result=setNames(meta, pop)) %>% pull(result) %>% 
    lapply(function(x){data.frame(estimate=x$b[1],
                                  pval=x$pval,
                                  ci_lower=x$ci.lb,
                                  ci_upper=x$ci.ub)}) %>%
    bind_rows(.id = "source")
  df.net[[n.gene]] <- data.frame(out.ensemble, gene=g)
  
}

saveRDS(df.net, "Ensemble_ILR-306025.RDS")

