
pack_R <- c("dplyr","ggplot2","ggrepel","umap", "edgeR",
            "RColorBrewer", "pheatmap", "tidyverse",
            "igraph", "ForceAtlas2", "biomaRt", "fmsb",
            "org.Hs.eg.db", "ggpubr", "RANN", "ForceAtlas2",
            "clusterProfiler", "RANN", "dbscan")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(1)


# Load ####

# Annotation of individuals
anno <- read.csv("data/Wellness_barcodes.txt", sep="\t", header=T)
anno <- anno[anno$Sample.type=="Helblod",]

# Proteomics
protein_data <- read.csv("data/wellness_norm_final_794_wj_alternative_anno.txt", sep="\t", header=T)
protein_data$subject <- anno$Subject[match(protein_data$sample, anno$Wellness.id)]
protein_data$id <- paste0(protein_data$subject, ":", protein_data$visit)
sample <- protein_data$id

protein.long <- protein_data[, c("id", "protein_name", "NPX")]
protein <- pivot_wider(protein.long, names_from=protein_name, values_from=NPX, values_fn=sum) %>% as.data.frame()
rownames(protein) <- protein$id
protein <- protein[,-1]

protein <- makeX(protein, na.impute = T)
rm(protein.long)

# metadata
metadata <- read.csv("data/complete.clinical.data.wellness.t2d.txt", sep="\t", header=T)
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

# Color palette ####
pop.palette <- c(
  "Monocytes_classical" = "#C24100", "Monocytes_intermediate"  = "#E66101", "Monocytes_nonclassical"  = "#FDB863",  
  "Basophils" = "#FFD92F", "Myeloid_DC" = "#E6AB02", "Plasmacytoid_DC" = "#B8A200",   
  "B_Cells" = "dodgerblue4", "Naive_B_Cells" = "dodgerblue2",  
  "Naive_CD4" = "#33A02C", "Central_Memory_CD4" = "#66C2A5","Effector_Memory_CD4" = "#1B9E77", "TEMRA_CD4" = "#41AE76",  
  "Naive_CD8" = "#FB9A99", "Central_Memory_CD8" = "#E31A1C", "Effector_Memory_CD8" = "#B2182B", "TEMRA_CD8" = "#67001F",  
  "NK_Cells" = "#CC78BC",  
  "ILCs" = "#B3B3B3"   
)

# Load network ####
cell.gene <- readRDS("Ensemble_ILR-306025.RDS") %>% bind_rows() %>% 
  filter(estimate>0) %>% 
  #mutate(pval.bon=pval*n()) %>% filter(pval.bon<0.05) %>% 
  mutate(pval.BH=p.adjust(pval, method="BH")) %>% filter(pval.BH<0.05) %>% 
  filter(gene %in% colnames(rna)) 
colnames(cell.gene)[1] <- "pop"
table(cell.gene$pop) %>% sort()


# Examples pop-gene ####

shared.sample <- intersect(rownames(cytof.group), rownames(rna.log))
df.plot <- data.frame(pop=c(cytof.group[shared.sample, "B_Cells"],
                            cytof.group[shared.sample, "B_Cells"],
                            cytof.group[shared.sample, "B_Cells"],
                            cytof.group[shared.sample, "NK_Cells"],
                            cytof.group[shared.sample, "NK_Cells"],
                            cytof.group[shared.sample, "NK_Cells"],
                            cytof.group[shared.sample, "Monocytes_classical"],
                            cytof.group[shared.sample, "Monocytes_classical"],
                            cytof.group[shared.sample, "Monocytes_classical"]),
                      gene=c(rna.log[shared.sample, "BLK"],
                             rna.log[shared.sample, "BANK1"],
                             rna.log[shared.sample, "CD79A"],
                             rna.log[shared.sample, "GZMB"],
                             rna.log[shared.sample, "GNLY"],
                             rna.log[shared.sample, "KLRC1"],
                             rna.log[shared.sample, "CLEC4E"],
                             rna.log[shared.sample, "CD14"],
                             rna.log[shared.sample, "CD163"]),
                      group=c(rep("BLK", length(shared.sample)),
                              rep("BANK1", length(shared.sample)),
                              rep("CD79A", length(shared.sample)),
                              rep("GZMB", length(shared.sample)),
                              rep("GNLY", length(shared.sample)),
                              rep("KLRC1", length(shared.sample)),
                              rep("CLEC4E", length(shared.sample)),
                              rep("CD14", length(shared.sample)),
                              rep("CD163", length(shared.sample))))

p <- ggplot(df.plot, aes(x=pop, y=gene)) +
  facet_wrap(~group, axes = "all", ncol=3, scales="free") +
  geom_point(shape=21, fill="gray78") +
  geom_smooth(method="lm",
              color = "tomato",
              fill = "tomato", 
              linetype=2) +
  stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top",
           size=3) +
  theme_classic() +
  theme(strip.background = element_blank(),
    strip.text = element_text(face = "plain"))
p

pdf("Network-examples.pdf", width=6, height=6)
print(p)
dev.off()

# pop-protein associations ####
shared.sample <- intersect(rownames(protein), rownames(cytof.group))
pop.prot <- cor(cytof.group[shared.sample,], protein[shared.sample,], method="spearman", use="complete.obs") %>%
  as.data.frame() %>% rownames_to_column("pop") %>% reshape2::melt(id.vars="pop")
pop.prot$variable <- as.character(pop.prot$variable)
pop.prot.summ <- pop.prot %>% group_by(pop) %>% arrange(desc(value)) %>% summarise(top.prot=list(head(variable, n=10)))


# ICC per pop ####
gene.unique <- cell.gene %>% group_by(gene) %>% summarise(n=n()) %>% filter(n==1) %>% pull(gene)

df.plot <- cell.gene %>% mutate(icc=icc.df$icc.icc[match(cell.gene$gene, icc.df$var)]) %>% filter(gene %in% gene.unique)
pop.ord <- df.plot %>% group_by(pop) %>% summarise(median_icc=median(icc)) %>% arrange(desc(median_icc)) %>% pull(pop)
df.plot$pop <- factor(df.plot$pop, levels=pop.ord)

#retain only pop wtih at least 20 genes
MIN_N_GENES <- 20
pop.keep <- cell.gene %>% group_by(pop) %>% summarise(n.gene=n()) %>% filter(n.gene>=MIN_N_GENES) %>% pull(pop) %>% unique()
cell.gene <- cell.gene %>% filter(pop %in% pop.keep)

top.gene <- df.plot %>% group_by(pop) %>% arrange(desc(icc)) %>% summarise(g=list(head(gene, n=5))) %>% pull(g) %>% unlist()
df.plot$label <- ifelse(df.plot$gene %in% top.gene, df.plot$gene, NA)

n.gene <- df.plot %>% group_by(pop) %>% summarise(n=n()) 
p <- ggplot(df.plot, aes(y=pop, x=icc)) +
  geom_quasirandom(aes(color=pop),alpha=0.5) +
  geom_boxplot(outlier.shape = NA, fill=NA) +
  geom_vline(xintercept=median(df.plot[!duplicated(df.plot$gene),]$icc), linetype=2) +
  geom_text(data=n.gene, aes(y=pop, x=1, label=paste0("n=", n)), hjust=0, size=4) +
  geom_text(aes(label=label), size=2, color="black") +
  scale_color_manual(values=pop.palette) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1)) +
  xlim(0,1.2) +
  ylab("") + xlab("ICC") +
  coord_flip()
p

pdf("ICC-per-pop-0711.pdf", width=8, height=3)
print(p)
dev.off()

tmp <- df.plot %>% group_by(pop) %>% arrange(desc(icc)) %>% summarise(g=list(head(gene, n=5)))
tmp$g <- as.character(tmp$g)
tmp <- setNames(tmp$g, tmp$pop)

# Enrichment of ad-hoc units ####
gene.unit <- list(
  Bcell=cell.gene %>% filter(pop %in% c("Naive_B_Cells", "B_Cells")) %>% pull(gene) %>% unique(),
  Cytotoxic=cell.gene %>% filter(pop %in% c("TEMRA_CD8", "NK_Cells")) %>% pull(gene) %>% unique(),
  CD4=cell.gene %>% filter(pop %in% c("Naive_CD4", "Central_Memory_CD4", "Effector_Memory_CD4", "TEMRA_CD4")) %>% pull(gene) %>% unique(),
  CD8=cell.gene %>% filter(pop %in% c("Naive_CD8", "Central_Memory_CD8", "Effector_Memory_CD8")) %>% pull(gene) %>% unique(),
  Myeloid=cell.gene %>% filter(pop %in% c("Myeloid_DC", "Monocytes_classical")) %>% pull(gene) %>% unique(),
  MonoNonClass=cell.gene %>% filter(pop %in% c("Monocytes_intermediate", "Monocytes_nonclassical")) %>% pull(gene) %>% unique(),
  pDC=cell.gene %>% filter(pop == "Plasmacytoid_DC") %>% pull(gene) %>% unique(),
  basophils=cell.gene %>% filter(pop == "Basophils") %>% pull(gene) %>% unique()
)

pop.unit <- list(
  Bcell = c("Naive_B_Cells", "B_Cells"),
  Cytotoxic = c("TEMRA_CD8", "NK_Cells"),
  CD4 = c("Naive_CD4", "Central_Memory_CD4", "Effector_Memory_CD4", "TEMRA_CD4"),
  CD8 = c("Naive_CD8", "Central_Memory_CD8", "Effector_Memory_CD8"),
  Myeloid = c("Myeloid_DC", "Monocytes_classical"),
  MonoNonClass = c("Monocytes_intermediate", "Monocytes_nonclassical"),
  pDC = "Plasmacytoid_DC",
  basophils = "Basophils"
)

# GSEA ####
library("BiocParallel")
df.GSEA <- data.frame()
for (u in names(gene.unit)){
  g <- gene.unit[[u]]
  
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
  
  if (any(gsea_result@result$p.adjust<0.05)){
    df.GSEA <- rbind(df.GSEA, data.frame(gsea_result@result, pop=u))
    
    p <- dotplot(gsea_result, color="NES") + ggtitle(u)
    print(p)
    # cairo_pdf(paste0(u, "-GSEA.pdf"), width=8, height=4)
    # barplot(go_result)
    # dev.off()
  }
}

str <- df.GSEA$leading_edge %>% strsplit(", ") %>% unlist()
str <- gsub("%", "", str)
str <- gsub("tags=|list=|signal=", "", str)
df.GSEA$tags <- str[seq(1,length(str),by=3)] %>% as.numeric()
df.GSEA$list <- str[seq(2,length(str),by=3)] %>% as.numeric()
df.GSEA$signal <- str[seq(3,length(str),by=3)] %>% as.numeric()

#Plot
u.keep <- c("Bcell", "Cytotoxic", "Myeloid")
df.plot <- data.frame()
for (u in u.keep){
  df.u <- df.GSEA %>% filter(pop==u) %>% arrange(desc(tags)) %>% head(n=4)
  df.plot <- rbind(df.plot, df.u)
}

df.plot$Description <- factor(df.plot$Description, levels=df.plot$Description %>% rev())
df.plot$pop <- factor(df.plot$pop, c("Myeloid","Bcell", "Cytotoxic") %>% rev())
p <- ggplot(df.plot, aes(x=tags, y=Description)) +
  facet_wrap(~pop, ncol=1, scales = "free_y") +
  geom_point(aes(color=p.adjust, size=NES)) +
  scale_color_gradient(low="#660000", high="#FF9999") +
  theme_classic() +
  theme(
    strip.background = element_blank(),  # removes facet header background
    strip.text = element_text(color = "black"),
    legend.position="bottom"
  ) +
  ylab("")
p

pdf("GSEA-unit.pdf", height=5, width=7)
print(p)
dev.off()

# GSEA of very cell type ####
df.GSEA <- data.frame()
for (pop.n in pop.keep){
  g <- cell.gene %>% filter(pop==pop.n) %>% pull(gene)
  
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
  
  if (any(gsea_result@result$p.adjust<0.05)){
    df.GSEA <- rbind(df.GSEA, data.frame(gsea_result@result, pop=pop.n))
    
    p <- dotplot(gsea_result, color="NES") + ggtitle(pop.n)
    print(p)
  }
}

str <- df.GSEA$leading_edge %>% strsplit(", ") %>% unlist()
str <- gsub("%", "", str)
str <- gsub("tags=|list=|signal=", "", str)
df.GSEA$tags <- str[seq(1,length(str),by=3)] %>% as.numeric()
df.GSEA$list <- str[seq(2,length(str),by=3)] %>% as.numeric()
df.GSEA$signal <- str[seq(3,length(str),by=3)] %>% as.numeric()


df.plot <- data.frame()
for (u in pop.keep){
  df.u <- df.GSEA %>% filter(pop==u) %>% arrange(desc(tags)) %>% head(n=4)
  df.plot <- rbind(df.plot, df.u)
}

p <- ggplot(df.plot, aes(x=tags, y=Description)) +
  facet_wrap(~pop, ncol=1, scales = "free_y", axes="all") +
  geom_point(aes(color=p.adjust, size=NES)) +
  scale_color_gradient(low="#660000", high="#FF9999") +
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(color = "black"),
    legend.position="bottom"
  ) +
  ylab("")
p

pdf("GSEA-celltype.pdf", height=10, width=7)
print(p)
dev.off()


# GO of intersect of TEMRA CD8 and NK ####
g <- intersect(cell.gene$gene[cell.gene$pop=="TEMRA_CD8"], cell.gene$gene[cell.gene$pop=="NK_Cells"])
gene_entrez <- bitr(g, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_ids <- gene_entrez$ENTREZID

go_result <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
pdf("GO_NK_TEMRACD8.pdf", width=8, height=4)
barplot(go_result)
dev.off()

my_fisher_test(cell.gene$gene[cell.gene$pop=="TEMRA_CD8"], colnames(rna), cell.gene$gene[cell.gene$pop=="NK_Cells"])

# GO of each unit ####

fisher_GO <- function(g){
  #GO
  GO.obj <- enrichGO(
    gene = g, 
    OrgDb = org.Hs.eg.db, 
    ont = "BP",
    keyType = "SYMBOL", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05,
    universe = ensg.symbol$symbol %>% unique()
  ) 
  
  if (any(GO.obj@result$p.adjust<0.05)){
    df <- GO.obj@result 
    df <- df[order(df$p.adjust, decreasing=F),] %>% subset(p.adjust<0.05)
    
    df$GeneRatio <- sapply(df$GeneRatio, function(x) {
      parts <- as.numeric(unlist(strsplit(x, "/")))
      parts[1] / parts[2] })    
    df$BgRatio <- sapply(df$BgRatio, function(x) {
      parts <- as.numeric(unlist(strsplit(x, "/")))
      parts[1] / parts[2] })
    
    df$OR <- df$GeneRatio / df$BgRatio
    
    df$Description <- factor(df$Description, levels=df$Description[order(df$OR, decreasing=F)])
    p <- ggplot(df, aes(x=OR, y=Description, size=GeneRatio, color=p.adjust)) + 
      geom_point() + 
      geom_vline(xintercept=1, linetype=2) +
      scale_color_gradientn(colors=c("goldenrod", "gray78"), 
                            breaks=c(0,0.05,1),
                            limits=c(0,1),
                            values = scales::rescale(c(0,0.05,1))) +
      theme_classic2() +
      theme(text=element_text(size=16, family="sans"),
            plot.title = element_text(hjust = 0.5)) +
      xlab("Odds ratio") + ylab("")
  } else {
    p <- ggplot()
    df <- data.frame()
  }
  
  return(list(p=p, df=df))
}

df.GO <- data.frame()
for (unit.name in names(pop.unit)){
  print(unit.name)
  g <- cell.gene %>% filter(pop %in% pop.unit[[unit.name]]) %>% pull(gene) %>% unique()
  df <- fisher_GO(g)$df
  if (nrow(df)>0){
    df.GO <- rbind(df.GO, data.frame(df, unit=unit.name))
  }
}



# Supplementary Table 3 ####

df <- cell.gene[, c("pop", "gene", "estimate", "pval", "pval.BH")]
colnames(df) <- c("Population", "Gene", "Estimate", "P-value", "Adjusted p-value")

df2 <- df.GO[,c("ID", "Description", "pvalue", "p.adjust", "qvalue", "Count", "pop")]
colnames(df2) <- c("Pathway ID", "Pathway description", "P-value", "Adjusted p-value", "Q-value", "Number of genes", "Population")

wb <- createWorkbook()
addWorksheet(wb, "PBMC association network")
addWorksheet(wb, "Gene Ontology enrichment")
writeData(wb, "PBMC association network", df)
writeData(wb, "Gene Ontology enrichment", df2)

saveWorkbook(wb, "Supplementary Table 3.xlsx", overwrite = TRUE)


# Known markers ####
# Check DICE annotation
CD4naive.DICE <- read.csv("DICE/CD4_NAIVEtop_unique_SNPs.csv")$Gene.name
CD8naive.DICE <- read.csv("DICE/CD8_NAIVEtop_unique_SNPs.csv")$Gene.name
Bnaive.DICE <- read.csv("DICE/B_CELL_NAIVEtop_unique_SNPs.csv")$Gene.name
MonoClass.DICE <- read.csv("DICE/MONOCYTEStop_unique_SNPs.csv")$Gene.name
MonoNonClass.DICE <- read.csv("DICE/M2top_unique_SNPs.csv")$Gene.name
NK.DICE <- read.csv("DICE/NKtop_unique_SNPs.csv")$Gene.name

CD4naive.DICE %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
CD8naive.DICE %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
Bnaive.DICE %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
MonoClass.DICE %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
MonoNonClass.DICE %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
NK.DICE %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }


#Check Azimuth annotation, https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC
CD4naive.Azi <- c("TCF7", "CD4", "CCR7", "IL7R", "FHIT", "LEF1", "MAL", "NOSIP", "LDHB", "PIK3IP1")
CD4TEM.Azi <- c("IL7R", "CCL5", "FYB1", "GZMK", "IL32", "GZMA", "KLRB1", "TRAC", "LTB", "AQP3")
CD4TCM.Azi <- c("IL7R", "TMSB10", "CD4", "ITGB1", "LTB", "TRAC", "AQP3", "LDHB", "IL32", "MAL")
CD8naive.Azi <- c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1")
CD8TEM.Azi <- c("CCL5", "GZMH", "CD8A", "TRAC", "KLRD1", "NKG7", "GZMK", "CST7", "CD8B", "TRGC2")
CD8TCM.Azi <- c("CD8B", "ANXA1", "CD8A", "KRT1", "LINC02446", "YBX3", "IL7R", "TRAC", "NELL2", "LDHB")
Bnaive.Azi <- c("IGHM", "IGHD", "CD79A", "IL4R", "MS4A1", "CXCR4", "BTG1", "TCL1A", "CD79B", "YBX3")
Bmemory.Azi <- c("MS4A1", "COCH", "AIM2", "BANK1", "SSPN", "CD79A", "TEX9", "RALGPS2", "TNFRSF13C", "LINC01781",	"IGHA2", "MZB1", "TNFRSF17", "DERL3", "TXNDC5", "TNFRSF13B", "POU2AF1", "CPNE5", "HRASLS2", "NT5DC2")
mDC.Azi <- c("CLEC9A", "DNASE1L3", "C1orf54", "IDO1", "CLNK", "CADM1", "FLT3", "ENPP1", "XCR1", "NDRG2", "FCER1A", "HLA-DQA1", "CLEC10A", "CD1C", "ENHO", "PLD4", "GSN", "SLC38A1", "NDRG2", "AFF3")
pDC.Azi <- c("ITM2C", "PLD4", "SERPINF1", "LILRA4", "IL3RA", "TPM2", "MZB1", "SPIB", "IRF4", "SMPD3")
NK.Azi <- c("GNLY", "TYROBP", "NKG7", "FCER1G", "GZMB", "TRDC", "PRF1", "FGFBP2", "SPON2", "KLRF1","XCL2", "FCER1G", "SPINK2", "TRDC", "KLRC1", "XCL1", "SPTSSB", "PPP1R9A", "NCAM1", "TNFRSF11A")
MonoClass.Azi <- c("S100A9", "CTSS", "S100A8", "LYZ", "VCAN", "S100A12", "IL1B", "CD14", "G0S2", "FCN1")
MonoNonClass.Azi <- c("CDKN1C", "FCGR3A", "PTPRC", "LST1", "IER5", "MS4A7", "RHOC", "IFITM3", "AIF1", "HES4")
ILC.Azi <- c("KIT", "TRDC", "TTLL10", "LINC01229", "SOX4", "KLRB1", "TNFRSF18", "TNFRSF4", "IL1R1", "HPGDS")

CD4naive.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
CD4TEM.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
CD4TCM.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
CD8naive.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
CD8TEM.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
CD8TCM.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
Bnaive.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
Bmemory.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
mDC.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
pDC.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
NK.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
MonoClass.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
MonoNonClass.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }
ILC.Azi %>% match(cell.gene$gene) %>% na.omit() %>% { cell.gene[., "pop"] }


#combined markers
CD4naive <- c(CD4naive.DICE, CD4naive.Azi) %>% unique() %>% intersect(colnames(rna))
CD4TEM <- CD4TEM.Azi %>% intersect(colnames(rna))
CD4TCM <- CD4TCM.Azi %>% intersect(colnames(rna))
CD8naive <- c(CD8naive.DICE, CD8naive.Azi) %>% unique() %>% intersect(colnames(rna))
CD8TEM <- CD8TEM.Azi %>% intersect(colnames(rna))
CD8TCM <- CD8TCM.Azi %>% intersect(colnames(rna))
Bnaive <- c(Bnaive.DICE, Bnaive.Azi) %>% unique() %>% intersect(colnames(rna))
Bmemory <- Bmemory.Azi %>% intersect(colnames(rna))
mDC <- mDC.Azi %>% intersect(colnames(rna))
pDC <- pDC.Azi %>% intersect(colnames(rna))
NK <- c(NK.DICE, NK.Azi) %>% unique() %>% intersect(colnames(rna))
MonoClass <- c(MonoClass.DICE, MonoClass.Azi) %>% unique() %>% intersect(colnames(rna))
MonoNonClass <- c(MonoNonClass.DICE, MonoNonClass.Azi) %>% unique() %>% intersect(colnames(rna))



Bcell.marker <- c(Bnaive, Bmemory) %>% unique()
Cytotoxic.marker <- NK %>% unique()
CD4.marker <- c(CD4naive, CD4TEM, CD4TCM) %>% unique()
CD8.marker <- c(CD8naive, CD8TEM, CD8TCM) %>% unique()
Myeloid.marker <- c(MonoClass, mDC) %>% unique()
MonoNonClass.marker <- MonoNonClass %>% unique()
Plasmacytoid_DC.marker <- pDC %>% unique()


#intersect with cell-gene network
intersect(cell.gene %>% filter(pop=="Naive_B_Cells") %>% pull(gene), Bnaive)
intersect(cell.gene %>% filter(pop=="B_Cells") %>% pull(gene), Bmemory)
intersect(cell.gene %>% filter(pop=="Naive_CD4") %>% pull(gene), CD4naive)
intersect(cell.gene %>% filter(pop=="Central_Memory_CD4") %>% pull(gene), CD4TCM)
intersect(cell.gene %>% filter(pop=="Effector_Memory_CD4") %>% pull(gene), CD4TEM)
intersect(cell.gene %>% filter(pop=="Naive_CD8") %>% pull(gene), CD8naive)
intersect(cell.gene %>% filter(pop=="Central_Memory_CD8") %>% pull(gene), CD8TCM)
intersect(cell.gene %>% filter(pop=="NK_Cells") %>% pull(gene), NK)
intersect(cell.gene %>% filter(pop=="Monocytes_classical") %>% pull(gene), MonoClass)
intersect(cell.gene %>% filter(pop=="Monocytes_nonclassical") %>% pull(gene), MonoNonClass)
intersect(cell.gene %>% filter(pop=="Plasmacytoid_DC") %>% pull(gene), pDC)
intersect(cell.gene %>% filter(pop=="Myeloid_DC") %>% pull(gene), mDC)


df.plot <- data.frame(gene=Bcell.marker, pop="Bcell.marker") %>% 
  rbind( data.frame(gene=Cytotoxic.marker, pop="Cytotoxic.marker") ) %>% 
  rbind( data.frame(gene=CD4.marker, pop="CD4.marker") ) %>% 
  rbind( data.frame(gene=CD8.marker, pop="CD8.marker") ) %>% 
  rbind( data.frame(gene=Myeloid.unit, pop="Myeloid.unit") ) %>% 
  rbind( data.frame(gene=MonoNonClass.unit, pop="MonoNonClass.unit") ) %>% 
  rbind( data.frame(gene=Plasmacytoid_DC.unit, pop="Plasmacytoid_DC.unit") ) 
df.plot$gene.label <- paste0(df.plot$gene, ";", df.plot$pop)


for (unit.name in names(pop.unit)){
  df.plot[,unit.name] <- df.plot$gene %in% (cell.gene %>% filter(pop %in% pop.unit[[unit.name]]) %>% pull(gene))
}


g.ord <- df.plot$gene.label
df.plot <- reshape2::melt(df.plot, id.vars=c("gene", "pop", "gene.label")) %>% filter(value==TRUE)

df.plot$gene.label <- factor(df.plot$gene.label, levels=g.ord)
p <- ggplot(df.plot, aes(y = gene, x = variable)) +
  facet_wrap( ~pop, scales = "free" ) +
  geom_point(shape=15) +
  theme_classic() +
  theme(axis.text.y = element_text(size=4), strip.background = element_blank()) +
  xlab("") + ylab("")
p

pdf("marker-genes.pdf", height=10, width=6)
print(p)
dev.off()


pop.palette <- c("monocytes"="sienna2", "CD8pos"="#E31A1C", "NK_cells"="turquoise", 
                 "CD8pos_naive"="#FB9A99", "CD4pos_naive"="#B2DF8A",             
                 "basophils"="goldenrod2", 
                 "B_cells"="dodgerblue4",               
                 "Dendritic_cells"="bisque3", "Naive_B_cells"="dodgerblue2", 
                 "CD4pos"="#33A02C", "Other"="gray78")

# HPA validation ####
source('C:/Users/albze08/Desktop/postDoc/functions/my_fisher_test.R')
enrich_HPA <- function(g){
  g.ensg <- ensg.symbol$ensg[ensg.symbol$symbol %in% g]
  
  f.Bnaive <- my_fisher_test(g.ensg, ensg.symbol$ensg %>% unique() %>% na.omit(), Bnaive$Ensembl)
  f.Bmemory <- my_fisher_test(g.ensg, ensg.symbol$ensg %>% unique() %>% na.omit(), Bmemory$Ensembl)
  f.CD4naive <- my_fisher_test(g.ensg, ensg.symbol$ensg %>% unique() %>% na.omit(), cd4_naive$Ensembl)
  f.CD8naive <- my_fisher_test(g.ensg, ensg.symbol$ensg %>% unique() %>% na.omit(), cd8_naive$Ensembl)
  f.CD4memory <- my_fisher_test(g.ensg, ensg.symbol$ensg %>% unique() %>% na.omit(), cd4_memory$Ensembl)
  f.CD8memory <- my_fisher_test(g.ensg, ensg.symbol$ensg %>% unique() %>% na.omit(), cd8_memory$Ensembl)
  f.monocyte <- my_fisher_test(g.ensg, ensg.symbol$ensg %>% unique() %>% na.omit(), monocytes$Ensembl)
  f.DC <- my_fisher_test(g.ensg, ensg.symbol$ensg %>% unique() %>% na.omit(), DC$Ensembl)
  f.basophil <- my_fisher_test(g.ensg, ensg.symbol$ensg %>% unique() %>% na.omit(), basophil$Ensembl)
  
  overlap.HPA <- data.frame(pop=c("monocytes",  "Dendritic_cells", "basophils", 
                                  "B_cells", "Naive_B_cells", "CD4pos_naive", 
                                  "CD4pos",  "CD8pos", "CD8pos_naive"),
                            OR=c(f.monocyte$odds_ratio, f.DC$odds_ratio, f.basophil$odds_ratio, 
                                 f.Bmemory$odds_ratio, f.Bnaive$odds_ratio, f.CD4naive$odds_ratio, 
                                 f.CD4memory$odds_ratio, f.CD8memory$odds_ratio, f.CD8naive$odds_ratio),
                            n=c(f.monocyte$n_genes, f.DC$n_genes, f.basophil$n_genes, 
                                f.Bmemory$n_genes, f.Bnaive$n_genes, f.CD4naive$n_genes, 
                                f.CD4memory$n_genes, f.CD8memory$n_genes, f.CD8naive$n_genes),
                            pval=c(f.monocyte$p_value, f.DC$p_value, f.basophil$p_value, 
                                   f.Bmemory$p_value, f.Bnaive$p_value, f.CD4naive$p_value, 
                                   f.CD4memory$p_value, f.CD8memory$p_value, f.CD8naive$p_value),
                            n.HPA=c(nrow(monocytes), nrow(DC), nrow(basophil), 
                                    nrow(Bmemory), nrow(Bnaive), nrow(cd4_naive), 
                                    nrow(cd4_memory), nrow(cd8_memory), nrow(cd8_naive)))
  
  overlap.HPA$adj.pval <- p.adjust(overlap.HPA$pval, method="BH")
  
  overlap.HPA$pop <- factor(overlap.HPA$pop, levels=overlap.HPA$pop[order(overlap.HPA$OR, decreasing=F)])
  overlap.HPA$GeneRatio <- overlap.HPA$n/overlap.HPA$n.HPA
  p <- ggplot(overlap.HPA, aes(x=OR, y=pop, size=GeneRatio, color=pval)) + 
    geom_point() + 
    geom_vline(xintercept=1, linetype=2) +
    scale_color_gradientn(colors=c("goldenrod", "gray78"), 
                          breaks=c(0,0.05,1),
                          limits=c(0,1),
                          values = scales::rescale(c(0,0.05,1))) +
    theme_classic2() +
    theme(text=element_text(size=14, family="sans"),
          plot.title = element_text(hjust = 0.5)) +
    xlab("Odds ratio") + ylab("")
  p
  return(list(p=p, df=overlap.HPA))
}


df.HPA.module <- data.frame()
for (m in names(unit.pop)){
  print(m)
  g <- cell.gene %>% filter(pop %in% unit.pop[[m]]) %>% pull(gene) %>% unique()
  
  #HPA cell type
  out.HPA <- enrich_HPA(g)
  p1 <- out.HPA$p + ggtitle(m)
  df.HPA.module <- rbind(df.HPA.module, data.frame(out.HPA$df, cluster=m))
  
}
df.HPA.module$adj.pval <- p.adjust(df.HPA.module$pval, method="BH")

p <- ggplot(df.HPA.module %>% filter(adj.pval<0.05), aes(y=pop,x=cluster)) +
  geom_point(aes(color=-log10(adj.pval), size=OR)) + 
  #scale_fill_gradientn(colors=c("darkred", "tomato", "gray", "gray88"), values=c(0,0.05, 0.1,1), breaks=0.05, labels="0.05") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0), legend.position="bottom") +
  theme(text=element_text(size=10, family="Arial"), plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("") + ggtitle("HPA enrichment")
p

cairo_pdf("pdfs/HPA-enrich-031025.pdf", width=4, height=6)
print(p)
dev.off()


#HBA examples ####
hba.rna <- read.table("HPA/rna_immune_cell.tsv/rna_immune_cell.tsv", sep="\t", header=T) %>% filter(Immune.cell != "total PBMC")

pop.hba.ord <- c("naive B-cell", "memory B-cell",
                 "naive CD4 T-cell", "memory CD4 T-cell", 
                 "naive CD8 T-cell", "memory CD8 T-cell",
                 "NK-cell",
                 "T-reg", "MAIT T-cell", "gdT-cell",
                 "classical monocyte", "intermediate monocyte", "non-classical monocyte",
                 "myeloid DC", "plasmacytoid DC",
                 "basophil", "eosinophil", "neutrophil")

palette.hba <- c("naive B-cell"="#1c86eeff", "memory B-cell"="#104e8bff",
                 "naive CD4 T-cell"="#b2df8aff", "memory CD4 T-cell"="#33a02cff", 
                 "naive CD8 T-cell"="#fb9a99ff", "memory CD8 T-cell"="#e31a1cff",
                 "NK-cell"="#40e0d0ff",
                 "T-reg"="gray77", "MAIT T-cell"="gray77", "gdT-cell"="gray77",
                 "classical monocyte"="#ee7942ff", "intermediate monocyte"="#ee7942ff", "non-classical monocyte"="#ee7942ff",
                 "myeloid DC"="#fffacdff", "plasmacytoid DC"="#fffacdff",
                 "basophil"="#eeb422ff", "eosinophil"="gray77", "neutrophil"="gray77")

df.plot <- hba.rna
df.plot$Immune.cell <- factor(df.plot$Immune.cell, levels=pop.hba.ord)

gene.vec <- c("CD79A", "BANK1", "BLK", "FCRLA",
              "NKG7", "GNLY", "GZMB", "PRF1", "KLRD1", "KLRC1",
              "CD14", "CLEC4E", "CLEC4D", "CD163", "S100A9")

for (g in gene.vec){
  p <- ggplot(df.plot %>% filter(Gene.name==g), aes(x=Immune.cell, y=nTPM, fill=Immune.cell)) + 
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values=palette.hba) +
    theme_classic() + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
    xlab("") + ylab(g) + theme(legend.position = "none")
  
  #pdf(paste0("pdfs/",g, "-HBA.pdf"), width=5, height=3)
  print(p)
  #dev.off()
}

ggplot(df.plot %>% filter(`Gene.name`=="BANK1"), aes(x=Immune.cell, y=nTPM)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_classic() + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  xlab("") + ylab("BANK1")

ggplot(df.plot %>% filter(`Gene.name`=="XCL2"), aes(x=Immune.cell, y=nTPM)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_classic() + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  xlab("") + ylab("XCL2")

ggplot(df.plot %>% filter(`Gene.name`=="VCAN"), aes(x=Immune.cell, y=nTPM)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_classic() + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  xlab("") + ylab("VCAN")



# Clustering of samples ####

#PCA each unit
condense_PCA <- function(df, Npca=50){
  pca_result <- prcomp(df, center = F, scale. = F)
  df_pca <- pca_result$x[, 1:Npca]
  return(df_pca)
}

da_forcluster <- data.frame()
for (u in names(gene.unit)){
  gene.u <- rna.log[,gene.unit[[u]]] %>% condense_PCA(Npca=20) %>% scale() %>% as.data.frame()
  da_forcluster <- rbind(da_forcluster, gene.u)
}

#Prepare input data
gene.main <- cell.gene  %>% pull(gene) %>% unique()
da_forcluster <- rna.log[, gene.main] %>% scale()

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
#km_527 <- igraph::cluster_louvain(g_knn_527, resolution=0.5) #kn=20 
km_527 <- igraph::cluster_fast_greedy(g_knn_527) #kn=20 

community_527 <- km_527$membership
names(community_527) <- km_527$names

primary_cluster <- data.frame(community_527)
primary_cluster$sample <- rownames(primary_cluster)
colnames(primary_cluster)[1] <- "primary_cluster"

table(primary_cluster$primary_cluster)


# CRP individual
samples.crp <- c("3920:1","3920:2","3920:3","3920:4","3920:5","3920:6")
metadata$CRP[match(samples.crp, metadata$id)]
primary_cluster$primary_cluster[match(samples.crp, primary_cluster$sample)]

df <- data.frame(crp=metadata$CRP[match(samples.crp, metadata$id)], visit=1:6)


pdf("3920-CRP.pdf", height=1.5, width=4)
ggplot(df, aes(x=visit, y=crp)) +
  geom_point() + geom_line() +
  theme_classic()
dev.off()

# Sex distribution
primary_cluster_sex <- primary_cluster
primary_cluster_sex$sex <- metadata$Gender[match(primary_cluster$sample, metadata$id)]
primary_cluster_sex %>% group_by(primary_cluster) %>% reframe(table(sex))

df.plot <- data.frame(n=c(primary_cluster_sex %>% filter(primary_cluster=="1") %>% filter(sex=="f") %>% nrow(),
                          primary_cluster_sex %>% filter(primary_cluster=="1") %>% filter(sex=="m") %>% nrow(),
                          primary_cluster_sex %>% filter(primary_cluster=="2") %>% filter(sex=="f") %>% nrow(),
                          primary_cluster_sex %>% filter(primary_cluster=="2") %>% filter(sex=="m") %>% nrow(),
                          primary_cluster_sex %>% filter(primary_cluster=="3") %>% filter(sex=="f") %>% nrow(),
                          primary_cluster_sex %>% filter(primary_cluster=="3") %>% filter(sex=="m") %>% nrow()),
                      cluster=rep(c("A", "B", "C"), each=2),
                      sex=rep(c("f", "m"), 3))
df.plot$prc <- NA
df.plot$prc[df.plot$cluster=="A"] <- 100* df.plot$n[df.plot$cluster=="A"] / sum(df.plot$n[df.plot$cluster=="A"])
df.plot$prc[df.plot$cluster=="B"] <- 100* df.plot$n[df.plot$cluster=="B"] / sum(df.plot$n[df.plot$cluster=="B"])
df.plot$prc[df.plot$cluster=="C"] <- 100* df.plot$n[df.plot$cluster=="C"] / sum(df.plot$n[df.plot$cluster=="C"])

p <- ggplot(df.plot, aes(x=cluster, y=prc, fill=sex)) + 
  geom_bar(stat="identity", position="stack") +
  geom_text(aes(label=round(prc))) +
  theme_classic()
p

pdf("sex-cluster.pdf", width=3, height=3)
print(p)
dev.off()


#Stability across visits

cluster.per.ind <- primary_cluster %>% mutate(ind=gsub("\\:.*", "", sample)) %>% 
  group_by(ind) %>% summarise(tab.cluster=list(table(primary_cluster))) 
n.freq.label <- rep(NA, nrow(cluster.per.ind))
for (n in 1:nrow(cluster.per.ind)){
  n.freq.label[n] <- max(cluster.per.ind$tab.cluster[[n]])
}

# Something about representative genes ####
df <- rna.log[,top.gene] %>% rownames_to_column(var="sample")
df$visit <- gsub(".*\\:", "", df$sample)
df$ind <- gsub("\\:.*", "", df$sample)
df <- reshape2::melt(df, id.vars=c("sample", "ind", "visit"))

ggplot(df, aes(x=visit, y=value)) +
  facet_wrap(~variable) +
  geom_line(aes(group=ind)) +
  theme_classic()

d <- left_join(primary_cluster, df, by="sample") %>% group_by(variable) %>% mutate(value=rank(value))
d$primary_cluster <- as.character(d$primary_cluster)
comparisons <- list(c(1,2), c(1,3), c(2,3))

ggplot(d, aes(x=primary_cluster, y=value)) +
  facet_wrap(~variable, nrow=4) +
  geom_boxplot() +
  stat_compare_means(method="wilcox.test", comparisons=comparisons, size=3, label="p.signif") +
  theme_classic()  +
  theme(strip.background = element_blank(),     
        panel.background = element_blank(),     # Tar bort bakgrun
        panel.border = element_blank()  ) 



#t-SNE ####
tsne_df <- Rtsne(da_forcluster, dims = 2, perplexity = 5, verbose = F, max_iter = 500)
df.plot <- tsne_df$Y[,1:2] %>% as.data.frame()
rownames(df.plot) <- rownames(da_forcluster)
colnames(df.plot)[1:2] <- c("tSNE1", "tSNE2")

df.plot$primary_cluster <- paste0("cluster", primary_cluster$primary_cluster[match(rownames(df.plot), primary_cluster$sample)] )


#Plot
p <- ggplot(df.plot, aes(x=tSNE1, y=tSNE2)) + 
  geom_point(size=2, aes(color=primary_cluster)) + 
  scale_color_manual(values=c(cluster1="#E69F00", cluster2="#56B4E9", cluster3="#009E73")) + 
  theme_classic() + theme(legend.position="none") +
  ylab("")
p

pdf("clust-tSNE-171125.pdf", height=4, width=4)
print(p)
dev.off()

table(metadata$Gender[match(primary_cluster$sample[primary_cluster$primary_cluster==1], metadata$id)])
table(metadata$Gender[match(primary_cluster$sample[primary_cluster$primary_cluster==2], metadata$id)])
table(metadata$Gender[match(primary_cluster$sample[primary_cluster$primary_cluster==3], metadata$id)])

#t-SNE, color by CD4:CD8
idx.CD4 <- grep("CD4", ignore.case = T, colnames(cytof.group))
idx.CD8 <- grep("CD8", ignore.case = T, colnames(cytof.group))
ord <- match(rownames(df.plot), rownames(cytof.group))
df.plot$CD4_CD8 <- rowSums(cytof.group[ord, idx.CD4]) / rowSums(cytof.group[ord, idx.CD8])
p0 <- ggplot(df.plot) + 
  geom_point(data=df.plot %>% filter(is.na(CD4_CD8)), size=2, aes(x=tSNE1, y=tSNE2, color=CD4_CD8)) + 
  geom_point(data=df.plot %>% filter(!is.na(CD4_CD8)), size=2, aes(x=tSNE1, y=tSNE2, color=(CD4_CD8))) + 
  scale_color_gradient(low="#F1DABF", high="#000500", na.value = "grey88") + 
  theme_classic() + theme(legend.position="bottom") +
  ylab("") + theme(legend.position = "none")
p0


#t-SNE, color by frequencies
ord <- match(rownames(df.plot), rownames(cytof.group))
df.plot$Bcells <- rowSums(cytof.group[ord, unit.pop$Bcell])
df.plot$Cytotoxic <- rowSums(cytof.group[ord, unit.pop$Cytotoxic])
df.plot$Myeloid <- rowSums(cytof.group[ord, unit.pop$Myeloid])

p1 <- ggplot(df.plot) + 
  geom_point(data=df.plot %>% filter(is.na(Bcells)), size=2, aes(x=tSNE1, y=tSNE2, color=Bcells)) + 
  geom_point(data=df.plot %>% filter(!is.na(Bcells)), size=2, aes(x=tSNE1, y=tSNE2, color=(Bcells))) + 
  scale_color_gradient(low="#F1DABF", high="#000500", na.value = "grey88") + 
  theme_classic() + theme(legend.position="bottom") +
  ylab("") + theme(legend.position = "none")
p1

p2 <- ggplot(df.plot) + 
  geom_point(data=df.plot %>% filter(is.na(Cytotoxic)), size=2, aes(x=tSNE1, y=tSNE2, color=Cytotoxic)) + 
  geom_point(data=df.plot %>% filter(!is.na(Cytotoxic)), size=2, aes(x=tSNE1, y=tSNE2, color=(Cytotoxic))) + 
  scale_color_gradient(low="#F1DABF", high="#000500", na.value = "grey88") + 
  theme_classic() + theme(legend.position="bottom") +
  ylab("") + theme(legend.position = "none")
p2


p3 <- ggplot(df.plot) + 
  geom_point(data=df.plot %>% filter(is.na(Myeloid)), size=2, aes(x=tSNE1, y=tSNE2, color=Myeloid)) + 
  geom_point(data=df.plot %>% filter(!is.na(Myeloid)), size=2, aes(x=tSNE1, y=tSNE2, color=(Myeloid))) + 
  scale_color_gradient(low="#F1DABF", high="#000500", na.value = "grey88") + 
  theme_classic() + theme(legend.position="bottom") +
  ylab("") + theme(legend.position = "none")
p3

pdf("t-SNE-sub.pdf", height=3, width=12)
ggarrange(p0,p1,p2,p3, nrow=1)
dev.off()

# Alluvial plot ####
str <- strsplit(primary_cluster$sample, ":") %>% unlist()
ind <- str[seq(1,length(str),2)]
visit <- str[seq(2,length(str),2)]

group.uniq <- unique(primary_cluster$primary_cluster)
visit.uniq <- unique(visit)

visit.list <- vector(mode="list", length(visit.uniq))
for (n in 1:length(visit.uniq)){
  visit.list[[n]] <- group.uniq
}
names(visit.list) <- paste0("visit", visit.uniq)

#all combinations across the visits
all.comb <- expand.grid(visit.list)  

#group across visits, per individual
df.ind <- data.frame(subject=unique(ind))
df.ind <- cbind(df.ind, matrix(NA, nrow(df.ind), length(visit.uniq)))
for (n in 1:nrow(df.ind)){
  for (m in 1:length(visit.uniq)){
    idx <- which(ind==df.ind$subject[n] & visit==as.character(m))
    if (length(idx)>0){
      df.ind[n,m+1] <- paste0("cluster", primary_cluster$primary_cluster[idx])
    }
  }
}
colnames(df.ind)[2:ncol(df.ind)] <- paste0("visit", 1:length(visit.uniq))

df.plot <- data.frame()
for (n in 1:6){
  x <- paste0("visit",n) %>% rep(nrow(df.ind)) #ith visit
  node <- df.ind[,n+1] #status at ith visit
  
  if (n<6){ #next visit and next status
    next_x <- paste0("visit",n+1) %>% rep(nrow(df.ind)) #ith visit
    next_node <- df.ind[,n+2]
  } else {
    next_x <- rep(NA, nrow(df.ind))
    next_node <- rep(NA, nrow(df.ind))
  }
  
  df.n <- data.frame(x=x, node=node, next_x=next_x, next_node=next_node)
  df.plot <- rbind(df.plot, df.n)
}
df.plot <- df.plot[which(!is.na(df.plot$node)),]


df.plot$x <- factor(df.plot$x, levels=paste0("visit", 1:6))
df.plot$next_x <- factor(df.plot$next_x, levels=paste0("visit", 1:6))
p <- ggplot(df.plot, aes(x = x, 
                         next_x = next_x, 
                         node = node, 
                         next_node = next_node,
                         fill = node %>% as.character(),
                         label=node)) +
  geom_alluvial(flow.alpha = .6, width=0.2) +
  scale_fill_manual(values=c(cluster1="#E69F00", cluster2="#56B4E9", cluster3="#009E73")) + 
  theme_bw() +  
  theme(text=element_text(size=12, family="Arial"), plot.title = element_text(hjust = 0.5), legend.position="none") + 
  xlab("") + ylab("Individual")
p

cairo_pdf("sankey-clustering.pdf", height=2, width=4)
print(p)
dev.off()


# Plot units' frequencies in clusters ####
unit.pop <- list(
  Bcell=c("Naive_B_Cells", "B_Cells"),
  Cytotoxic=c("NK_Cells", "TEMRA_CD8"),
  CD4=c("Naive_CD4", "Central_Memory_CD4", "Effector_Memory_CD4", "TEMRA_CD4"),
  CD8=c("Naive_CD8", "Central_Memory_CD8", "Effector_Memory_CD8"),
  Myeloid=c("Monocytes_classical", "Myeloid_DC"),
  MonoNonClass=c("Monocytes_intermediate", "Monocytes_nonclassical"),
  pDC=c("Plasmacytoid_DC"),
  basophils="Basophils"
  )

cytof.unit <- matrix(0, nrow(cytof.group), length(unit.pop))
rownames(cytof.unit) <- rownames(cytof.group)
colnames(cytof.unit) <- names(unit.pop)
for (u in names(unit.pop)){
  if (length(unit.pop[[u]])>1){
    cytof.unit[,u] <- rowSums(cytof.group[,unit.pop[[u]]])
  } else {
    cytof.unit[,u] <- cytof.group[,unit.pop[[u]]]
  }
}
cytof.unit <- data.frame(cytof.unit)

df.plot <- cytof.unit %>% mutate(across(everything(), ~ 100*rank(.) / n())) %>% 
  merge(primary_cluster, by.x="row.names", by.y="sample") %>% subset(select=-Row.names) %>%
  mutate(primary_cluster=as.character(primary_cluster)) %>%
  reshape2::melt(id.vars="primary_cluster") 

#df.plot$variable <- factor(df.plot$variable, levels=c("NK_Cells", "TEMRA_CD8", "B_Cells", "Monocytes_classical"))
comparisons <- list(c(1,2), c(1,3), c(2,3))
p.freq <- ggplot(df.plot, aes(x=primary_cluster, y=value)) +
  facet_wrap(~variable, nrow=1) +
  geom_boxplot() +
  theme(strip.text = element_blank(), 
        strip.background = element_blank()) +
  stat_compare_means(method="wilcox.test", comparisons=comparisons, size=3, label="p.signif") +
  theme_classic()  +
  theme(strip.background = element_blank(),     
        panel.background = element_blank(),     
        panel.border = element_blank()) 
p.freq

pdf("clustering-freq.pdf", height=6, width=8)
print(p)
dev.off()


# Compare CD4:CD8 ratio ####
df.plot <- primary_cluster %>% filter(sample %in% rownames(cytof))
idx.CD4 <- grep("CD4", ignore.case = T, colnames(cytof.group))
idx.CD8 <- grep("CD8", ignore.case = T, colnames(cytof.group))
ord <- match(df.plot$sample, rownames(cytof.group))
df.plot$CD4_CD8 <- rowSums(cytof.group[ord, idx.CD4]) / rowSums(cytof.group[ord, idx.CD8])

p1 <- ggplot(df.plot, aes(x=primary_cluster, y=CD4_CD8)) +
  geom_quasirandom(shape=21, aes(fill=primary_cluster)) +
  geom_boxplot(fill=NA, outlier.shape=NA) +
  scale_fill_manual(values=c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73")) + 
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")), label="p.signif") +
  scale_y_continuous(breaks = 1:9) +     # <–– adds yticks 1, 2, 3, 4, 5
  theme_classic() + theme(legend.position = "none")
p1

pdf("CD4_CD8.pdf", width=3, height=4)
print(p1)
dev.off()


# Compare innate vs adaptive ####
df.plot <- primary_cluster %>% filter(sample %in% rownames(cytof))
idx.adaptive <- grep("CD4|CD8|B_Cells", ignore.case = T, colnames(cytof.group))
idx.innate <- setdiff(1:ncol(cytof.group), idx.adaptive)
ord <- match(df.plot$sample, rownames(cytof.group))
df.plot$innate_adaptive <- rowSums(cytof.group[ord, idx.innate]) / rowSums(cytof.group[ord, idx.adaptive])

p2 <- ggplot(df.plot, aes(x=primary_cluster, y=innate_adaptive)) +
  geom_quasirandom(shape=21, aes(fill=primary_cluster)) +
  geom_boxplot(fill=NA, outlier.shape=NA) +
  scale_fill_manual(values=c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73")) + 
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")), label="p.signif") +
  theme_classic() + theme(legend.position = "none")
p2

pdf("innate_adaptive.pdf", width=3, height=4)
print(p2)
dev.off()

pdf("CD4_CD8_innate_adaptive.pdf", width=2, height=4)
ggarrange(p1,p2, ncol=1)
dev.off()

# Spider of immmune frequencies ####
library("fmsb")

#Radar of all immune frequencies
cytof.rank <- 100*apply(cytof.group,2,rank)/nrow(cytof.group)

sample.A <- primary_cluster$sample[primary_cluster$primary_cluster=="1"] %>% intersect(rownames(cytof.rank))
sample.B <- primary_cluster$sample[primary_cluster$primary_cluster=="2"] %>% intersect(rownames(cytof.rank))
sample.C <- primary_cluster$sample[primary_cluster$primary_cluster=="3"] %>% intersect(rownames(cytof.rank))

df.plot <- as.data.frame(matrix(NA, 3, ncol(cytof.rank)))
colnames(df.plot) <- colnames(cytof.rank)
for (n in 1:ncol(cytof.rank)){
  df.plot[1,n] <- cytof.rank[sample.A, n] %>% median(na.rm=T)
  df.plot[2,n] <- cytof.rank[sample.B, n] %>% median(na.rm=T)
  df.plot[3,n] <- cytof.rank[sample.C, n] %>% median(na.rm=T)
}
df.plot <- rbind(
  100,  # max
  0,  # min
  df.plot
)

#order pop
pop.ord <- c(grep("B_cells", value=T, ignore.case = T, colnames(cytof.group)),
             grep("CD4", value=T, ignore.case = T, colnames(cytof.group)),
             grep("CD8", value=T, ignore.case = T, colnames(cytof.group)),
             grep("NK", value=T, ignore.case = T, colnames(cytof.group)),
             grep("ILC", value=T, ignore.case = T, colnames(cytof.group)),
             grep("monocytes", value=T, ignore.case = T, colnames(cytof.group)),
             grep("myeloid", value=T, ignore.case = T, colnames(cytof.group)),
             grep("plasmacytoid", value=T, ignore.case = T, colnames(cytof.group)),
             grep("basophils", value=T, ignore.case = T, colnames(cytof.group)
             )
) %>% rev()
df.plot <- df.plot[,pop.ord]



pdf("radar-freq-major.pdf", height=5, width=5)
radarchart(df.plot,
           maxmin=T,
           axistype = 1,
           pcol = c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73"),      # line color
           pfcol = scales::alpha(c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73"), 0.2),  # fill color
           plty = 1,            # line width
           cglcol = "black",     # grid color
           cglty = 1,           # grid line type
           axislabcol = "#141414c9",
           cglwd = 0.8)
dev.off()


#Radar of all immune frequencies
cytof.rank <- 100*apply(cytof.nonnegative,2,rank)/nrow(cytof.nonnegative)

sample.A <- primary_cluster$sample[primary_cluster$primary_cluster=="1"] %>% intersect(rownames(cytof.rank))
sample.B <- primary_cluster$sample[primary_cluster$primary_cluster=="2"] %>% intersect(rownames(cytof.rank))
sample.C <- primary_cluster$sample[primary_cluster$primary_cluster=="3"] %>% intersect(rownames(cytof.rank))

df.plot <- as.data.frame(matrix(NA, 3, ncol(cytof.rank)))
colnames(df.plot) <- colnames(cytof.rank)
for (n in 1:ncol(cytof.rank)){
  df.plot[1,n] <- cytof.rank[sample.A, n] %>% median(na.rm=T)
  df.plot[2,n] <- cytof.rank[sample.B, n] %>% median(na.rm=T)
  df.plot[3,n] <- cytof.rank[sample.C, n] %>% median(na.rm=T)
}
df.plot <- rbind(
  100,  # max
  0,  # min
  df.plot
)

#order pop
pop.ord <- c(grep("CD4", value=T, ignore.case = T, colnames(cytof)),
             grep("CD8", value=T, ignore.case = T, colnames(cytof)),
             grep("B_cells|plasmablasts", value=T, ignore.case = T, colnames(cytof)),
             grep("CD56", value=T, ignore.case = T, colnames(cytof)),
             grep("monocytes", value=T, ignore.case = T, colnames(cytof)),
             grep("mDC", value=T, ignore.case = T, colnames(cytof)),
             grep("pDC", value=T, ignore.case = T, colnames(cytof)),
             grep("basophils", value=T, ignore.case = T, colnames(cytof)),
             grep("ILC", value=T, ignore.case = T, colnames(cytof))
)

df.plot <- df.plot[,pop.ord]
colnames(df.plot) <- paste0(1:ncol(df.plot))        

     
#replace celltype label with a number
pdf("radar-freq-all.pdf", height=5, width=5)
radarchart(df.plot,
           maxmin=T,
           axistype = 1,
           pcol = c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73"),      # line color
           pfcol = scales::alpha(c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73"), 0.2),  # fill color
           plty = 1,            # line width
           cglcol = "black",     # grid color
           cglty = 1,           # grid line type
           axislabcol = "#141414c9",
           cglwd = 0.8)
dev.off()


#Radar of immune frequencies
cytof.unit.rank <- 100*apply(cytof.unit,2,rank)/nrow(cytof.unit)

sample.A <- primary_cluster$sample[primary_cluster$primary_cluster=="1"] %>% intersect(rownames(cytof.unit.rank))
sample.B <- primary_cluster$sample[primary_cluster$primary_cluster=="2"] %>% intersect(rownames(cytof.unit.rank))
sample.C <- primary_cluster$sample[primary_cluster$primary_cluster=="3"] %>% intersect(rownames(cytof.unit.rank))
df.plot <-  data.frame(A=c(100,0,
                           cytof.unit.rank[sample.A, "Bcell"] %>% median(na.rm=T),
                           cytof.unit.rank[sample.A, "Cytotoxic"] %>% median(na.rm=T),
                           cytof.unit.rank[sample.A, "Myeloid"] %>% median(na.rm=T)),
                       B=c(100,0,
                           cytof.unit.rank[sample.B, "Bcell"] %>% median(na.rm=T),
                           cytof.unit.rank[sample.B, "Cytotoxic"] %>% median(na.rm=T),
                           cytof.unit.rank[sample.B, "Myeloid"] %>% median(na.rm=T)),
                       C=c(100,0,
                           cytof.unit.rank[sample.C, "Bcell"] %>% median(na.rm=T),
                           cytof.unit.rank[sample.C, "Cytotoxic"] %>% median(na.rm=T),
                           cytof.unit.rank[sample.C, "Myeloid"] %>% median(na.rm=T))
)

pdf("radar-freq.pdf", height=5, width=5)
radarchart(df.plot,
           maxmin=T,
           axistype = 1,
           pcol = c("#104e8bff", "#40e0d0ff", "#ee7942ff"),      # line color
           pfcol = scales::alpha(c("#104e8bff", "#40e0d0ff", "#ee7942ff"), 0.3),  # fill color
           plty = 0,            # line width
           cglcol = "black",     # grid color
           cglty = 1,           # grid line type
           axislabcol = "#141414c9",
           cglwd = 0.8)
dev.off()
           


#Radar of marker genes ####
rna.rank <- 100*apply(rna.log,2,rank)/nrow(rna.log)

sample.A <- primary_cluster$sample[primary_cluster$primary_cluster=="1"] %>% intersect(rownames(rna.rank))
sample.B <- primary_cluster$sample[primary_cluster$primary_cluster=="2"] %>% intersect(rownames(rna.rank))
sample.C <- primary_cluster$sample[primary_cluster$primary_cluster=="3"] %>% intersect(rownames(rna.rank))
df.plot <-  data.frame(A=c(100,0,
                           rna.rank[sample.A, "BLK"] %>% median(na.rm=T),
                           rna.rank[sample.A, "GZMH"] %>% median(na.rm=T),
                           rna.rank[sample.A, "CES1"] %>% median(na.rm=T)),
                       B=c(100,0,
                           rna.rank[sample.B, "BLK"] %>% median(na.rm=T),
                           rna.rank[sample.B, "GZMH"] %>% median(na.rm=T),
                           rna.rank[sample.B, "CES1"] %>% median(na.rm=T)),
                       C=c(100,0,
                           rna.rank[sample.C, "BLK"] %>% median(na.rm=T),
                           rna.rank[sample.C, "GZMH"] %>% median(na.rm=T),
                           rna.rank[sample.C, "CES1"] %>% median(na.rm=T))
)

pdf("radar-rna.pdf", height=5, width=5)
radarchart(df.plot,
           maxmin=T,
           axistype = 1,
           pcol = c("#4A8FE7", "#44E5E7", "#D89D6A"),      # line color
           pfcol = scales::alpha(c("#4A8FE7", "#44E5E7", "#D89D6A"), 0.3),  # fill color
           plty = 0,            # line width
           cglcol = "black",     # grid color
           cglty = 1,           # grid line type
           axislabcol = "#141414c9",
           cglwd = 0.8)
dev.off()


           




# Limma on RNA ####
library("limma")
primary_cluster$primary_cluster <- as.character(primary_cluster$primary_cluster)
design <- model.matrix(~0+primary_cluster, primary_cluster)

expr <- rna.log[, unique(cell.gene$gene)]
fit <- lmFit(t(expr), design)

contrast.matrix <- makeContrasts(
  cluster1 = primary_cluster1 - 0.5*(primary_cluster2 + primary_cluster3),
  cluster2 = primary_cluster2 - 0.5*(primary_cluster1 + primary_cluster3),
  cluster3 = primary_cluster3 - 0.5*(primary_cluster2 + primary_cluster1),
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)

cluster1.DEG <- topTable(fit2, coef = "cluster1", number=Inf)
cluster2.DEG <- topTable(fit2, coef = "cluster2", number=Inf)
cluster3.DEG <- topTable(fit2, coef = "cluster3", number=Inf)

#Unique unit genes
gene.unit.uniq <- gene.unit
MIN_SIZE_MODULE <- 20
for (cl in names(gene.unit.uniq)){
  for (cl2 in names(gene.unit.uniq)){
    if (cl != cl2){
      g1 <- setdiff(gene.unit[[cl]], gene.unit[[cl2]])
      g2 <- setdiff(gene.unit[[cl2]], gene.unit[[cl]])
      gene.unit.uniq[[cl]] <- g1
      gene.unit.uniq[[cl2]] <- g2
    }
  }
}

# Plot
p.list <- list()
p.list.simple <- list()
df.simple <- data.frame()
for (u in names(gene.unit)){
  g <- gene.unit.uniq[[u]]
  
  df.plot <- 
    cluster1.DEG %>% rownames_to_column("gene") %>% filter(gene %in% g) %>% filter(adj.P.Val<0.05) %>% mutate(cluster="cluster1") %>% rbind(
    cluster2.DEG %>% rownames_to_column("gene") %>% filter(gene %in% g) %>% filter(adj.P.Val<0.05) %>% mutate(cluster="cluster2")) %>% rbind(
    cluster3.DEG %>% rownames_to_column("gene") %>% filter(gene %in% g) %>% filter(adj.P.Val<0.05) %>% mutate(cluster="cluster3")) 
  
  #Divide in intervals
  NUM_BINS <- 15
  break_vec <- seq(0, max(df.plot$logFC %>% abs()), length.out=NUM_BINS)
  
  df.plot <- df.plot %>%
    mutate(interval = cut(logFC %>% abs(), breaks = break_vec)) %>% as.data.frame()
  df.plot$interval_num <- df.plot$interval %>% as.numeric()
  
  #Plot
  df.plot$group <- ifelse(df.plot$logFC>0, "positive", "negative")
  df <- df.plot %>% group_by(cluster,interval_num,group) %>% summarise(n=n())
  p.list[[u]] <- ggplot(df) +
    facet_wrap(~cluster, nrow=1) +
    geom_bar(data=df %>% filter(group=="positive"), stat="identity", position="dodge", aes(x=interval_num, y=n), fill="tomato") +
    geom_bar(data=df %>% filter(group=="negative"), stat="identity", position="dodge", aes(x=interval_num, y=-n), fill="skyblue") +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(strip.text = element_blank(), 
          strip.background = element_blank()) +
    xlim(0,NUM_BINS-1) +
    xlab("") + ylab("") + ggtitle(u)
  
  #Plot simplified
  df <- df.plot %>% group_by(cluster, group) %>% summarise(n=n()) %>% mutate(n_signed=ifelse(group=="negative", -n, n))
  p.list.simple[[u]] <- ggplot(df) +
    geom_bar(aes(x=cluster, y=n_signed, fill=group), stat="identity", position="stack") +
    geom_hline(yintercept=0, linetype=2) +
    theme_classic() + theme(legend.position="none") +
    xlab("") + ylab("") + ggtitle(u) +
    scale_fill_manual(values=c(positive="skyblue", negative="tomato"))
  
  df.simple <- rbind(df.simple, data.frame(df, pop=u))
}
p.expr <- ggarrange(plotlist = p.list.simple, nrow=1)

p <- ggarrange(p.freq, p.expr, ncol=1, heights=c(4,4))
p


p.expr <- ggplot(df.simple) +
  facet_wrap(~pop, nrow=1, scales = "free_y") +
  geom_bar(aes(x=cluster, y=n_signed, fill=group), stat="identity", position="stack") +
  geom_hline(yintercept=0, linetype=2) +
  theme_classic() + theme(legend.position="none") +
  theme(strip.background = element_blank(),     
        panel.background = element_blank(),     
        panel.border = element_blank()) +
  xlab("") + ylab("") +
  scale_fill_manual(values=c(positive="skyblue", negative="tomato"))
p.expr

pdf("expr-clust.pdf", width=12, height=2)
print(p.expr)
dev.off()

pdf("freq-clust.pdf", width=8, height=1.5)
print(p.freq)
dev.off()



# Limma on proteins ####
library("limma")
primary_cluster$primary_cluster <- as.character(primary_cluster$primary_cluster)
design <- model.matrix(~0+primary_cluster, primary_cluster)

expr <- protein[rownames(primary_cluster),]
fit <- lmFit(t(expr), design)

contrast.matrix <- makeContrasts(
  cluster1 = primary_cluster1 - 0.5*(primary_cluster2 + primary_cluster3),
  cluster2 = primary_cluster2 - 0.5*(primary_cluster1 + primary_cluster3),
  cluster3 = primary_cluster3 - 0.5*(primary_cluster2 + primary_cluster1),
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)

cluster1.DEP <- topTable(fit2, coef = "cluster1", number=Inf)
cluster2.DEP <- topTable(fit2, coef = "cluster2", number=Inf)
cluster3.DEP <- topTable(fit2, coef = "cluster3", number=Inf)




#Summarise DEPs
DEP.summ <- data.frame(upA=sum(cluster1.DEP$adj.P.Val<0.05 & cluster1.DEP$logFC>0),
                       upB=sum(cluster2.DEP$adj.P.Val<0.05 & cluster1.DEP$logFC>0),
                       upC=sum(cluster3.DEP$adj.P.Val<0.05 & cluster1.DEP$logFC>0))
DEP.summ$notSign <- nrow(cluster1.DEP) - DEP.summ$upA - DEP.summ$upB - DEP.summ$upC

# GO
p1 <- cluster1.DEP %>% rownames_to_column("protein") %>% filter(adj.P.Val<0.05) %>% 
  filter(logFC>0) %>% arrange(desc(logFC)) %>% pull(protein) 
p2 <- cluster2.DEP %>% rownames_to_column("protein") %>% filter(adj.P.Val<0.05) %>% 
  filter(logFC>0) %>% arrange(desc(logFC)) %>% pull(protein) 
p3 <- cluster3.DEP %>% rownames_to_column("protein") %>% filter(adj.P.Val<0.05) %>% 
  filter(logFC>0) %>% arrange(desc(logFC)) %>% pull(protein) 

universe_entrez <- bitr(colnames(protein), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
universe_ids <- universe_entrez$ENTREZID

go_result_1 <- go_fun(p1, universe_ids)
go_result_2 <- go_fun(p2, universe_ids)
go_result_3 <- go_fun(p3, universe_ids)

pdf("GO-DEP-A.pdf", width=6, height=3)
barplot(go_result_1) + ggtitle("Cluster A")
dev.off()

pdf("GO-DEP-B.pdf", width=6, height=3)
barplot(go_result_2) + ggtitle("Cluster B")
dev.off()

barplot(go_result_3) + ggtitle("Cluster C")


# Spider of DEPs ####
protein.rank <- 100*apply(protein,2,rank)/nrow(protein)

sample.A <- primary_cluster$sample[primary_cluster$primary_cluster=="1"] %>% intersect(rownames(protein))
sample.B <- primary_cluster$sample[primary_cluster$primary_cluster=="2"] %>% intersect(rownames(protein)) 
sample.C <- primary_cluster$sample[primary_cluster$primary_cluster=="3"] %>% intersect(rownames(protein))

df.plot <- as.data.frame(matrix(NA, 3, ncol(protein.rank)))
colnames(df.plot) <- colnames(protein.rank)
for (n in 1:ncol(protein.rank)){
  df.plot[1,n] <- protein.rank[sample.A, n] %>% median(na.rm=T)
  df.plot[2,n] <- protein.rank[sample.B, n] %>% median(na.rm=T)
  df.plot[3,n] <- protein.rank[sample.C, n] %>% median(na.rm=T)
}
df.plot <- rbind(
  100,  # max
  0,  # min
  df.plot
)


#order pop
DEP <- rbind(cluster1.DEP %>% rownames_to_column("protein") %>% mutate(cluster="A"),
             cluster2.DEP %>% rownames_to_column("protein") %>% mutate(cluster="B"),
             cluster3.DEP %>% rownames_to_column("protein") %>% mutate(cluster="C")) %>%
  filter(logFC>0) %>% arrange(adj.P.Val) %>% head(n=30) %>% arrange(cluster) %>% pull(protein) %>% rev()

df.plot <- df.plot[,DEP]


pdf("radar-protein-1711.pdf", height=5, width=5)
radarchart(df.plot,
           maxmin=T,
           axistype = 1,
           pcol = c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73"),      # line color
           pfcol = scales::alpha(c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73"), 0.2),  # fill color
           plty = 1,            # line width
           cglcol = "black",     # grid color
           cglty = 1,           # grid line type
           axislabcol = "#141414c9",
           cglwd = 0.8)
dev.off()



# Pathway enrichment ####
g1 <- cluster1.DEG %>% rownames_to_column("gene") %>% filter(adj.P.Val<0.05) %>% 
  filter(logFC>0) %>% arrange(desc(logFC)) %>% pull(gene) %>% head(n=200)
g2 <- cluster2.DEG %>% rownames_to_column("gene") %>% filter(adj.P.Val<0.05) %>% 
  filter(logFC>0) %>% arrange(desc(logFC)) %>% pull(gene) %>% head(n=200)
g3 <- cluster3.DEG %>% rownames_to_column("gene") %>% filter(adj.P.Val<0.05) %>% 
  filter(logFC>0) %>% arrange(desc(logFC)) %>% pull(gene) %>% head(n=200)

universe_entrez <- bitr(colnames(rna.log), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
universe_ids <- universe_entrez$ENTREZID

go_fun <- function(g, universe_ids){
  gene_entrez <- bitr(g, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  entrez_ids <- gene_entrez$ENTREZID
  
  go_result <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                        pvalueCutoff = 1, universe=universe_ids)
  
  return(go_result)
}

go_result_1 <- go_fun(g1, universe_ids)
go_result_2 <- go_fun(g2, universe_ids)
go_result_3 <- go_fun(g3, universe_ids)

barplot(go_result_1) + ggtitle("Cluster A")
barplot(go_result_2) + ggtitle("Cluster B")
barplot(go_result_3) + ggtitle("Cluster C")

go_tot <- rbind(
  go_result_1 %>% mutate(group="A"),
 go_result_2 %>% mutate(group="B"),
 go_result_3 %>% mutate(group="C")
 )

df.plot <- go_result_1@result %>% arrange(p.adjust) %>% head(n=8) %>% mutate(cluster="A") %>%
  rbind(go_result_2@result %>% arrange(p.adjust) %>% head(n=8) %>% mutate(cluster="B")) %>%
  rbind(go_result_3@result %>% arrange(p.adjust) %>% head(n=8) %>% mutate(cluster="C"))

df.plot$Desc.id <- paste0(df.plot$Description, "_", df.plot$cluster)
df.plot$Desc.id <- factor(df.plot$Desc.id, levels=df.plot$Desc.id[order(df.plot$Count)])
p <- ggplot(df.plot, aes(x=Count, y=Desc.id, fill=p.adjust)) + 
  facet_wrap(~cluster, scale="free_y", ncol=1) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_gradient(low="#660000", high="#FF9999") +
  theme_classic() + theme(legend.position = "bottom") + xlab("") + ylab("")
p

pdf("GO_cluster_1011.pdf", height=6, width=10)
print(p)
dev.off()


#Volcano of DEGs ####
plot_volcano <- function(df.plot){
df.plot$group <- "ns"
df.plot$group[df.plot$adj.P.Val<0.05 & df.plot$logFC>0] <-"up"
df.plot$group[df.plot$adj.P.Val<0.05 & df.plot$logFC<0] <-"down"

df.plot <- df.plot %>% arrange(desc(logFC))
df.plot$label <- NA
df.plot$label[1:30] <- rownames(df.plot)

p <- ggplot(df.plot, aes(x=logFC, y=-log10(adj.P.Val), fill=group)) +
  geom_point(shape=21) +
  geom_vline(xintercept=0, linetype=2) +
  geom_hline(yintercept=-log10(0.05), linetype=2) +
  geom_text(aes(label=label), size=2) +
  scale_fill_manual(values=c(ns="gray77", up="#606C38", down="#BC6C25")) +
  theme_classic() + theme(legend.position = "none")
return(p)
}
pdf("volcano_A.pdf", width=6, height=4)
plot_volcano(cluster1.DEG)
dev.off()

pdf("volcano_B.pdf", width=6, height=4)
plot_volcano(cluster2.DEG)
dev.off()

pdf("volcano_C.pdf", width=6, height=4)
plot_volcano(cluster3.DEG)
dev.off()


# Plot clinical phenotype in clusters ####
df.plot <- metadata %>%
  subset(select=-c(SCAPISorIGT_id, subject, Study, Visitdate, id, visit, Birthdate, subject_id,
                   Weight, Waist, Height, Hip, Bioimp_fat, Bioimp_muscle, Bioimp_bone,
                   Hb, WBC, Plt, RBC, Hct, MCV, MCH, MCHC, Neut, Lymph, Mono, Eos, Baso)) %>% 
  mutate(Gender=ifelse(Gender=="m", 1, 0)) %>%
  mutate(across(everything(), ~ 100*rank(.) / n())) %>% 
  merge(primary_cluster, by.x="row.names", by.y="sample") %>% subset(select=-Row.names) %>%
  mutate(primary_cluster=as.character(primary_cluster)) %>%
  reshape2::melt(id.vars="primary_cluster") 

comparisons <- list(c(1,2), c(1,3), c(2,3))
p <- ggplot(df.plot, aes(x=primary_cluster, y=value)) +
  facet_wrap(~variable, nrow=3) +
  geom_boxplot() +
  theme(strip.text = element_blank(), 
        strip.background = element_blank()) +
  stat_compare_means(method="wilcox.test", comparisons=comparisons, size=3, label="p.signif") +
  theme_classic()  +
  theme(strip.background = element_blank(),     
        panel.background = element_blank(),     
        panel.border = element_blank()  ) 
p

pdf("clustering-clinical.pdf", height=6, width=8)
print(p)
dev.off()


#Radar of clinical differences ####
clinical <- metadata %>%
  subset(select=-c(Gender, SCAPISorIGT_id, subject, Study, Visitdate, id, visit, Birthdate, subject_id,
                   Hb, WBC, Plt, RBC, Hct, MCV, MCH, MCHC, Neut, Lymph, Mono, Eos, Baso,
                   Housing_change, Housing_current, Address_change, MaritalStatus_change, MaritalStatus_current, ShareHousehold_change, 
                   ShareHousehold_current, Employment_change, Employment_current, Tobacco_change, Tobacco_current, PerceivedHealth,       
                   Stress, PhysicalActivity, SedentaryTime_hours, SedentaryTime_minutes, SedentaryTime_unknown, TravelAbroad,          
                   TravelAbroad_country, Animals, Animals_type,
                   Calcitriol, Calcidiol, Health_Other, Common_cold_influenzae, NSAID_painmed, Bp_med, Lipid_med,             
                   Antibiotics_med, Diab_med, Med_details, Smoking, Smoking_hours)) 

clinical.rank <- 100*apply(clinical,2,rank)/nrow(clinical)

sample.A <- primary_cluster$sample[primary_cluster$primary_cluster=="1"] 
sample.B <- primary_cluster$sample[primary_cluster$primary_cluster=="2"] 
sample.C <- primary_cluster$sample[primary_cluster$primary_cluster=="3"] 

df.plot <- as.data.frame(matrix(NA, 3, ncol(clinical.rank)))
colnames(df.plot) <- colnames(clinical.rank)
for (n in 1:ncol(clinical.rank)){
  df.plot[1,n] <- clinical.rank[sample.A, n] %>% median(na.rm=T)
  df.plot[2,n] <- clinical.rank[sample.B, n] %>% median(na.rm=T)
  df.plot[3,n] <- clinical.rank[sample.C, n] %>% median(na.rm=T)
}
df.plot <- rbind(
  100,  # max
  0,  # min
  df.plot
)


pop.ord <- c("HbA1c", "ApoB", "Chol", "LDL", "Age_at_Visit", "ALAT", "GGT", "ApoB.apoA1", "Cap_Gluc",     
             "Bioimp_bone", "Gluc", "TG", "BMI", "CRP", "TNT", "DBP",          
             "Waist", "Bioimp_muscle", "SBP", "Height", "Weight", "Crea", "Urate", "CystC",        
             "ApoA1", "HDL", "Hip", "Bioimp_fat", "ProBNP")
  
df.plot <- df.plot[,pop.ord]

pdf("radar-clinical-1711.pdf", height=5, width=5)
radarchart(df.plot,
           maxmin=T,
           axistype = 1,
           pcol = c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73"),      # line color
           pfcol = scales::alpha(c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73"), 0.2),  # fill color
           plty = 1,            # line width
           cglcol = "black",     # grid color
           cglty = 1,           # grid line type
           axislabcol = "#141414c9",
           cglwd = 0.8)
dev.off()



#lm for clinical and lifestyle
library("emmeans")

clinical <- metadata %>%
  subset(select=-c(SCAPISorIGT_id, subject, Study, Visitdate, id, visit, Birthdate, subject_id,
                   Hb, WBC, Plt, RBC, Hct, MCV, MCH, MCHC, Neut, Lymph, Mono, Eos, Baso,
                   Housing_change, Housing_current, Address_change, MaritalStatus_change, MaritalStatus_current, ShareHousehold_change, 
                   ShareHousehold_current, Employment_change, Employment_current, Tobacco_change, Tobacco_current, PerceivedHealth,       
                   Stress, PhysicalActivity, SedentaryTime_hours, SedentaryTime_minutes, SedentaryTime_unknown, TravelAbroad,          
                   TravelAbroad_country, Animals, Animals_type,
                   Calcitriol, Calcidiol, Health_Other, Common_cold_influenzae, NSAID_painmed, Bp_med, Lipid_med,             
                   Antibiotics_med, Diab_med, Med_details, Smoking, Smoking_hours)) %>% 
  mutate(Gender=ifelse(Gender=="m", 1, 0))

lm.clinical <- data.frame()
for (n in 1:ncol(clinical)){
  
  df <- primary_cluster
  df$primary_cluster <- as.character(df$primary_cluster)
  df$x <- clinical[match(df$sample, metadata$id), n]
  
  fit <- lm(x ~ primary_cluster, df)
  lm.clinical <- rbind(lm.clinical,
                       emmeans(fit, specs = "primary_cluster") %>%
                         contrast(method = "del.eff", adjust = "none") %>% as.data.frame() %>% mutate(var=colnames(clinical)[n]))
  
}
lm.clinical$adj.pval <- p.adjust(lm.clinical$p.value,  method="BH")
lm.clinical$cluster <- NA
lm.clinical$cluster[lm.clinical$contrast == "primary_cluster1 effect"] <- "clusterA"
lm.clinical$cluster[lm.clinical$contrast == "primary_cluster2 effect"] <- "clusterB"
lm.clinical$cluster[lm.clinical$contrast == "primary_cluster3 effect"] <- "clusterC"

# summary clinical differences
clinical.summ <- data.frame(upA=sum(lm.clinical$adj.pval<0.05 & lm.clinical$cluster=="clusterA" & lm.clinical$estimate>0),
                       upB=sum(lm.clinical$adj.pval<0.05 & lm.clinical$cluster=="clusterB" & lm.clinical$estimate>0),
                       upC=sum(lm.clinical$adj.pval<0.05 & lm.clinical$cluster=="clusterC" & lm.clinical$estimate>0))
clinical.summ$notSign <- ncol(clinical) - clinical.summ$upA - clinical.summ$upB - clinical.summ$upC



#Write Suppl Table with lm ####
df <- lm.clinical %>% arrange(adj.pval)
df$contrast[df$contrast=="primary_cluster1 effect"] <- "cluster A vs others"
df$contrast[df$contrast=="primary_cluster2 effect"] <- "cluster B vs others"
df$contrast[df$contrast=="primary_cluster3 effect"] <- "cluster C vs others"

wb <- createWorkbook()
addWorksheet(wb, "DPEA, cluster A vs others")
addWorksheet(wb, "DPEA, cluster B vs others")
addWorksheet(wb, "DPEA, cluster C vs others")
addWorksheet(wb, "Clinical variables")
writeData(wb, "DPEA, cluster A vs others", cluster1.DEP %>% rownames_to_column("Protein"))
writeData(wb, "DPEA, cluster B vs others", cluster2.DEP %>% rownames_to_column("Protein"))
writeData(wb, "DPEA, cluster C vs others", cluster3.DEP %>% rownames_to_column("Protein"))
writeData(wb, "Clinical variables", df)

saveWorkbook(wb, "Supplementary Table 4.xlsx", overwrite = TRUE)


