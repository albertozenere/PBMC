pack_R <- c("dplyr","ggplot2","ggrepel","umap", "edgeR",
            "RColorBrewer", "pheatmap", "tidyverse",
            "igraph", "biomaRt", "ggbeeswarm",
            "org.Hs.eg.db", "ggpubr", 
            "clusterProfiler",
            "igraph", "readxl")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}
set.seed(1)

# Load metadata ####
metadata <- read.csv("data/S3_Wellness_visit1_6_Clin_200109.txt", sep="\t", header=T)
metadata$visit <- gsub("Visit ", "", metadata$VisitName)
metadata$VisitName <- NULL
metadata$subject_id <- gsub("1-", "", metadata$subject)
metadata$id <- paste0(metadata$subject_id, ":", metadata$visit)
rownames(metadata) <- metadata$id


# Plot some stats ####
ind.uniq <- unique(metadata$subject_id)
color.UMAP <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] %>% sample(length(ind.uniq))
color.UMAP <- setNames(color.UMAP, ind.uniq)

df.plot.sex <- table(metadata$Gender) %>% as.data.frame() %>% 
  mutate(stack="yes") %>% mutate(prc=round(100*Freq/sum(Freq), digits=2)) %>%
  mutate(y=c(400,200))

p1 <- ggplot(df.plot.sex, aes(x=stack, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack") +
  geom_text(aes(label=prc, y=y)) +
  theme_classic() + theme(legend.position="none")
p1

p2 <- ggplot(metadata %>% filter(visit==1), aes(x=1, y=Age_at_Visit, fill=Gender)) +
  geom_quasirandom(size=2, shape=21, width=0.2) +
  geom_boxplot(fill=NA, outlier.shape=NA, width=0.3) +
  scale_fill_manual(values=c(m="#2daaaeff", f="#c93e34ff")) +
  theme_classic() + theme(legend.position="none")
p2

pdf("age-baseline.pdf", height=3, width=2)
print(p2)
dev.off()

ind.ord <- metadata %>% group_by(subject_id) %>% summarise(mean_BMI=mean(BMI), id_vec=list(id)) %>% arrange(mean_BMI) %>% pull(subject_id) 
df.plot.bmi <- metadata
df.plot.bmi$subject_id <- factor(df.plot.bmi$subject_id, levels=ind.ord)

p3 <- ggplot(df.plot.bmi, aes(x=subject_id, y=BMI %>% log(), fill=subject_id)) +
  geom_quasirandom(size=2, shape=21) +
  scale_fill_manual(values=color.UMAP) +
  theme_classic() + theme(legend.position="none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  xlab("samples")
p3

pdf("BMI-visits.pdf", height=3, width=9)
print(p3)
dev.off()

ind.ord <- metadata %>% group_by(subject_id) %>% summarise(mean_crp=mean(CRP), id_vec=list(id)) %>% arrange(mean_crp) %>% pull(subject_id) 
df.plot.crp <- metadata
df.plot.crp$subject_id <- factor(df.plot.crp$subject_id, levels=ind.ord)

p3 <- ggplot(df.plot.crp, aes(x=subject_id, y=CRP, fill=subject_id)) +
  geom_quasirandom(size=2, shape=21) +
  scale_y_log10() +
  scale_fill_manual(values=color.UMAP) +
  theme_classic() + theme(legend.position="none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  xlab("samples")
p3

pdf("CRP-visits.pdf", height=3, width=6)
print(p3)
dev.off()


p4 <- ggplot(metadata, aes(x=1, y=HDL)) +
  geom_quasirandom(size=2, alpha=.5) +
  geom_boxplot(fill=NA, outlier.shape=NA) +
  theme_classic()
p4

p5 <- ggplot(metadata, aes(x=1, y=SBP)) +
  geom_quasirandom(size=2, alpha=.5) +
  geom_boxplot(fill=NA, outlier.shape=NA) +
  theme_classic()
p5

p <- ggarrange(p1,p2,p3,p4,p5, nrow=1)
p

pdf("data-overview.pdf", height=3, width=10)
print(p)
dev.off()
  
# Ridge of frequencies ####
library("ggridges")
library("scales")
pop.palette <- c(NK_Cells="turquoise",           
                 Monocytes_classical="sienna2", Monocytes_intermediate="sienna2", Monocytes_nonclassical="sienna2",
                 Plasmacytoid_DC="bisque3", Myeloid_DC="bisque3",
                 Basophils="goldenrod2",             
                 Naive_B_Cells="dodgerblue2",
                 B_Cells="dodgerblue4",
                 Naive_CD4="#B2DF8A",  
                 Central_Memory_CD4="#33A02C", Effector_Memory_CD4="#33A02C", TEMRA_CD4="#33A02C",
                 Naive_CD8="#FB9A99",    
                 Central_Memory_CD8="#E31A1C", Effector_Memory_CD8="#E31A1C", TEMRA_CD8="#E31A1C",
                 ILCs="gray78")

df.plot <- cytof.group%>% reshape2::melt()

pop.ord <- df.plot %>% group_by(variable) %>% summarise(median_freq=median(value)) %>% 
  arrange(median_freq) %>% pull(variable) %>% as.character()
df.plot$variable <- factor(df.plot$variable, levels=pop.ord)
p <- ggplot(df.plot %>% filter(value > 1e-5),
            aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_manual(values = pop.palette) +
  scale_x_continuous(
    trans = 'log10',
    breaks = as.numeric(outer(1:9, 10^(-5:0), `*`)),  
    labels = function(x) {
      labs <- rep("", length(x))
      powers <- 10^(-5:0)
      labs[x %in% powers] <- parse(text = paste0("10^", log10(powers[x %in% powers])))
      labs
    }
  ) +
  annotation_logticks(sides = "b") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Log10 frequency") + ylab("")
p

pdf("freq-major.pdf", width=5, height = 6)
print(p)
dev.off()
