install.packages("reshape2")
install.packages('heatmaply')


library(tidyverse) 
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(reshape2)
library(heatmaply)

group <- factor(targets$group)
group
design <- model.matrix(~0 + group)
design
colnames(design) <- levels(group)
colnames

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)

contrast.matrix <- makeContrasts(blooming = flower - fruit,
                                 levels=design)
contrast.matrix
fits <- contrasts.fit(fit, contrast.matrix)
fits
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")  
myTopHits.df

vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
head(results)
#myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
#dim(myTopHits)
#head(myTopHits)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)
nrow(diffGenes.df)

write_tsv(diffGenes.df,"DiffGenes.txt")
heatmaply(diffGenes.df[2:11], 
          xlab = "Samples", ylab = "DEGs", 
          main = "DEGs in cutaneous leishmaniasis",
          scale = "column",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.0000001,
          titleX = T,
          titleY = T,
          hide_colorbar = TRUE,
          branches_lwd = 0.1,
          label_names = c("Gene", "Sample:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          labCol = colnames(diffGenes.df)[2:11],
          labRow = diffGenes.df$geneID,
          heatmap_layers = theme(axis.line=element_blank())
)


