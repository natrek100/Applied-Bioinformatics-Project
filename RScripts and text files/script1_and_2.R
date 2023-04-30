Sys.unsetenv("R_LIBS_USER")
dir.create("RLibrary")
.libPaths()
.libPaths(paste(getwd(), "RLibrary", sep="/"))
setRepositories()

install.packages('BiocManager')
install.packages('tidyverse')
install.packages("tximport")
install.packages('ensembldb')
install.packages('rhdf5')
install.packages('datapasta')
BiocManager::install('tximport')

library(rhdf5)
library(tidyverse)
library(tximport)
library(ensembldb) 
library(biomaRt)

targets <- read_tsv("studydesign.txt")
path <- file.path(targets$sra_accession, "abundance.tsv")
all(file.exists(path))

listMarts(host="plants.ensembl.org")
myMart <- useMart(biomart="plants_mart", host="plants.ensembl.org")
lyc.anno <- useMart(biomart="plants_mart", dataset="slycopersicum_eg_gene", host="plants.ensembl.org")
lyc.attributes <- listAttributes(lyc.anno)

Tx.lyc <- getBM(attributes=c('ensembl_transcript_id',
                             'ensembl_gene_id', 'description'),
                mart = lyc.anno)
Tx.lyc <- as_tibble(Tx.lyc)
view(Tx.lyc)
Tx.lyc <- dplyr::rename(Tx.lyc, target_id = ensembl_transcript_id,
                        gene_name = ensembl_gene_id)

Tx.lyc <- dplyr::select(Tx.lyc, "target_id", "gene_name")
view(Tx.lyc)

Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx.lyc, 
                     txOut = FALSE, #How does the result change if this =FALSE vs =TRUE?
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = FALSE)

Txi_gene

install.packages("edgeR")
install.packages('matrixStats')
install.packages('cowplot')

library(edgeR)
library(matrixStats)
library(cowplot)
install.packages(tibble)

myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
colSums(myTPM)
colSums(myCounts)

targets
sampleLabels <- targets$sample
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))
head(myTPM.stats)
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=1, size=3)

ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=1, size=2) +
  geom_smooth(method=lm) +
  geom_hex(show.legend = T, bins=20) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized data",
       caption="DIYtranscriptomics - Spring 2020") +
  # theme_classic() +
  # theme_dark() 
  theme_bw()

myDGEList <- DGEList(myCounts)
myDGEList
save(myDGEList, file = "myDGEList")
load(file = "myDGEList")
cpm <- cpm(myDGEList)
colSums(cpm)
log2.cpm <- cpm(myDGEList, log=TRUE)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
log2.cpm.df
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, 
                                  cols = FL01:FR03,
                                  names_to = "samples", 
                                  values_to = "expression") 
log2.cpm.df.pivot

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# what do you think of the distribution of this data?
# Try using coord_flip() at the end of the ggplot code

table(rowSums(myDGEList$counts==0)==10)
keepers <- rowSums(cpm>1)>=5
view(keepers)
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered)
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df,
                                           cols = FL01:FR03,
                                           names_to = "samples",
                                           values_to = "expression")
log2.cpm.filtered.df.pivot


p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = FL01:FR03,
                                                names_to = "samples",
                                                values_to = "expression")

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) + #each ggplot can be saved as an object! (p1, p2, p3)
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()


plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)
print("Step 2 complete!")
