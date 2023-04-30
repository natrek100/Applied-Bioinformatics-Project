install.packages('DT')
install.packages('plotly')
install.packages('gt')

library(DT) 
library(plotly)
library(gt)

targets
group <- targets$group
group <- factor(group)
group
log2.cpm.filtered.norm.df

distance <- dist(t(log2.cpm.filtered.norm), method = "maximum")
clusters <- hclust(distance, method = "complete")
plot(clusters, labels=sampleLabels)

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pca.res
ls(pca.res)
summary(pca.res)
pca.res$rotation 
pca.res$x
screeplot(pca.res)
pc.var<-pca.res$sdev^2 
pc.per<-round(pc.var/sum(pc.var)*100, 1) 
pc.per

pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels) + # color = group
  geom_point(size=4) +
  # geom_label() +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

pca.res.df <- pca.res$x[,1:4] %>%
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)
view(pca.res.df)

pca.pivot <- pivot_longer(pca.res.df, 
                          cols = PC1:PC4, 
                          names_to = "PC", 
                          values_to = "loadings")
pca.pivot
ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

mydata.df <- log2.cpm.filtered.norm.df %>% 
  mutate(flower.AVG = (FL01 + FL02 + FL03 )/3,
         fruit.AVG = (FR01 + FR02 + FR03)/3,
         #now make columns comparing each of the averages above that you're interested in
         LogFC = (fruit.AVG - flower.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

mydata.df

datatable(mydata.df[,c(1:10)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         #dom = "Blfrtip", 
                         #buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))
