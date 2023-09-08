#BiocManager::install("GEOquery")
library(dplyr)
library(tidyverse)
library(GEOquery)

#Read in data
geo_data <- read.csv("/Users/maitreepatel/Desktop/R/GSE183947_fpkm.csv")
#checking dimensions of dataframe
dim(geo_data)
geo_data[1:1,]
#getting metadata or clinical data
#rows are genes
#columns are sample (tumour and normal) 
#series file on NCBI can be obtained using GEOquery package
gse <- getGEO(GEO = "GSE183947",
              GSEMatrix = TRUE)

#getting phenodata of the first element in the gse list
metadata <- pData(phenoData(gse[[1]]))
#metadata dataframe has 40 columns
#we only require few
#and performing additional operations

metadata.modified <- metadata %>%
  select(c(1, 10, 11, 17)) %>%
  #renaming column characteristic_ch1 to tissue
  rename(tissue = characteristics_ch1) %>% 
  #renaming column characteristic_ch1.1 to metastasis
  rename(metastasis = characteristics_ch1.1) %>%
  #global substitution of "tissue: " with "" (removing it)
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

#looking at gene expression data
#head(geo_data)
#first column corresponds to gene
#converting to long format
#reshaping data
data.long <- geo_data %>%
  rename(gene = X) %>%
  #converting wide to long format
  gather(key = "samples", 
         value = "FPKM",
         -gene) #removing gene column
head(data.long)

#adding or matching metadata to gene expression data
#join using sample as ID
data.long <- data.long %>%
  #. to input data.long
  left_join(.,
            metadata.modified,
            by = c("samples" = "description"))

#exploring data
#extracting expression of two genes  
#compare expression in tumour versus normal
data.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue)  %>%
  summarise(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(mean_FPKM)

#Visualizing expression data
library(ggplot2)

#barplot
data.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., 
         aes(x = samples,
             y = FPKM,
             fill = tissue)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90))

#density plot
data.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(aes(x = FPKM,
             fill = tissue)) +
  geom_density(alpha = 0.3)

#box plot
data.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(aes(x = metastasis,
             y = FPKM)) +
  geom_boxplot()
#need to apply statistical test to test significance in the difference
#or use of violin plot
data.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(aes(x = metastasis,
             y = FPKM)) +
  geom_violin()

#scatterplot
data.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(aes(x = BRCA1,
             y = BRCA2)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = FALSE)
#perform a correlation test to test for significance

data.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(aes(x = BRCA1,
             y = BRCA2,
             color = tissue)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = FALSE)
#we see different results in groups

#heatmap
#used to compare expresion of multiple genes

genes.of.interest <- c("BRCA1", "BRCA2", "TP53", "ALK", "MYCN")

p <- data.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(aes(x = samples,
             y = gene,
             fill = FPKM)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient(low = "white",
                      high = "red")

#saving plots
ggsave(p, filename = "heatmap.pdf", width = 10, height = 8)