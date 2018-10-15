# Data Import --------------------------------------------------------------------------------------------------------------------------------
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if (!require("gtools")) install.packages("gtools")
library(gtools)
if (!require("reshape2")) install.packages("reshape2")
library(reshape2)

SNP <- read_tsv("https://raw.githubusercontent.com/dkohlhase/BCB546X_R_Assignment/master/snp_position.txt")
geno <- read_tsv("https://raw.githubusercontent.com/dkohlhase/BCB546X_R_Assignment/master/fang_et_al_genotypes.txt")

# Data Inspection --------------------------------------------------------------------------------------------------------------------------------
SNP
str(SNP)
geno
str(geno)

# Data Processing --------------------------------------------------------------------------------------------------------------------------------
# Parse id, chrom loc, and nuc loc from whole SNP file
sSNP <- SNP[ ,c(1,3,4)]

  # Full Dataset --------------------------------------------------------------------------------------------------------------------------------
tgeno <- as.data.frame(t(geno))
names(tgeno) <- lapply(tgeno[1, ], as.character)
tgeno <- tgeno[-1,]
head(tgeno[,1:5])
tgeno <- rownames_to_column(tgeno, var="SNP_ID")
head(tgeno[,1:5])

mgeno <- merge(sSNP,tgeno, by.x="SNP_ID", by.y="SNP_ID", all = TRUE)
head(mgeno[,1:6])

  # Maize --------------------------------------------------------------------------------------------------------------------------------
# Parse the maize information from whole genotype file
# Maize = ZMMIL, ZMMLR, and ZMMMR
maize <- filter(geno, Group == "ZMMIL" | Group == "ZMMLR" | Group == "ZMMMR")
maize_names <- maize[,c(1:3)]
maize_names$full <- apply(maize_names, 1, paste, collapse = ", ") 

# Transpose maize files
tmaize <- as.data.frame(t(maize[,c(-1:-3)]))
colnames(tmaize) <- maize_names$full
rnames <- as.data.frame(rownames(tmaize))
rownames(tmaize) <- NULL
tmaize <- cbind(rnames,tmaize)
head(tmaize[,1:3])  #Progress check
colnames(tmaize)[1] <- "SNP_ID"

# Join SNP and geno
mmaize <- inner_join(x = sSNP, y = tmaize, by = "SNP_ID")

# Sort by increasing 'Position'
mmaize <- mmaize[mixedorder(mmaize$Position),]
#head(mmaize[,1:5], n = 10)   #Progress check

# Create new files for each 'Chromosome'
maize_chrom01_increase <- filter(mmaize, Chromosome == "1") %>% write_csv(path = "maize_chrom01_increase.csv")
maize_chrom02_increase <- filter(mmaize, Chromosome == "2") %>% write_csv(path = "maize_chrom02_increase.csv")
maize_chrom03_increase <- filter(mmaize, Chromosome == "3") %>% write_csv(path = "maize_chrom03_increase.csv")
maize_chrom04_increase <- filter(mmaize, Chromosome == "4") %>% write_csv(path = "maize_chrom04_increase.csv")
maize_chrom05_increase <- filter(mmaize, Chromosome == "5") %>% write_csv(path = "maize_chrom05_increase.csv")
maize_chrom06_increase <- filter(mmaize, Chromosome == "6") %>% write_csv(path = "maize_chrom06_increase.csv")
maize_chrom07_increase <- filter(mmaize, Chromosome == "7") %>% write_csv(path = "maize_chrom07_increase.csv")
maize_chrom08_increase <- filter(mmaize, Chromosome == "8") %>% write_csv(path = "maize_chrom08_increase.csv")
maize_chrom09_increase <- filter(mmaize, Chromosome == "9") %>% write_csv(path = "maize_chrom09_increase.csv")
maize_chrom10_increase <- filter(mmaize, Chromosome == "10") %>% write_csv(path = "maize_chrom10_increase.csv")

# Sort by decreasing 'Position'
mmaize <- mmaize[mixedorder(mmaize$Position, decreasing = T),]

# Replace ?/? with -/-
mmaize <- apply(X = mmaize, MARGIN = 2, FUN = as.character)
mmaize[mmaize == "?/?"] <- "-/-"
mmaize <- as.data.frame(mmaize)

# Create new files for each 'Chromosome'
maize_chrom01_decrease <- filter(mmaize, Chromosome == "1") %>% write_csv(path = "maize_chrom01_decrease.csv")
maize_chrom02_decrease <- filter(mmaize, Chromosome == "2") %>% write_csv(path = "maize_chrom02_decrease.csv")
maize_chrom03_decrease <- filter(mmaize, Chromosome == "3") %>% write_csv(path = "maize_chrom03_decrease.csv")
maize_chrom04_decrease <- filter(mmaize, Chromosome == "4") %>% write_csv(path = "maize_chrom04_decrease.csv")
maize_chrom05_decrease <- filter(mmaize, Chromosome == "5") %>% write_csv(path = "maize_chrom05_decrease.csv")
maize_chrom06_decrease <- filter(mmaize, Chromosome == "6") %>% write_csv(path = "maize_chrom06_decrease.csv")
maize_chrom07_decrease <- filter(mmaize, Chromosome == "7") %>% write_csv(path = "maize_chrom07_decrease.csv")
maize_chrom08_decrease <- filter(mmaize, Chromosome == "8") %>% write_csv(path = "maize_chrom08_decrease.csv")
maize_chrom09_decrease <- filter(mmaize, Chromosome == "9") %>% write_csv(path = "maize_chrom09_decrease.csv")
maize_chrom10_decrease <- filter(mmaize, Chromosome == "10") %>% write_csv(path = "maize_chrom10_decrease.csv")

  # Teosinte --------------------------------------------------------------------------------------------------------------------------------
# Parse the teosinte information from whole genotype file
# Teosinte = ZMPBA, ZMPIL, and ZMPJA

teosinte <- filter(geno, Group == "ZMPBA" | Group == "ZMPIL" | Group == "ZMPJA")
teosinte_names <- teosinte[,c(1:3)]
teosinte_names$full <- apply(teosinte_names, 1, paste, collapse = ", ") 

# Transpose teosinte files
tteosinte <- as.data.frame(t(teosinte[,c(-1:-3)]))
colnames(tteosinte) <- teosinte_names$full
rnames <- as.data.frame(rownames(tteosinte))
rownames(tteosinte) <- NULL
tteosinte <- cbind(rnames,tteosinte)
#head(tteosinte[,1:3])   #Progress check
colnames(tteosinte)[1] <- "SNP_ID"

# Join SNP and geno
mteosinte <- inner_join(x = sSNP, y = tteosinte, by = "SNP_ID")

# Sort by increasing 'Position'
mteosinte <- mteosinte[mixedorder(mteosinte$Position),]
#head(mteosinte[,1:5], n = 10)  #Progress check

# Create new files for each 'Chromosome'
teosinte_chrom01_increase <- filter(mteosinte, Chromosome == "1") %>% write_csv(path = "teosinte_chrom01_increase.csv")
teosinte_chrom02_increase <- filter(mteosinte, Chromosome == "2") %>% write_csv(path = "teosinte_chrom02_increase.csv")
teosinte_chrom03_increase <- filter(mteosinte, Chromosome == "3") %>% write_csv(path = "teosinte_chrom03_increase.csv")
teosinte_chrom04_increase <- filter(mteosinte, Chromosome == "4") %>% write_csv(path = "teosinte_chrom04_increase.csv")
teosinte_chrom05_increase <- filter(mteosinte, Chromosome == "5") %>% write_csv(path = "teosinte_chrom05_increase.csv")
teosinte_chrom06_increase <- filter(mteosinte, Chromosome == "6") %>% write_csv(path = "teosinte_chrom06_increase.csv")
teosinte_chrom07_increase <- filter(mteosinte, Chromosome == "7") %>% write_csv(path = "teosinte_chrom07_increase.csv")
teosinte_chrom08_increase <- filter(mteosinte, Chromosome == "8") %>% write_csv(path = "teosinte_chrom08_increase.csv")
teosinte_chrom09_increase <- filter(mteosinte, Chromosome == "9") %>% write_csv(path = "teosinte_chrom09_increase.csv")
teosinte_chrom10_increase <- filter(mteosinte, Chromosome == "10") %>% write_csv(path = "teosinte_chrom10_increase.csv")

# Sort by decreasing 'Position'
mteosinte <- mteosinte[mixedorder(mteosinte$Position, decreasing = T),]

# Replace ?/? with -/-
mteosinte <- apply(X = mteosinte, MARGIN = 2, FUN = as.character)
mteosinte[mteosinte == "?/?"] <- "-/-"
mteosinte <- as.data.frame(mteosinte)

# Create new files for each 'Chromosome'
teosinte_chrom01_decrease <- filter(mteosinte, Chromosome == "1") %>% write_csv(path = "teosinte_chrom01_decrease.csv")
teosinte_chrom02_decrease <- filter(mteosinte, Chromosome == "2") %>% write_csv(path = "teosinte_chrom02_decrease.csv")
teosinte_chrom03_decrease <- filter(mteosinte, Chromosome == "3") %>% write_csv(path = "teosinte_chrom03_decrease.csv")
teosinte_chrom04_decrease <- filter(mteosinte, Chromosome == "4") %>% write_csv(path = "teosinte_chrom04_decrease.csv")
teosinte_chrom05_decrease <- filter(mteosinte, Chromosome == "5") %>% write_csv(path = "teosinte_chrom05_decrease.csv")
teosinte_chrom06_decrease <- filter(mteosinte, Chromosome == "6") %>% write_csv(path = "teosinte_chrom06_decrease.csv")
teosinte_chrom07_decrease <- filter(mteosinte, Chromosome == "7") %>% write_csv(path = "teosinte_chrom07_decrease.csv")
teosinte_chrom08_decrease <- filter(mteosinte, Chromosome == "8") %>% write_csv(path = "teosinte_chrom08_decrease.csv")
teosinte_chrom09_decrease <- filter(mteosinte, Chromosome == "9") %>% write_csv(path = "teosinte_chrom09_decrease.csv")
teosinte_chrom10_decrease <- filter(mteosinte, Chromosome == "10") %>% write_csv(path = "teosinte_chrom10_decrease.csv")

  # Zygosity --------------------------------------------------------------------------------------------------------------------------------
# Create a new column indicating whether SNP is homozygous
# Melt the original genotype file
zygotes_long <- filter(geno, Group == "ZMMIL" | Group == "ZMMLR" | Group == "ZMMMR" | Group == "ZMPBA" | Group == "ZMPIL" | Group == "ZMPJA")
zygotes <- melt(zygotes_long, measure.vars = colnames(geno)[4:986])
colnames(zygotes)[4:5] <- c("SNP_ID", "Homozygous")
colnames(zygotes)
# Change all homozygous SNPs to TRUE
zygotes[zygotes == "A/A"] <- TRUE
zygotes[zygotes == "C/C"] <- TRUE
zygotes[zygotes == "G/G"] <- TRUE
zygotes[zygotes == "T/T"] <- TRUE
# Change all heterozygotes to FALSE
zygotes[zygotes == "A/C"] <- FALSE
zygotes[zygotes == "A/G"] <- FALSE
zygotes[zygotes == "A/T"] <- FALSE
zygotes[zygotes == "C/G"] <- FALSE
zygotes[zygotes == "C/T"] <- FALSE
zygotes[zygotes == "G/T"] <- FALSE
# Change all missing values to NA
zygotes[zygotes == "?/?"] <- NA
# Sort the dataframe using "Group" and "Species_ID"
zygotes <- arrange(zygotes, Sample_ID, Group)

# Data Visualization --------------------------------------------------------------------------------------------------------------------------------
# SNPs per chromosome
ggplot(data = mgeno) +
  geom_bar(mapping = aes(x = Chromosome)) +
  scale_x_discrete(limits=c(1:10, "unknown", "multiple")) +
  ggtitle(label = "SNPs per chromosome") +
  xlab(label = "Chromosome") +
  ylab(label = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
    ) 
ggsave(filename = "SNPs per chromosome.png", device = "png")

# SNPs per Group
ggplot(data = geno) + 
  geom_bar(mapping = aes(x = Group)) +
  ggtitle(label = "SNPs per Group") +
  ylab("Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 12)
    )
ggsave(filename = "SNPs per Group.png", device = "png")

# Total SNP Zygosity
ggplot(data = zygotes) +
  geom_bar(mapping = aes(x = Homozygous, fill = Homozygous), stat = "count") +
  ggtitle(label = "Total SNP Zygosity") +
  ylab(label = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none"
    )
ggsave(filename = "Total SNP Zygosity.png", device = "png")

# SNP Zygosity by Sample_ID
ggplot(data = zygotes) +
  geom_bar(mapping = aes(x = Sample_ID, fill = Homozygous), stat = "count") +
  ggtitle(label = "SNP Zygosity by Ordered Sample_ID") +
  ylab(label = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
    )
ggsave(filename = "SNP Zygosity by Ordered Sample_ID.png", device = "png")

# SNP Zygosity by Group
ggplot(data = zygotes) +
  geom_bar(mapping = aes(x = Group, fill = Homozygous), stat = "count") +
  ggtitle(label = "SNP Zygosity by Group") +
  xlab(label = "Chromosome") +
  ylab(label = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
    ) 
ggsave(filename = "SNP Zygosity by Group.png", device = "png")

# Position Adjustment for SNP Zygosity by Group
ggplot(data = zygotes) + 
  geom_bar(mapping = aes(x = Group, fill = Homozygous), position = "fill") +
  ggtitle(label = "Position Adjustment for SNP Zygosity by Group") +
  xlab(label = "Chromosome") +
  ylab(label = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
    )
ggsave(filename = "Position Adjustment for SNP Zygosity by Group.png", device = "png")

# Additional visualization
# Distribution of SNPs across chromosomes
df <- mgeno %>% 
    mutate(Dist_Bin = cut(as.numeric(Position), breaks = 20))

df1 <- subset(df, Chromosome != "unknown")
df2 <- subset(df1,  Chromosome != "NA")
df3 <- subset(df2, Chromosome != "multiple")

ggplot(data = df3) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 315, hjust = -0.02, size = 5)
  ) +
  facet_wrap(~Chromosome, scales = "free", nrow = 2) +
  geom_bar(mapping = aes(x = Dist_Bin), stat = "count") +
  ggtitle(label = "Distribution of SNPs across chromosomes") +
  xlab(label = "Position Binned") +
  ylab(label = "Number of SNPs")
ggsave(filename = "Distribution of SNPs across chromosomes.png", device = "png")