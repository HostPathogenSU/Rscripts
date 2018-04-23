###########################
#replicate reproducibility
#proteomics
###########################
#
#boxplots
#


df <- read.table("misc.txt", sep = "\t", header = T, strip.white = T)


ggplot(df, aes(df$sample, log(df$Intensity))) +
  geom_boxplot(aes( fill = as.factor(df$sample))) +
  scale_fill_brewer(palette = 'BuGn') +
  theme_bw() +
  labs(fill = "replicates") +
  ylab("Log2(LFQ Intensity)") +
  xlab("S507 replicates") +
  ggtitle("S507/S5527 inter and intragroup variation") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=90, hjust=TRUE))
  

####
#correlation heatmap
# require reshape2, ggplot, stats, hmisc, scales
####

df <- read.table("test.txt", sep = "\t", header = T, strip.white = T)

df$Protein.IDs <- NULL
df$qval <- NULL
df$LogFC <- NULL
df$av_mut <- NULL
df$avWT <- NULL
cormatrix = rcorr(as.matrix(df), type='pearson')

abbreviateSTR <- function(value, prefix){  # format string more concisely
  lst = c()
  for (item in value) {
    if (is.nan(item) || is.na(item)) { # if item is NaN return empty string
      lst <- c(lst, '')
      next
    }
    item <- round(item, 2) # round to two digits
    if (item == 0) { # if rounding results in 0 clarify
      item = '<.01'
    }
    item <- as.character(item)
    item <- sub("(^[0])+", "", item)    # remove leading 0: 0.05 -> .05
    item <- sub("(^-[0])+", "-", item)  # remove leading -0: -0.05 -> -.05
    lst <- c(lst, paste(prefix, item, sep = ""))
  }
  return(lst)
}

d <- df

cormatrix = rcorr(as.matrix(d), type='pearson')
cordata = melt(cormatrix$r)
cordata$labelr = abbreviateSTR(melt(cormatrix$r)$value, 'r')
cordata$labelP = abbreviateSTR(melt(cormatrix$P)$value, 'P')
cordata$label = paste(cordata$labelr, "\n", 
                      cordata$labelP, sep = "")
cordata$strike = ""
cordata$strike[cormatrix$P > 0.05] = "X"

txtsize <- par('din')[2] / 2

cordata$value <- round(cordata$value, digits = 3)

ggplot(cordata, aes(x=X1, y=X2, fill=value)) +
  geom_tile() + 
  scale_fill_gradientn(colours = c("blue", "cyan4", "white"), name = "Pearson", limits = c(0.5, 1)) +
  theme(axis.text.x = element_text(angle=90, hjust=TRUE)) +
  xlab("") + ylab("") + 
  geom_text(label=cordata$value, size=txtsize * 0.8, color="grey9") +
  ggtitle("Pearson correlation between CDC1551, dPPE38 and compliment CFS replicates") +
  theme(plot.title = element_text(hjust = 0.5))

#################
#pheatmap
################
require(pheatmap)
hm <- read.table("heatmap.txt", sep = "\t", header = T)
head(hm)
rnames <- hm[,1]
hm.matrix <- data.matrix(hm[,2:ncol(hm)])
rownames(hm.matrix) <- rnames
pheatmap(hm.matrix, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         cellwidth = 20, cellheight = 6.5, cluster_cols = T, scale = "row", fontsize  = 6, treeheight_col = 2, fontsize_col = 8)

pheatmap(hm.matrix, 
         treeheight_col = 0, fontsize = 6,
         color = c("Blue", "deepskyblue2","red"), 
         breaks = c(-8,-6, 0, 5))

pheatmap(hm.matrix)
############################################################################################################33
#volcanoplot
##############################################################################################################
# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)

#
#this data requires gene names, fold change and q-vals already computed
#read in csv file
results = read.csv("test.csv", header = T)

#Complicated colors motherfucker
results = mutate(results, 
                 sig = ifelse(results$padj < 0.05 & abs(results$log2FoldChange) < 1, "q-value < 0.05, FC < 1",
                              ifelse(abs(results$log2FoldChange) > 1 & results$padj > 0.05, "q-value > 0.05, FC > 1",
                                     ifelse(abs(results$log2FoldChange) > 1 & results$padj < 0.05, "significant", "non-significant"))))


#plot volcano plot

ggplot(results, aes(log2FoldChange, -log10(pvalue))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  geom_point(aes(col=sig), alpha = 0.5) +
  scale_color_manual(values=c("non-significant" = "grey9",
                              "q-value < 0.05, FC < 1" = "darkblue", 
                              "q-value > 0.05, FC > 1" = "gold4",
                              "significant" = "firebrick")) +
  xlim(c(-3,3)) +
  ylim(c(0,8)) +
  geom_label_repel(show.legend = FALSE,
                   label=ifelse(abs(results$log2FoldChange) > 1 & results$padj < 0.05, 
                                as.character(results$Gene), '')) +
  theme_bw() +
  ggtitle("mytitle") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(aes(xintercept=-1), lty = 2) +
  geom_hline(aes(yintercept= 1.30103), lty = 2) +
  geom_vline(aes(xintercept = 1), lty= 2)


####################################################################################################################################################
#Full spectral count analysis
#requires filtered dataset with normalised counts
#will calculate fold change pvals and q-vals using beta binomial with Benjamini Hochberg correction
#####################################################################################################################33
#spectral count for real NB i think the ibb should be used but whatevs we stick with what they do
#require ibb http://www.oncoproteomics.nl/software/BetaBinomial.html
#require stats
#require ggplot
#require ggrepel
#input is uniport ID, gene names, controls and tests

#read the data frame, contains uniprot[1], genenames[2], reps [3-...]
df_sc <- read.table("test.txt", header = TRUE, strip.white = TRUE, sep = "\t")

#perform the test
#bb <- bb.test(df_sc[, 3:6], colSums(df_sc[, 3:6]), c(rep("Mvu", 2), rep("dEspG1", 2)), n.threads = -1)

#modify the brackets for your own reps this data set had two reps starting at col3 and ending at col6
bb <- ibb.test(df_sc[, 3:6], colSums(df_sc[, 3:6]), c(rep("Mvu", 2), rep("dEspG1", 2)), n.threads = -1)

#calculate averages for fold change calculation [rename]
df_sc$Mvu_av <- rowMeans(na.rm = FALSE, subset(df_sc, select = c("Mvu_1", "Mvu_2")))
df_sc$dEspG1_av <- rowMeans(na.rm = FALSE, subset(df_sc, select = c("dEspG1_1", "dEspG1_2")))

#calculate the LOG2 fold change
df_sc$log.fold.change <- log2((df_sc$dEspG1_av/df_sc$Mvu_av))

# add p-vals to dataframe
df_sc$p.value <- bb$p.value

#do a benjamini hotchberg multiple comparisons correction & add to dataframe
df_sc$BH.correction <- p.adjust(df_sc$p.value, method = "BH")
df_sc$log.q.value <- c(-log10(df_sc$BH.correction))

#write to a dataframe for later use
write.table(df_sc, file = "spectralcount_MvuvdEspG1.txt", sep = "\t", row.names = FALSE)

#rework the tables in excell, give the #name a minus 1000 & inf +1000 this is presence absence dataset
#I know this is a pain, deal with it

#now read it in
df_sc <- read.table("spectralcount_MvuvdEspG1.txt", header = TRUE, sep = "\t", fill = T)

#set the colours
threshold = as.factor(abs(df_sc$log.fold.change) > 2 & df_sc$log.q.value > 1.30103)


#yay sweet sweet volcano plots
#remember to rename to fit your own dataframe, axis and titles
ggplot(data=df_sc, aes(x=df_sc$log.fold.change, y=df_sc$log.q.value, colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  scale_color_manual(values = c("grey9", "darkred"),
                     labels = c("non significant", "significant")) +
  #geom_label_repel(show.legend = FALSE, aes(label=ifelse(df_sc$log.q.value > 1.30103 & abs(df_sc$log.fold.change) > 2,as.character(df_sc$Gene.names),''))) +
  geom_text_repel(show.legend = F, aes(label = df_sc$target)) +
  theme_classic() +
  #theme(legend.position = "none") +
  xlim(c(-10,10)) + ylim(c(0,4)) +  #match the axis to your data                                                    
  labs(x ="Log2 Fold change", y = "-Log10 q-value", col = "q-value < 0.05, fold change > 2")   +
  ggtitle("Significant differential regulation Mvu vs dEspG1 in GE") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(aes(xintercept=-2), lty = 2) +
  geom_hline(aes(yintercept= 1.30103), lty = 2) +
  geom_vline(aes(xintercept = 2), lty= 2) +
  geom_vline(aes(xintercept = -10), lty= 1) +
  geom_vline(aes(xintercept = 10), lty= 1)

#presence absence plot

#separate into only presence/absence from previous dataframe
#make a new column that will contain all the spectral counts and use averages
#make the absent in sample negative for plotting in the average counts. 

df_pa <- read.table("PA.txt", strip.white = TRUE, header = TRUE, sep = "\t")

presence.absence = ifelse(df_pa$spectral.count < 0 & df_pa$BH.correction < 0.05, "absent in dEccA1", 
                          ifelse(df_pa$spectral.count > 0 & df_pa$BH.correction < 0.05, "present in dEccA1", "non significant"))

#remeber to change names here
ggplot(data = df_pa) +
  geom_point(aes(df_pa$log.q.value, df_pa$spectral.count, colour = presence.absence), alpha = 0.5) +
  scale_color_manual(values = c("non significant" = "grey9", "present in dEccA1" = "firebrick", "absent in dEccA1" = "deepskyblue2")) +
  geom_label_repel(show.legend = FALSE, 
                   aes(x = df_pa$log.q.value, y = df_pa$spectral.count,
                       label=ifelse(df_pa$log.q.value > 1.30103,
                                    as.character(df_sc$Gene.names),''))) +
  xlab("-log q-value") +
  ylab("Average spectral counts") +
  ggtitle("Presence/Absence of proteins between dEccA1 and WT in GE") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

#######################################
#Display peptides
########################################


df <- read.table("misc.txt", sep = "\t", header = T, strip.white = T)

#make the pep values nice
df$PEP <- round(-log10(df$PEP), digits = 3)

ggplot(df, aes(df$sample, df$log.Intenstity.)) +
  geom_histogram(stat = "identity", aes(fill = df$name)) + 
  scale_fill_brewer(palette = "BuGn", name = "peptides", type = "qual", direction = -1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=TRUE)) +
  ylab("Log2(Intensity)") +
  ggtitle("Peptides mapping to Rv2623_8 fusion protein") +
  theme(plot.title = element_text(hjust = 0.5))

############################################################################################
#Gene ontology vis
##########################################################################################

df <- read.table("misc.txt", sep = "\t", strip.white = T, header = T)

ggplot(df, aes(x = reorder(Term, p.value), y = -log10(df$p.value))) +
  geom_bar(stat = "identity", fill = "firebrick") +
  theme_classic() +
  guides(fill = FALSE) +
  xlab(NULL) +
  ylab("-Log10(P-Value)") +
  theme(axis.text.x = element_text(angle = 70, hjust = T, size = 12)) +
  ggtitle("Enriched GO terms in S5527 WCL") +
  theme(plot.title = element_text(hjust = 0.5))


###########################################################################################
#General t-test.
##########################################################################################
# require stats

#read in dataframe, 
#dataframe contains an identifier(gene name), observations for each rep with 3 reps and two
#groups
df <- read.table("", sep = "\t", header = T, strip.white = T)

#perform the t.test and save that bitch
#we loop this motherfucker!
# df[,2:7] is the data that should be analysed, it thus excludes col1
#t.test(x[1:3],x[4:6]) is the columns from the data extracted in df[,2:7]
#x[1:3] is the reps for the first group and x[4:6] is the reps for the second group
t.result <- apply(df[,2:7], 1, function (x) t.test(x[1:3],x[4:6],paired=TRUE))

#adding a p-val column
df$p_value <- unlist(lapply(t.result, function(x) x$p.value))

#If this is multiple testing so correct, please dear God dont forget to correct
# we are gonna use Benjamini Hochberg cause he is my boy
#As always add to the dataframe as well otherwise why bother
df$fdr <- p.adjust(df$p_value, method = "BH")

##########################################################################################
#correllogram
#ggally, hmisc, ggplot2, reshape2
#########################################################################################
df <- read.table("df7_correlelogram_timecourse.txt", sep = "\t", header = T, strip.white = T)
df$ProteinID <- NULL
df.melt <- melt(df)

library("GGally")
data(iris)
ggpairs(iris[, 1:4])


########################################################################################
#LIMMA three groups
#######################################################################################
#reads dataframe like so:
#proteinID, followed by reps, this is standard output from perseus

df <-  read.delim("PG-with-2-U-LIMMA.txt", 
                  row.names = 1, 
                  stringsAsFactors = F)[c(1:9)] # takes the first 9 columns

require(reshape2)
require(dplyr)
require(limma)
require(ggplot2)
require(ggpubr)

#make a object containing levels for later
meta.dat <- data.frame(Group = c(rep("CDC1551", 3), rep("Comp", 3), rep("dPPE38", 3)))
rownames(meta.dat) <- colnames(df)

#make df for limma
f.df <- factor(meta.dat$Group)
design <- model.matrix(~0+f.df)
colnames(design) <- levels(f.df)
#fit model
fit <- lmFit(df, design)

#pairwise comparison between three groups
#make the matrix
cont.matrix <- makeContrasts(dPPE38vWt="dPPE38-CDC1551",
                             dPPE38vComp="dPPE38-Comp",
                             CompvWt="Comp-CDC1551", levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

#this is the final results array
fit2 <- eBayes(fit2)

#use this to make venn diagram
results <- decideTests(fit2, adjust.method = "BH")

#volcano plots

#make a df containing all the neccesary shit
#1 = dPPE38vWt
#2 = dPPE38vComp
#3 = CompvWt

volc.plot.1 <- data.frame(ID = names(fit2$coefficients[,1]),
                          p = fit2$p.value[,1],
                          p.adj = p.adjust(fit2$p.value[,1], "BH"),
                          EffectSize = fit2$coefficients[,1],
                          threshold = as.factor(-log10(p.adjust(fit2$p.value[,1], "BH")) <= -log10(0.01)),
                          comparison = 1)

volc.plot.1 = mutate(volc.plot.1, 
                     sig = ifelse(volc.plot.1$p.adj < 0.05 & abs(volc.plot.1$EffectSize) < 1, "q-value < 0.05, -1 > FC < 1",
                                  ifelse(volc.plot.1$EffectSize > 1 & volc.plot.1$p.adj < 0.05, "upregulated",
                                         ifelse(volc.plot.1$EffectSize < -1 & volc.plot.1$p.adj < 0.05, "downregulated", "non-significant"))))

n.sig.down.1  = sum(volc.plot.1$p.adj < 0.05 & volc.plot.1$EffectSize < -1 ) 
n.sig.up.1  = sum(volc.plot.1$p.adj < 0.05 & volc.plot.1$EffectSize > 1 ) 

volcplotA <- ggplot(volc.plot.1, aes(x=EffectSize, y=-log10(p))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_classic() +
  geom_point(aes(col=sig), alpha = 0.5, size = 1.9) +
  scale_color_manual(values=c("non-significant" = "grey9",
                              "downregulated" = "blue", 
                              "q-value < 0.01, -1 > FC < 1" = "orange",
                              "upregulated" = "red")) +
  geom_text(aes(x = -4.9, y=10, label=n.sig.down.1)) +
  geom_text(aes(x = 4.9, y=10, label=n.sig.up.1)) +
  #ggtitle("A") +
  geom_vline(aes(xintercept=-1), lty = 5, colour = "grey") +
  geom_hline(aes(yintercept= -log10(0.05)), lty = 5, colour = "grey") +
  geom_vline(aes(xintercept = 1), lty = 5, colour = "grey") +
  theme(legend.position = "top", legend.title = element_blank()) +
  theme(title = element_text(face = "bold"))


volc.plot.2 <- data.frame(ID = names(fit2$coefficients[,2]),
                          p = fit2$p.value[,2],
                          p.adj = p.adjust(fit2$p.value[,2], "BH"),
                          EffectSize = fit2$coefficients[,2],
                          comparison = 2)

volc.plot.2 = mutate(volc.plot.2, 
                     sig = ifelse(volc.plot.2$p.adj < 0.05 & abs(volc.plot.2$EffectSize) < 1, "q-value < 0.05, -1 > FC < 1",
                                  ifelse(volc.plot.2$EffectSize > 1 & volc.plot.2$p.adj < 0.05, "upregulated",
                                         ifelse(volc.plot.2$EffectSize < -1 & volc.plot.2$p.adj < 0.05, "downregulated", "non-significant"))))

n.sig.down.2  = sum(volc.plot.2$p.adj < 0.05 & volc.plot.2$EffectSize < -1 ) 
n.sig.up.2  = sum(volc.plot.2$p.adj < 0.05 & volc.plot.2$EffectSize > 1 ) 

volcplotB <- ggplot(volc.plot.2, aes(x=EffectSize, y=-log10(p))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_classic() +
  geom_point(aes(col=sig), alpha = 0.5, size = 1.9) +
  scale_color_manual(values=c("non-significant" = "grey9",
                              "downregulated" = "blue", 
                              "q-value < 0.01, -1 > FC < 1" = "orange",
                              "upregulated" = "red")) +
  geom_text(aes(x = -4.9, y=10, label=n.sig.down.2)) +
  geom_text(aes(x = 4.9, y=10, label=n.sig.up.2)) +
  #ggtitle("B") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(aes(xintercept=-1), lty = 5, colour = "grey") +
  geom_hline(aes(yintercept= -log10(0.05)), lty = 5, colour = "grey") +
  geom_vline(aes(xintercept = 1), lty = 5, colour = "grey") +
  theme(legend.position = "top", legend.title = element_blank()) +
  theme(title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0))


volc.plot.3 <- data.frame(ID = names(fit2$coefficients[,3]),
                          p = fit2$p.value[,3],
                          p.adj = p.adjust(fit2$p.value[,3], "BH"),
                          EffectSize = fit2$coefficients[,3],
                          comparison = 3)

volc.plot.3 = mutate(volc.plot.3, 
                     sig = ifelse(abs(volc.plot.3$EffectSize) < 1 & volc.plot.3$p.adj < 0.05, "q-value < 0.01, -1 > FC < 1",
                                  ifelse(volc.plot.3$EffectSize > 1 & volc.plot.3$p.adj < 0.01, "upregulated",
                                         ifelse(volc.plot.3$EffectSize < -1 & volc.plot.3$p.adj < 0.05, "downregulated", "non-significant"))))

n.sig.down.3  = sum(volc.plot.3$p.adj < 0.05 & volc.plot.3$EffectSize < -1 ) 
n.sig.up.3  = sum(volc.plot.3$p.adj < 0.05 & volc.plot.3$EffectSize > 1 ) 
n.sig.total = sum(volc.plot.3$p.adj < 0.01 & abs(volc.plot.3$EffectSize) > 1)

volcplotC <- ggplot(volc.plot.3, aes(x=EffectSize, y=-log10(p))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_classic() +
  geom_point(aes(col=sig), alpha = 0.5, size = 1.9) +
  scale_color_manual(values=c("non-significant" = "grey9",
                              "downregulated" = "blue", 
                              "q-value < 0.01, -1 > FC < 1" = "orange",
                              "upregulated" = "red")) +
  geom_text(aes(x = -4.9, y=10, label=n.sig.down.3)) +
  geom_text(aes(x = 4.9, y=10, label=n.sig.up.3)) +
  #ggtitle("C") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(aes(xintercept=-1), lty = 5, colour = "grey") +
  geom_hline(aes(yintercept= -log10(1e-3)), lty = 5, colour = "grey") +
  geom_vline(aes(xintercept = 1), lty = 5, colour = "grey") +
  theme(legend.position = "top", legend.title = element_blank()) +
  theme(title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0))


#ggarrange(volcplotA, volcplotB, volcplotC, ncol = 2, nrow = 2, common.legend = T, legend = "top")

#write tables
write.table(volc.plot.1, file = "dPPE38vWt_Limma.txt", sep = "\t", row.names = F)
write.table(volc.plot.2, file = "dPPE38vComp_Limma.txt", sep = "\t", row.names = F)
write.table(volc.plot.3, file = "CompvWt_Limma.txt", sep = "\t", row.names = F)

############################################################################################################
#PTXQC
############################################################################################################
require(PTXQC)


txt_folder = "C:\\Users\\jgallant\\Desktop\\BAL_TMT-label\\label_check_2\\combined\\txt"


r = createReport(txt_folder)

#############################################################################################################
