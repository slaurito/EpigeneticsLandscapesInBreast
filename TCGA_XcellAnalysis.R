
###########################################################################

##Complete bioinformatic analysis corresponding to the manuscript 
#"Left-Right Asymmetry in the Breast TME is associated with distinct 
#Bioelectric and Epigenetic Tumor Profiles" by Real et al. 

#Load Libraries
library(TCGAbiolinks)
library(dplyr)
library(EDASeq)
library('biomaRt')
library(ggpubr)
library(ggplot2)
library(dplyr)
library(scales)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(stringr)
library (xCell)
library(ggfortify)
library(reshape2)
library(limma)
library(edgeR)
library(GSVA)
library(tidyr)

#Download and Prepare Clinical Data
Subtypes<-TCGAquery_subtype(tumor="BRCA")

Query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab")

GDCdownload(Query)
ClinicAll <- GDCprepare(Query)  
Clinical <- ClinicAll$clinical_patient_brca

#Select IDCs only
ClinicalIDC <- Clinical[Clinical$histological_type =="Infiltrating Ductal Carcinoma",]


#Download and prepare RNAseq Data

QueryRNA <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  sample.type= c("Primary Tumor"))
GDCdownload(QueryRNA)
RnaData<-GDCprepare(QueryRNA )

#Outliers Analysis
CorOutliers <- TCGAanalyze_Preprocessing(RnaData)


# normalization of genes
DataNorm <- TCGAanalyze_Normalization(
  tabDF =CorOutliers , 
  geneInfo =  geneInfoHT,
  method = 'geneLength'
)

#Tumor purity filter
CpeTable <- Tumor.purity

CpeTable <- CpeTable[CpeTable$Cancer.type=="BRCA",]

ColsNum <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")


CpeTable[ColsNum] <- lapply(CpeTable[ColsNum], 
                            function(col) as.numeric(gsub(",",".", col)))
#Remove NAs
CpeTable <- CpeTable[!is.na(CpeTable$CPE),]

#Select tumors with purity below 70%
CpeTableFilt <- CpeTable[CpeTable$CPE <= 0.7,]

CpeTableFilt$Sample.ID <- as.character(CpeTableFilt$Sample.ID)

CpeTableFilt$barcode <-substr(CpeTableFilt$Sample.ID,1, nchar(CpeTableFilt$Sample.ID)-4)

#Match samples in purity table with samples in IDCs table
IdcCpeFilt<- ClinicalIDC[ClinicalIDC$bcr_patient_barcode %in% CpeTableFilt$barcode,]

# Remove Male samples
IdcCpeFilt <- IdcCpeFilt[IdcCpeFilt$gender =="FEMALE",]


##############   NOT   RUN ##############

#Another option for selection by tumor purity#

TumorFiltUp<-TCGAtumor_purity(colnames(DataNorm),
                            absolute=0,
                            ihc=0,
                            lump=0,
                            estimate=0,
                            cpe  = 0.7)

TumorFiltUp <- data.frame(TumorFiltUp$pure_barcodes)
TumorFiltUp$patient <- NA

colnames(TumorFiltUp) <- c("ID", "patient")

TumorFiltUp$patient<-substr(TumorFiltUp$ID,1, nchar(TumorFiltUp$ID)-12)

TumorFiltUp$barcode<-substr(TumorFiltUp$ID,1, nchar(TumorFiltUp$ID)-16)

FinalFiltered<- ClinicalIDC[!(ClinicalIDC$bcr_patient_barcode %in% TumorFiltUp$barcode),]

###############################################################################



#Match HGNC symbol/Ensembl code
EnsList <- as.data.frame(rownames(DataNorm))

Mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

GList <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                         "hgnc_symbol"),
               values=EnsList,mart= Mart)


#Data Formatting

RawData <- merge(DataNorm, GList, by.x='row.names', by.y=1)
RawData <- RawData[,-c(1)]

RawData <-  aggregate(.~hgnc_symbol, RawData, mean)
rownames(RawData) <- RawData$hgnc_symbol
RawData <- RawData[,-c(1)]

#Assign Left or Right side for all samples
BreastSite <- dplyr::select(ClinicalIDC,c(2,91))
rownames(BreastSite) <- BreastSite$bcr_patient_barcode
BreastSite$SiteDico<-NA


# Convert Anatomical Tumor Site to dicotomized variable (Left/Right) 


BreastSite$SiteDico <- str_extract(BreastSite$anatomic_neoplasm_subdivision,
                                   "(\\w+)")



Samples <- data.frame(ID= colnames(RawData))
Samples$PatID <- substr(Samples$ID, 1, nchar(Samples$ID) -16)
BreastSite <- merge (Samples, BreastSite, by.x=2, by.y=1)

BreastSite <- BreastSite[!duplicated(BreastSite$PatID),]

colnames(BreastSite) <- c("PatID", "ID", "AnatomicSubdivision", "SiteDico")

FilteredData <- RawData[,colnames(RawData) %in% BreastSite$ID]


###############################################################################

#xCell data

###############################################################################


xCellData <- xCellAnalysis(FilteredData)

xCellDataT <- data.frame(t(xCellData))

################################################################################

#####Assign ImmuneScore to Left/Right samples

################################################################################

LeftSamples <- BreastSite$ID[BreastSite$SiteDico =="Left"] 
RightSamples <- BreastSite$ID[BreastSite$SiteDico == "Right"]

ImmunoLeft <- xCellDataT[rownames(xCellDataT) %in% LeftSamples,]
ImmunoRight <- xCellDataT[rownames(xCellDataT) %in% RightSamples,]

ImmunoMeanLeft <- mean(ImmunoLeft$ImmuneScore)
ImmunoMeanRight <- mean(ImmunoRight$ImmuneScore)

###############################################################################

### BOXPLOT

################################################################################


###Data Formatting for graphical display


AllScores <- xCellDataT[, 65:67]

AllScores$Site <- NA
AllScores$Site[rownames(AllScores) %in% LeftSamples] <- "Left"
AllScores$Site[rownames(AllScores) %in% RightSamples] <- "Right"
AllScores$PatId <- substr(rownames(AllScores),1, 
                          nchar(rownames(AllScores))-12)

AllScores$Cpe <- CpeTable$CPE[match(AllScores$PatId,CpeTable$Sample.ID)]


AllScoresPiled <- melt(AllScores, value.name="Site")
colnames(AllScoresPiled) <- c("Site", "Barcode", "ScoreType", "Score")
AllScoresPiled <- na.omit(AllScoresPiled)

#Boxplot
ggplot(AllScoresPiled, aes(x=ScoreType, y=Score, color=Site ))+ 
   geom_boxplot()+
   stat_compare_means(method = "wilcox.test")

################################################################################

#Subsetting for H&E inference model

###############################################################################
TopRight <-  ImmunoRight[order(ImmunoRight$StromaScore, 
                               decreasing = TRUE)[1:10], ]

DownLeft <- ImmunoLeft[order(ImmunoLeft$StromaScore,
                             decreasing = F)[c(1:10)],]

cat(rownames(TopRight), sep = ",")
cat(rownames(DownLeft), sep = ",")


CpeTable$Sample.ID <-as.character(CpeTable$Sample.ID)
CpeTable$ShortId <- substr(CpeTable$Sample.ID, 1, 
                           nchar(CpeTable$Sample.ID) -4)

BreastSite$Cpe <- CpeTable$CPE[match(BreastSite$PatID,
                              CpeTable$ShortId)]

TopRight$PatId <- BreastSite$PatID[match(rownames(TopRight), 
                                  BreastSite$ID)]
DownLeft$PatId <- BreastSite$PatID[match(rownames(DownLeft), 
                                  BreastSite$ID)]

TopRight$Site <- BreastSite$SiteDico[match(rownames(TopRight), 
                               BreastSite$ID)]

DownLeft$Site <- BreastSite$SiteDico[match(rownames(DownLeft), 
                                        BreastSite$ID)]
HoverSamples <- rbind(DownLeft,TopRight)

HoverSamples$Cpe <- BreastSite$Cpe[match(HoverSamples$PatId,
                                    BreastSite$PatID)]
###############################################################################

### Cell type Counts per samples from HoverNet (external Data)

###############################################################################
HoverNetRes <- read.delim("c:/route/to/CountsPerSample.tsv")

HoverNetRes$PatId <-substr(HoverNetRes$Sample,1, 
                             nchar(HoverNetRes$Sample)-48)


StromalTotalData <- merge(HoverNetRes, HoverSamples, by.x=11, by.y=68)
StromalTotalData <- StromalTotalData[,-c(79)]


write.table(StromalTotalData,
            file = "StromalTotalData.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


###############################################################################

#####               Purity Correlation Analysis

###############################################################################


#### CPE Score vs Stromal Relative Abundance (HoverNet)

#Pearson correlation 
CorTest <- cor.test(StromalTotalData$Cpe, StromalTotalData$RelStroma, 
                    method = "pearson")


#ScatterPlot
ggplot(StromalTotalData, aes(x = Cpe, y = RelStroma)) +
  geom_point(size = 3, alpha = 0.7, color = "#1f77b4") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  stat_cor(method = "pearson",
           label.x.npc = "left",   
           label.y.npc = "top",    
           size = 4) +             
  theme_classic(base_size = 14) +
  labs(
    x = "Tumor purity (CPE)",
    y = "Relative stroma (RelStroma)",
    title = "Correlation between stromal content and tumor purity"
  )


### CPE Score vs StromaScore(xCell)

# Pearson Correlation
CorTest2 <- cor.test(StromalTotalData$Cpe, StromalTotalData$StromaScore, 
                     method = "pearson")


#Pearson correlation for all IDCs CPE vs StromalScore-XcellBased

#ScatterPlot
ggplot(StromalTotalData, aes(x = Cpe, y = StromaScore)) +
  geom_point(size = 3, alpha = 0.7, color = "#1f77b4") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  stat_cor(method = "pearson",
           label.x.npc = "left",   
           label.y.npc = "bottom",    
           size = 4) +             
  theme_classic(base_size = 14) +
  labs(
    x = "Tumor purity (CPE)",
    y = "StromaScore xCell-based",
    title = "Correlation between stromal content and tumor purity"
  )



# correlación de Pearson
CorTestXcellVCpe <- cor.test(AllScores$Cpe, AllScores$StromaScore, 
                             method = "pearson")
#ScatterPlot
ggplot(AllScores, aes(x = Cpe, y = StromaScore)) +
  geom_point(size = 3, alpha = 0.7, color = "#1f77b4") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  stat_cor(method = "pearson",
           label.x.npc = "left",   # posición relativa en el gráfico
           label.y.npc = "top",    # arriba del todo
           size = 4) +             # tamaño de letra
  theme_classic(base_size = 14) +
  labs(
    x = "Tumor purity (CPE)",
    y = "StromaScore xCell based",
    title = "Correlation between stromal content and tumor purity"
  )



###############################################################################
####  Differential Gene Expression (Purity Ranges 330 Ion-channel genes)
###############################################################################

CpeTablePures <- CpeTable[CpeTable$CPE >= 0.7,]

#Match samples in purity table with samples in IDCs table
IdcPures<- ClinicalIDC[ClinicalIDC$bcr_patient_barcode %in% CpeTablePures$ShortId,]

# Remove Male samples
IdcPures <- IdcPures[IdcPures$gender =="FEMALE",]


CpeIdcPures <- CpeTablePures[CpeTablePures$ShortId %in% IdcPures$bcr_patient_barcode,]


CpeIdcPures$ID <- Samples$ID[match(CpeIdcPures$ShortId, Samples$PatID)]


# 0.9 - 1
PuresHigh <- subset(CpeIdcPures, CPE >= 0.9 & CPE <= 1)

# 0.8 - 0.89
PuresMid <- subset(CpeIdcPures, CPE >= 0.8 & CPE < 0.9)

# 0.7 - 0.79
PuresLow <- subset(CpeIdcPures, CPE >= 0.7 & CPE < 0.8)

# 330 ICH genes 
genes <- unlist(strsplit("HTR3A HTR3B HTR3C HTR3D HTR3E ASIC1 ASIC2 ASIC3 ASIC4
ASIC5 ANO1 ANO2 ANO3 ANO4 ANO5 ANO6 ANO7 ANO8 ANO9 ANO10 MIP AQP1 AQP2 AQP3 
AQP4 AQP5 AQP6 AQP7 AQP8 AQP9 AQP10 AQP11 AQP12A AQP12B BEST1 BEST2 BEST3 BEST4 
CACNG1 CACNG2 CACNG3 CACNG4 CACNG5 CACNG6 CACNG7 CACNG8 CACNA1A CACNA1B CACNA1C 
CACNA1D CACNA1E CACNA1F CACNA1G CACNA1H CACNA1I CACNA1S CACNA2D1 CACNA2D2 
CACNA2D3 CACNA2D4 CACNB1 CACNB2 CACNB3 CACNB4 CATSPERB CATSPERD CATSPERE 
CATSPERG CATSPERZ CATSPER1 CATSPER2 CATSPER3 CATSPER4 CFTR CLIC1 CLIC2 CLIC3 
CLIC4 CLIC5 CLIC6 CLCNKA CLCNKB BSND CLCN1 CLCN2 CLCN3 CLCN4 CLCN5 CLCN6 CLCN7 
CHRNA1 CHRNA2 CHRNA3 CHRNA4 CHRNA5 CHRNA6 CHRNA7 CHRNA9 CHRNA10 CHRNB1 CHRNB2 
CHRNB3 CHRNB4 CHRND CHRNE CHRNG CNGA1 CNGA2 CNGA3 CNGA4 CNGB1 CNGB3 HCN1 HCN2 
HCN3 HCN4 GABRA1 GABRA2 GABRA3 GABRA4 GABRA5 GABRA6 GABRB1 GABRB2 GABRB3 GABRD 
GABRE GABRG1 GABRG2 GABRG3 GABRP GABRQ GABRR1 GABRR2 GABRR3 GJA1 GJA3 GJA4 GJA5 
GJA6P GJA8 GJA9 GJA10 GJB1 GJB2 GJB3 GJB4 GJB5 GJB6 GJB7 GJC1 GJC2 GJC3 GJD2 
GJD3 GJD4 GJE1 GRIA1 GRIA2 GRIA3 GRIA4 GRID1 GRID2 GRIK1 GRIK2 GRIK3 GRIK4 
GRIK5 GRIN1 GRIN2A GRIN2B GRIN2C GRIN2D GRIN3A GRIN3B GLRA1 GLRA2 GLRA3 GLRA4 
GLRB HVCN1 ITPR1 ITPR2 ITPR3 KCNMA1 KCNN1 KCNN2 KCNN3 KCNN4 KCNU1 KCNJ1 KCNJ2 
KCNJ3 KCNJ4 KCNJ5 KCNJ6 KCNJ8 KCNJ9 KCNJ10 KCNJ11 KCNJ12 KCNJ13 KCNJ14 KCNJ15 
KCNJ16 KCNJ18 KCNT1 KCNT2 KCNK1 KCNK2 KCNK3 KCNK4 KCNK5 KCNK6 KCNK7 KCNK9 
KCNK10 KCNK12 KCNK13 KCNK15 KCNK16 KCNK17 KCNK18 KCNA1 KCNA2+ KCNA3 KCNA4 KCNA5 
KCNA6 KCNA7 KCNA10 KCNB1 KCNB2 KCNC1 KCNC2 KCNC3 KCNC4 KCND1 KCND2 KCND3 KCNF1 
KCNG1 KCNG2 KCNG3 KCNG4 KCNH1 KCNH2 KCNH3 KCNH4 KCNH5 KCNH6 KCNH7 KCNH8 KCNQ1 
KCNQ2 KCNQ3 KCNQ4 KCNQ5 KCNS1 KCNS2 KCNS3 KCNV1 KCNV2 P2RX1 P2RX2 P2RX3 P2RX4 
P2RX5 P2RX6 P2RX7 RYR1 RYR2 RYR3 SCNN1A SCNN1B SCNN1D SCNN1G NALCN SCN1A SCN2A 
SCN3A SCN4A SCN5A SCN8A SCN9A SCN10A SCN11A SCN1B SCN2B SCN3B SCN4B TRPA1 TRPC1 
TRPC2 TRPC3 TRPC4 TRPC5 TRPC6 TRPC7 MCOLN1 MCOLN2 MCOLN3 TRPM1 TRPM2 TRPM3 TRPM4 
TRPM5 TRPM6 TRPM7 TRPM8 PKD2 PKD2L1 PKD2L2 TRPV1 TRPV2 TRPV3 TRPV4 TRPV5 TRPV6 
TPCN1 TPCN2 VDAC1 VDAC2 VDAC3 LRRC8A LRRC8B LRRC8C LRRC8D LRRC8E ZACN", "\\s+"))


#RNA seq data for every interval
IchRnaData <- RawData[rownames(RawData) %in% genes,]

IchRnaData <- data.frame(t(IchRnaData))

IchRnaData$Site <- BreastSite$SiteDico[match(rownames(IchRnaData),
                                             BreastSite$ID)]

IchRnaData <- na.omit(IchRnaData) 

RnaHigh <- IchRnaData[rownames(IchRnaData) %in% PuresHigh$ID,]

RnaMid <- IchRnaData[rownames(IchRnaData) %in% PuresMid$ID,]

RnaLow <- IchRnaData[rownames(IchRnaData) %in% PuresLow$ID,]

#Metadata from intervals
MetaHigh <- data.frame( row.names = rownames(RnaHigh), Site = RnaHigh$Site )

MetaMid <- data.frame( row.names = rownames(RnaMid), Site = RnaMid$Site )

MetaLow <- data.frame( row.names = rownames(RnaLow), Site = RnaLow$Site )


#MAtrix conversion
CountsHigh <- RnaHigh[,!colnames(RnaHigh) %in% "Site"]

CountsHigh <- apply(CountsHigh, 2, as.numeric)

rownames(CountsHigh) <- rownames(MetaHigh) 

CountsMid <- RnaMid[,!colnames(RnaMid) %in% "Site"]

CountsMid <- apply(CountsMid, 2, as.numeric)

rownames(CountsMid) <- rownames(MetaMid)

CountsLow <- RnaLow[,!colnames(RnaLow) %in% "Site"]

CountsLow <- apply(CountsLow, 2, as.numeric)

rownames(CountsLow) <- rownames(MetaLow)

# DESeq object construction

ObjHigh <- DESeqDataSetFromMatrix( countData = t(CountsHigh),
                                   colData = MetaHigh,
                                   design = ~Site)


ObjMid <- DESeqDataSetFromMatrix( countData = t(CountsMid),
                                   colData = MetaMid,
                                   design = ~Site)

ObjLow <- DESeqDataSetFromMatrix( countData = t(CountsLow),
                                   colData = MetaLow,
                                   design = ~Site)
#######################
#DGE Analysis (Interval High)
######################3
DeHigh <- DESeq(ObjHigh)
ResHigh <- results(DeHigh, contrast = c("Site", "Left", "Right"))

ResHigh <- as.data.frame(ResHigh)

ResHigh <- ResHigh[order(ResHigh$pvalue), ]


EnhancedVolcano(ResHigh,
                lab = rownames(ResHigh),   
                x = 'log2FoldChange',        
                y = 'pvalue',                
                pCutoff = 0.05,              
                FCcutoff = 1.0,              
                pointSize = 2.0,
                labSize = 4.0,
                title = "Volcano plot Left vs Right",
                subtitle = "CPE Range 1-0.9",
                legendPosition = 'right')

#######################
#DGE Analysis (Interval Middle)
######################


DeMid <- DESeq(ObjMid)
ResMid <- results(DeMid, contrast = c("Site", "Left", "Right"))

ResMid <- as.data.frame(ResMid)

ResMid <- ResHigh[order(ResMid$pvalue), ]


EnhancedVolcano(ResMid,
                lab = rownames(ResMid),   
                x = 'log2FoldChange',        
                y = 'pvalue',                
                pCutoff = 0.05,              
                FCcutoff = 1.0,             
                pointSize = 2.0,
                labSize = 4.0,
                title = "Volcano plot Left vs Right",
                subtitle = "CPE Range 0.8-0.89",
                legendPosition = 'right')

#######################
#DGE Analysis (Interval Low)
######################


DeLow <- DESeq(ObjLow)
ResLow <- results(DeLow, contrast = c("Site", "Left", "Right"))

ResLow <- as.data.frame(ResLow)

ResLow <- ResLow[order(ResLow$pvalue), ]


EnhancedVolcano(ResLow,
                lab = rownames(ResLow),   
                x = 'log2FoldChange',       
                y = 'pvalue',                
                pCutoff = 0.05,              
                FCcutoff = 1.0,              
                pointSize = 2.0,
                labSize = 4.0,
                title = "Volcano plot Left vs Right",
                subtitle = "CPE Range 0.7-0.79",
                legendPosition = 'right')

#Save DGE Results
ResHigh$Gene <- rownames(ResHigh)
ResMid$Gene <- rownames(ResMid)
ResLow$Gene <- rownames(ResLow)
  
write.table(ResHigh,
            file = "DEGIntervalHigh.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(ResMid,
            file = "DEGIntervalMid.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


write.table(ResLow,
            file = "DEGIntervalLow.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


###############################################################################
##  DGE with the most 30% pure IDCs
###############################################################################

#Aquaporines and gap junction genes
AquaGapGenes <- c("AQP1", "AQP7", "AQP8", "AQP11", "GJA1", "GJB1", "GJB2", 
                  "GJC2", "GJD3", "GJA9") 

#Ich Genes

IchGenes <- c("CACNA1A", "CACNB2", "CACNA2D2", "KCNJ11", "SCN3A", "SCN3B")

# IDC with CPE >= 0.7 and Sub-dataframe (CpeIdcPures)

PureRnaData <- RawData[rownames(RawData) %in% genes,]

PureRnaData <- data.frame(t(PureRnaData))

PureRnaData$Site <- BreastSite$SiteDico[match(rownames(PureRnaData),
                                             BreastSite$ID)]

PureRnaData$PatId <-substr(rownames(PureRnaData),1, 
                           nchar(rownames(PureRnaData))-12)

PureRnaData <- PureRnaData[PureRnaData$PatId %in% CpeIdcPures$Sample.ID,]


PureRnaData <- na.omit(PureRnaData) 

PureRnaData <- PureRnaData[,-c(1)]

#Metadata from intervals
MetaPure <- data.frame( row.names = rownames(PureRnaData), 
                        Site = PureRnaData$Site )


#MAtrix conversion
CountsPure <- PureRnaData[,!colnames(PureRnaData) %in% "Site"]

CountsPure <- CountsPure[,!colnames(CountsPure) %in% "PatId"]

CountsPure <- apply(CountsPure, 2, as.numeric)

rownames(CountsPure) <- rownames(MetaPure) 

CountsPure <- as.matrix(CountsPure)

# Filter low expressed genes
Keep <- rowSums(CountsPure >= 10) > 0
CountsPure <- CountsPure[Keep, ]


# DESeq object construction

ObjPure <- DESeqDataSetFromMatrix( countData = round(t(CountsPure)),
                                   colData = MetaPure,
                                   design = ~Site)




DePure <- DESeq(ObjPure)
ResPure <- results(DePure, contrast = c("Site", "Left", "Right"))

ResPure <- as.data.frame(ResPure)

ResPure <- ResPure[order(ResPure$pvalue), ]


EnhancedVolcano(ResPure,
                lab = rownames(ResPure),   
                x = 'log2FoldChange',        
                y = 'pvalue',                
                pCutoff = 0.05,              
                FCcutoff = 1.0,              
                pointSize = 2.0,
                labSize = 4.0,
                title = "Volcano plot Left vs Right",
                subtitle = "CPE Range 1-0.9",
                legendPosition = 'right')

ResPure$Gene <- rownames(ResPure)

AquaPure <- ResPure[ResPure$Gene %in% AquaGapGenes,]
IchPure <- ResPure[ResPure$Gene %in% IchGenes,] 



################################################################################
####            PAM50Subtypes by side Left/Right
################################################################################


SubtypesSamples <- BreastSite

SubtypesSamples$Pam50 <-Subtypes$BRCA_Subtype_PAM50[match(SubtypesSamples$PatID,
                                                          Subtypes$patient)]

SubtypesSamples <- na.omit(SubtypesSamples)



# Proportions by subtype and location
ProporTable <- SubtypesSamples %>%
  group_by(Pam50, SiteDico) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Pam50) %>%
  mutate(prop = n / sum(n))

# Fisher test by subtype
PvalTable <- SubtypesSamples %>%
  group_by(Pam50) %>%
  summarise(p_value = {
    tab <- table(SiteDico)
    if(any(tab < 5)) {
      fisher.test(tab)$p.value
    } else {
      chisq.test(tab)$p.value
    }
  }) %>%
  mutate(p_label = paste0("p = ", signif(p_value, 3)),
         y_pos = 1.05)

# Plot
p <- ggplot(ProporTable, aes(x = Pam50, y = prop, fill = SiteDico)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_y_continuous(labels = percent_format(accuracy = 1), 
                     expand = expansion(c(0, 0.15))) +
  scale_fill_manual(values = c("Left" = "#E69F00", "Right" = "#56B4E9")) +
  labs(
    x = "PAM50 Subtype",
    y = "Proportion of samples (%)",
    fill = "Side"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.margin = margin(10, 10, 20, 10)
  ) +
  geom_text(data = PvalTable, aes(x = Pam50, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 3)

print(p)

###############################################################################

#### CAFs population comparison by Pam50 subtype

###############################################################################

## For total IDCs

SubtypesSamples$StromaScore <- AllScores$StromaScore[match(SubtypesSamples$ID,
                                                           rownames(AllScores))]
#CAF-related Signatures (Cords et al)

FileSign<- readxl::read_xlsx("C:/route/to/CafSign.xlsx")

###Formatting
CafSign <- lapply(FileSign$`TOP 6 GENES FROM Cords et al (Ref 28)`, function(x){
  strsplit(x, ",\\s*")[[1]] 
})
names(CafSign) <- FileSign$`CAF subtype`

CafRnaData <- RawData[,colnames(RawData) %in% rownames(AllScores)]

#GSVA signature analysis
CafScores <- gsva(as.matrix(CafRnaData), CafSign, method = "ssgsea")

CafScores <- data.frame(t(CafScores))

CafScores$Site <- AllScores$Site[match(rownames(CafScores), 
                                       rownames(AllScores))]
CafScores$Pam50 <- SubtypesSamples$Pam50[match(rownames(CafScores), 
                                               SubtypesSamples$ID)]

CafScores <- na.omit(CafScores)

colnames(CafScores) <- c("iCAF", "dCAF", "vCAF", "mCAF", "hspCAF", "rCAF", 
                         "tCAF", "Site", "Pam50")

CafScores$PatId <- substr(rownames(CafScores), 1,12)

#  CafScores, PAM50 subtype and Location Data to long format

CafScoresLong <- CafScores %>%
  pivot_longer(cols = c(iCAF, dCAF, vCAF, mCAF, hspCAF, rCAF,tCAF),
               names_to = "Signature",
               values_to = "Score")


CafSummary <- CafScoresLong %>%
  group_by(Pam50, Site, Signature) %>%
  summarise(MeanScore = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  group_by(Pam50, Site) %>%
  mutate(Proportion = MeanScore / sum(MeanScore))

CafTests <- CafScoresLong %>%
  group_by(Pam50, Signature) %>%
  summarise(
    p_value = tryCatch(t.test(Score ~ Site)$p.value, 
                       error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))

#Barplot

ggplot(CafSummary, aes(x = Site, y = Proportion, fill = Signature)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Pam50) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Proporción relativa de CAF", x = "Sitio") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  geom_text(
    data = CafTests,
    aes(x = 1.5, y = 1.05, label = signif),
    inherit.aes = FALSE
  )




CafSummary2 <- CafScoresLong %>%
  group_by(Pam50, Signature, Site) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    sd_score = sd(Score, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

#Shapiro Test
CafScoresLong %>% group_by(Pam50, Signature, Site) %>%
  summarise(shapiro_p = shapiro.test(Score)$p.value)


CafTests2 <- CafScoresLong %>%
  group_by(Pam50, Signature) %>%
  summarise(
    p_value = tryCatch(wilcox.test(Score ~ Site)$p.value, 
                       error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))


###############################################################################
####PLOT CODE FOR PUBLICATION
###############################################################################


# CAF scores by subtypes

MeanCafScoresLongDf <- CafScoresLong %>%
  group_by(Pam50, Site, Signature) %>%
  summarise(mean_CAF = mean(Score, na.rm = TRUE), .groups = "drop")

MedianCafScoreLongDf <- CafScoresLong %>%
  group_by(Pam50, Site, Signature) %>%
  summarise(median_CAF = median(Score, na.rm = TRUE), .groups = "drop")

#Boxplot

ggplot(CafScoresLong, aes(x = Signature, y = Score, fill = Site)) +
  geom_boxplot(position = position_dodge(0.8)) +
  facet_wrap(~Pam50) +
  labs(y = "CAF Score (GSVA)", x = "CAF", fill = "Location") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  geom_text(
    data = CafTests2,
    aes(x = Signature, y = max(CafScoresLong$Score) + 0.05, label = signif),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.8)
  )+
  scale_fill_manual(
    values = c(
      "Left"  = "#2CB1B5",
      "Right" = "#C51B7D"
    )
  )

###############################################################################

#Save df data

write.table(CafSummary2,
            file = "C:/route/to/save/CafSummary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


write.table(CafTests2,
            file = "C:/route/to/save/CafStatsV2.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

## For IDCs below 0.7 CPE

ImpSamples <- AllScores[AllScores$Cpe < 0.7,]

CafRnaDataImp <- RawData[,colnames(RawData) %in% rownames(ImpSamples)]

#GSVA signature analysis
CafScoresImp <- gsva(as.matrix(CafRnaDataImp), CafSign, method = "ssgsea")

CafScoresImp <- data.frame(t(CafScoresImp))

CafScoresImp$Site <- AllScores$Site[match(rownames(CafScoresImp), 
                                       rownames(AllScores))]
CafScoresImp$Pam50 <- SubtypesSamples$Pam50[match(rownames(CafScoresImp), 
                                               SubtypesSamples$ID)]

CafScoresImp <- na.omit(CafScoresImp)

colnames(CafScoresImp) <- c("iCAF", "dCAF", "vCAF", "mCAF", "hspCAF", "rCAF", 
                         "tCAF", "Site", "Pam50")

CafScoresImp$PatId <- substr(rownames(CafScoresImp), 1,12)

#  CafScores, PAM50 subtype and Location Data to long format

CafScoresLongImp <- CafScoresImp %>%
  pivot_longer(cols = c(iCAF, dCAF, vCAF, mCAF, hspCAF, rCAF,tCAF),
               names_to = "Signature",
               values_to = "Score")


CafSummaryImp <- CafScoresLongImp %>%
  group_by(Pam50, Site, Signature) %>%
  summarise(MeanScore = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  group_by(Pam50, Site) %>%
  mutate(Proportion = MeanScore / sum(MeanScore))

CafTestsImp <- CafScoresLongImp %>%
  group_by(Pam50, Signature) %>%
  summarise(
    p_value = tryCatch(t.test(Score ~ Site)$p.value, 
                       error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))


ggplot(CafSummaryImp, aes(x = Site, y = Proportion, fill = Signature)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Pam50) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Proporción relativa de CAF", x = "Sitio") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  geom_text(
    data = CafTestsImp,
    aes(x = 1.5, y = 1.05, label = signif),
    inherit.aes = FALSE
  )




CafSummaryImp2 <- CafScoresLongImp %>%
  group_by(Pam50, Signature, Site) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    sd_score = sd(Score, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )



CafTestsImp2 <- CafScoresLongImp %>%
  group_by(Pam50, Signature) %>%
  summarise(
    p_value = tryCatch(wilcox.test(Score ~ Site)$p.value, 
                       error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))


ggplot(CafScoresLongImp, aes(x = Signature, y = Score, fill = Site)) +
  geom_boxplot(position = position_dodge(0.8)) +
  facet_wrap(~Pam50) +
  labs(y = "CAF Score (GSVA)", x = "CAF", fill = "Location") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  geom_text(
    data = CafTestsImp2,
    aes(x = Signature, y = max(CafScoresLong$Score) + 0.05, label = signif),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.8)
  )


write.table(CafSummaryImp2,
            file = "C:/route/to/save/CafSummaryImpureIDC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)



write.table(CafTestsImp2,
            file = "C:/route/to/save/CafStatsImp.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


###############################################################################

## Total CAF signature Left vs Right by Pam50

###############################################################################
CafTotal <- CafScores %>%
  rowwise() %>%  
  mutate(CAF_total = sum(c(iCAF, dCAF, vCAF, mCAF, hspCAF, rCAF, tCAF), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(PatId, Site, Pam50, CAF_total)

#Statistical Data

CafTotalTests <- CafTotal %>%
  group_by(Pam50) %>%
  summarise(
    p_value = tryCatch(wilcox.test(CAF_total ~ Site)$p.value, 
                       error = function(e) NA_real_),.groups = "drop"
  ) %>%
  mutate(signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))



#############################################################♥

#### PLOT CODE FOR PUBLICATION#####

#############################################################

MeanDf <- CafTotal %>%
  group_by(Pam50, Site) %>%
  summarise(mean_CAF = mean(CAF_total, na.rm = TRUE), .groups = "drop")

MedianDf <- CafTotal %>%
  group_by(Pam50, Site) %>%
  summarise(median_CAF = median(CAF_total, na.rm = TRUE), .groups = "drop")


#Empty Plot

ggplot(CafTotal, aes(x = Site, y = CAF_total, fill = Site)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black") +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~Pam50) +
  labs(y = "CAF total", x = "Site") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "Left"  = "#2CB1B5",
      "Right" = "#C51B7D"
    )
  )




##############################################################################




write.table(CafTotalTests,
            file = "C:/route/to/save/CafTotalStats.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

###########################################################################

## Total CAF signature Left vs Right by Pam50 for Idc with CPE below 0.7

#############################################################################

#Formatting
CafTotalImp <- CafScoresImp %>%
  rowwise() %>%  
  mutate(CAF_total = sum(c(iCAF, dCAF, vCAF, mCAF, hspCAF, rCAF, tCAF), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(PatId, Site, Pam50, CAF_total)

#Statistical Data

CafTotalTestsImp <- CafTotalImp %>%
  group_by(Pam50) %>%
  summarise(
    p_value = tryCatch(wilcox.test(CAF_total ~ Site)$p.value, 
                       error = function(e) NA_real_),.groups = "drop"
  ) %>%
  mutate(signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))


#############################################################♥

#### PLOT CODE FOR PUBLICATION#####

#############################################################


MeanDfImp <- CafTotalImp %>%
  group_by(Pam50, Site) %>%
  summarise(mean_CAF = mean(CAF_total, na.rm = TRUE), .groups = "drop")

MedianDfImp <- CafTotalImp %>%
  group_by(Pam50, Site) %>%
  summarise(median_CAF = median(CAF_total, na.rm = TRUE), .groups = "drop")

#Empty Plot

ggplot(CafTotalImp, aes(x = Site, y = CAF_total, fill = Site)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black") +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~Pam50) +
  labs(y = "CAF total", x = "Site") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "Left"  = "#2CB1B5",
      "Right" = "#C51B7D"
    )
  )




##############################################################################


write.table(CafTotalTestsImp,
            file = "C:/route/to/save/CafTotalStatsImp.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)





##############################################################################

## STROMASCORE ALL IDCs BY PAM50 PLOT CODE FOR PUBLICATION

##############################################################################

MeanDfSubtypes <- SubtypesSamples %>%
  group_by(Pam50, SiteDico) %>%
  summarise(mean_score = mean(StromaScore, na.rm = TRUE), .groups = "drop")

MedianDfSubtypes <- SubtypesSamples %>%
  group_by(Pam50, SiteDico) %>%
  summarise(median_score = median(StromaScore, na.rm = TRUE), .groups = "drop")


#Empty Plot

ggplot(SubtypesSamples, aes(x = SiteDico, y = StromaScore, fill = SiteDico)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black") +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~Pam50) +
  labs(y = "Stroma Score", x = "Site") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "Left"  = "#2CB1B5",
      "Right" = "#C51B7D"
    )
  )


##############################################################################

## STROMASCORE  IDCs CPE BELOW 0.7 BY PAM50 PLOT CODE FOR PUBLICATION

##############################################################################

MeanDfSubtypesImp <- SubtypesSamples[SubtypesSamples$Cpe > 0.7, ] %>%
  group_by(Pam50, SiteDico) %>%
  summarise(mean_score = mean(StromaScore, na.rm = TRUE), .groups = "drop")

MedianDfSubtypesImp <- SubtypesSamples[SubtypesSamples$Cpe > 0.7, ] %>%
  group_by(Pam50, SiteDico) %>%
  summarise(median_score = median(StromaScore, na.rm = TRUE), .groups = "drop")


#Empty Plot

ggplot(SubtypesSamples[SubtypesSamples$Cpe > 0.7,], aes(x = SiteDico, y = StromaScore, fill = SiteDico)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black") +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~Pam50) +
  labs(y = "Stroma Score", x = "Site") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "Left"  = "#2CB1B5",
      "Right" = "#C51B7D"
    )
  )



XenaSamples <- BreastSite

XenaSamples$XenaId <- substr(XenaSamples$ID,1, 
                             nchar(XenaSamples$ID)-13)
###############################################################################

#ID List obtained from (Xena)

###############################################################################

XenaIdcData <- readxl::read_xlsx("c:/route/to/XenaDataIDCs.xlsx")
XenaIdcData <- XenaIdcData[,-c(1)]
colnames(XenaIdcData) <- c("Samples", "Type", "StromaScore")
XenaIdcData$Site <- XenaSamples$SiteDico[match(XenaIdcData$Samples, 
                                               XenaSamples$XenaId)]
XenaIdcData$CPE <- XenaSamples$Cpe[match(XenaIdcData$Samples, 
                                         XenaSamples$XenaId)]

XenaIdcDataCaf <- readxl::read_xlsx("c:/route/to/XenaIdcCafSign.xlsx")
XenaIdcDataCaf <- XenaIdcDataCaf[,-c(1)]
colnames(XenaIdcDataCaf) <- c("Samples", "Type", "CafScore")

XenaIdcData$CafScore <- XenaIdcDataCaf$CafScore[match(XenaIdcData$Samples, 
                                                      XenaIdcDataCaf$Samples)]


###############################################################################
#StromaScore ALL IDCs plot
###############################################################################

PWilcox <- wilcox.test(
  as.numeric(StromaScore) ~ Site,
  data = XenaIdcData
)$p.value

PTest <- t.test(
  as.numeric(StromaScore) ~ Site,
  data = XenaIdcData
)$p.value

PLabel  <- paste0("p = ", formatC(PWilcox, format = "f", digits = 3))
PLabel2  <- paste0("p = ", formatC(PTest, format = "f", digits = 3))

# Formatting to numeric and factors
XenaIdcData$Site <- factor(XenaIdcData$Site, levels = c("Left", "Right"))
XenaIdcData$StromaScore <- as.numeric(XenaIdcData$StromaScore)

# Dynamic Limits (graphical display)
y_max <- max(XenaIdcData$StromaScore, na.rm = TRUE)
y_bar <- y_max * 1.03
y_txt <- y_max * 1.06

ggplot(XenaIdcData, aes(x = Site, y = StromaScore, fill = Site)) +
  
  ## Violín plot
  geom_violin(
    trim = FALSE,
    color = "black",
    linewidth = 0.6
  ) +
  geom_boxplot(
    width = 0.08,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "Left"  = "#2CB1B5",
      "Right" = "#C51B7D"
    )
  ) +
  annotate(
    "segment",
    x = 1, xend = 2,
    y = y_bar, yend = y_bar,
    linewidth = 0.6
  ) +
  annotate(
    "segment",
    x = c(1, 2), xend = c(1, 2),
    y = y_bar, yend = y_bar - (y_max * 0.01),
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 1.5,
    y = y_txt,
    label = PLabel2,
    size = 4
  ) +
  labs(
    x = NULL,
    y = "Stromal Score RNA-seq based"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 12)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)))+
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.6
    ),
    panel.grid.minor.y = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

MeanIdcStromascoreDf <- XenaIdcData %>%
  group_by(Site) %>%
  summarise(mean_score = mean(StromaScore, na.rm = TRUE), .groups = "drop")

###############################################################################
#StromaScore  IDCs CPE 0.7 plot
###############################################################################

XenaIdcDataFilt <- XenaIdcData[XenaIdcData$CPE <= 0.7,]
XenaIdcDataFilt <- na.omit(XenaIdcDataFilt)

#WilcoxTest

PWilcoxFilt <- wilcox.test(
  as.numeric(StromaScore) ~ Site,
  data = XenaIdcDataFilt
)$p.value

#T test

PTestFilt <- t.test(
  as.numeric(StromaScore) ~ Site,
  data = XenaIdcDataFilt
)$p.value

PLabelFilt  <- paste0("p = ", formatC(PWilcoxFilt, format = "f", digits = 3))
PLabelFilt2  <-paste0("p = ", sprintf("%.3f", trunc(0.0068 * 1000) / 1000))


XenaIdcDataFilt$Site <- factor(XenaIdcDataFilt$Site, levels = c("Left", "Right"))
XenaIdcDataFilt$StromaScore <- as.numeric(XenaIdcDataFilt$StromaScore)


y_max <- max(XenaIdcDataFilt$StromaScore, na.rm = TRUE)
y_bar <- y_max * 1.03
y_txt <- y_max * 1.06

ggplot(XenaIdcDataFilt, aes(x = Site, y = StromaScore, fill = Site)) +
  geom_violin(
    trim = FALSE,
    color = "black",
    linewidth = 0.6
  ) +
  geom_boxplot(
    width = 0.08,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "Left"  = "#2CB1B5",
      "Right" = "#C51B7D"
    )
  ) +
  annotate(
    "segment",
    x = 1, xend = 2,
    y = y_bar, yend = y_bar,
    linewidth = 0.6
  ) +
  annotate(
    "segment",
    x = c(1, 2), xend = c(1, 2),
    y = y_bar, yend = y_bar - (y_max * 0.01),
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 1.5,
    y = y_txt,
    label = PLabelFilt2,
    size = 4
  ) +
  labs(
    x = NULL,
    y = "Stromal Score RNA-seq based"
  ) +
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 12)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)))+
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.6
    ),
    panel.grid.minor.y = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

MeanIdcStromascoreDfFilt <- XenaIdcDataFilt %>%
  group_by(Site) %>%
  summarise(mean_score = mean(StromaScore, na.rm = TRUE), .groups = "drop")
###############################################################################

#CafScore ALL IDCs plot

###############################################################################

PWilcoxCaf <- wilcox.test(
  as.numeric(CafScore) ~ Site,
  data = XenaIdcData
)$p.value

PTestCaf <- t.test(
  as.numeric(CafScore) ~ Site,
  data = XenaIdcData
)$p.value

PLabelCaf  <- paste0("p = ", formatC(PWilcoxCaf, format = "f", digits = 3))
PLabelCaf2  <- paste0("p = ", formatC(PTestCaf, format = "f", digits = 3))

XenaIdcData$Site <- factor(XenaIdcData$Site, levels = c("Left", "Right"))
XenaIdcData$CafScore <- as.numeric(XenaIdcData$CafScore)


y_max <- max(XenaIdcData$CafScore, na.rm = TRUE)
y_bar <- y_max * 1.03
y_txt <- y_max * 1.06

ggplot(XenaIdcData, aes(x = Site, y = CafScore, fill = Site)) +
  
  
  geom_violin(
    trim = FALSE,
    color = "black",
    linewidth = 0.6
  ) +
  
  
  geom_boxplot(
    width = 0.08,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.6
  ) +
  
  
  scale_fill_manual(
    values = c(
      "Left"  = "#2CB1B5",
      "Right" = "#C51B7D"
    )
  ) +
 
  annotate(
    "segment",
    x = 1, xend = 2,
    y = y_bar, yend = y_bar,
    linewidth = 0.6
  ) +
  
  
  annotate(
    "segment",
    x = c(1, 2), xend = c(1, 2),
    y = y_bar, yend = y_bar - (y_max * 0.01),
    linewidth = 0.6
  ) +
  
  
  annotate(
    "text",
    x = 1.5,
    y = y_txt,
    label = PLabelCaf2,
    size = 4
  ) +
  
  
  labs(
    x = NULL,
    y = "Fibroblast Gene Expression Signature"
  ) +
  
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 12)
  ) +
  
  
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)))+
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    
    
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    
    
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.6
    ),
    panel.grid.minor.y = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

MeanIdcCafScoreDf <- XenaIdcData %>%
  group_by(Site) %>%
  summarise(mean_score = mean(CafScore, na.rm = TRUE), .groups = "drop")

###############################################################################

#CafScore  IDCs CPE 0.7 plot

###############################################################################
PWilcoxCafFilt <- wilcox.test(
  as.numeric(CafScore) ~ Site,
  data = XenaIdcDataFilt
)$p.value

PTestCafFilt <- t.test(
  as.numeric(CafScore) ~ Site,
  data = XenaIdcDataFilt
)$p.value

PLabelFiltCaf  <- paste0("p = ", formatC(PWilcoxCafFilt, format = "f", digits = 3))
PLabelFiltCaf2  <-paste0("p = ", formatC(PTestCafFilt, format = "f", digits = 3))

XenaIdcDataFilt$Site <- factor(XenaIdcDataFilt$Site, levels = c("Left", "Right"))
XenaIdcDataFilt$CafScore <- as.numeric(XenaIdcDataFilt$CafScore)


y_max <- max(XenaIdcDataFilt$CafScore, na.rm = TRUE)
y_bar <- y_max * 1.03
y_txt <- y_max * 1.06

ggplot(XenaIdcDataFilt, aes(x = Site, y = CafScore, fill = Site)) +
  
  
  geom_violin(
    trim = FALSE,
    color = "black",
    linewidth = 0.6
  ) +
  
  
  geom_boxplot(
    width = 0.08,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.6
  ) +
  
  
  scale_fill_manual(
    values = c(
      "Left"  = "#2CB1B5",
      "Right" = "#C51B7D"
    )
  ) +
  
  annotate(
    "segment",
    x = 1, xend = 2,
    y = y_bar, yend = y_bar,
    linewidth = 0.6
  ) +
  
  
  annotate(
    "segment",
    x = c(1, 2), xend = c(1, 2),
    y = y_bar, yend = y_bar - (y_max * 0.01),
    linewidth = 0.6
  ) +
  
  
  annotate(
    "text",
    x = 1.5,
    y = y_txt,
    label = PLabelFiltCaf2,
    size = 4
  ) +
  
  
  labs(
    x = NULL,
    y = "Fibroblast Gene Expresion Signature"
  ) +
  
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 12)
  ) +
  
  
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)))+
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    
    
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    
    
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.6
    ),
    panel.grid.minor.y = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

MeanIdcCafDfFilt <- XenaIdcDataFilt %>%
  group_by(Site) %>%
  summarise(mean_score = mean(CafScore, na.rm = TRUE), .groups = "drop")
###############################################################################





