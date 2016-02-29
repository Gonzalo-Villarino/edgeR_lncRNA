rm(list=ls())
setwd("~/Documents/NCSU/RNAseq_CMM_BL_Unique/Analysis_YFPs/HtSeq_reverse/")
getwd()
library("edgeR")
library("VennDiagram")
library("ggplot2")
library("gridExtra")


####################################################################################################################################
# EDGER (using htseq-counts previously)
# edgeR is based on Gener.Linear Model (GML) - assuming that read counts are distributed according to the negative binomial distribution
####################################################################################################################################

# Import (to select prot. coding genes only) 
gene_anno = read.table("TAIR10.gene.list.txt",sep="\t")

# Import the 14 htseq-count tables to start Diff.Expr.Analysis
YFP_NEG_B1T1_HtseqCounts <- read.table("YFP_NEG_B1T1_HtseqCounts.txt",header=F )
YFP_NEG_B1T2_HtseqCounts <- read.table("YFP_NEG_B1T2_HtseqCounts.txt",header=F )
YFP_NEG_B2T1_HtseqCounts <- read.table("YFP_NEG_B2T1_HtseqCounts.txt",header=F )
YFP_NEG_B2T2_HtseqCounts <- read.table("YFP_NEG_B2T2_HtseqCounts.txt",header=F )
YFP_NEG_B3T1_HtseqCounts <- read.table("YFP_NEG_B3T1_HtseqCounts.txt",header=F )
YFP_NEG_B3T2_HtseqCounts <- read.table("YFP_NEG_B3T2_HtseqCounts.txt",header=F )
YFP_NEG_B4T1_HtseqCounts <- read.table("YFP_NEG_B4T1_HtseqCounts.txt",header=F )
YFP_NEG_B4T2_HtseqCounts <- read.table("YFP_NEG_B4T2_HtseqCounts.txt",header=F )

YFP_POS_B1T1_HtseqCounts <- read.table("YFP_POS_B1T1_HtseqCounts.txt",header=F )
YFP_POS_B1T2_HtseqCounts <- read.table("YFP_POS_B1T2_HtseqCounts.txt",header=F )
YFP_POS_B2T1_HtseqCounts <- read.table("YFP_POS_B2T1_HtseqCounts.txt",header=F )
YFP_POS_B2T2_HtseqCounts <- read.table("YFP_POS_B2T2_HtseqCounts.txt",header=F )
YFP_POS_B3T1_HtseqCounts <- read.table("YFP_POS_B3T1_HtseqCounts.txt",header=F )
YFP_POS_B3T2_HtseqCounts <- read.table("YFP_POS_B3T2_HtseqCounts.txt",header=F )


# Merge single htseq-count files into one:
files <- c(
        "YFP_NEG_B1T1_HtseqCounts.txt",
        "YFP_NEG_B1T2_HtseqCounts.txt",
        "YFP_NEG_B2T1_HtseqCounts.txt",
        "YFP_NEG_B2T2_HtseqCounts.txt",
        "YFP_NEG_B3T1_HtseqCounts.txt",
        "YFP_NEG_B3T2_HtseqCounts.txt",
        "YFP_NEG_B4T1_HtseqCounts.txt",
        "YFP_NEG_B4T2_HtseqCounts.txt",
        "YFP_POS_B1T1_HtseqCounts.txt",
        "YFP_POS_B1T2_HtseqCounts.txt",
        "YFP_POS_B2T1_HtseqCounts.txt",
        "YFP_POS_B2T2_HtseqCounts.txt",
        "YFP_POS_B3T1_HtseqCounts.txt",
        "YFP_POS_B3T2_HtseqCounts.txt"
)
c <- 1
for (filename in files) {
        if(c == 1){ # if it is the first file just read file
                htseqAllCount_BL = read.table(filename,sep="\t")
        }
        else{ # else merge the other files
                tmp = read.table(filename,sep="\t")
                names(tmp) = c(c,c+1)
                htseqAllCount_BL = merge(htseqAllCount_BL,tmp,by=c(1))
        }
        c = c+1
}
names(htseqAllCount_BL) = c("GeneID","YFP_NEG_B1T1","YFP_NEG_B1T2","YFP_NEG_B2T1","YFP_NEG_B2T2","YFP_NEG_B3T1",
                            "YFP_NEG_B3T2","YFP_NEG_B4T1","YFP_NEG_B4T2","YFP_POS_B1T1","YFP_POS_B1T2", "YFP_POS_B2T1",
                            "YFP_POS_B2T2","YFP_POS_B3T1","YFP_POS_B3T2")
names(htseqAllCount_BL)



########################################################################
# merge technical replicates and store it into a new table
i = 2
j = 16
while(i<16){
        htseqAllCount_BL[,j] = htseqAllCount_BL[,i]+htseqAllCount_BL[,i+1]
        i = i+2
        j = j+1
}

htseqAllCount_BL_tech = htseqAllCount_BL[,c(1,16:22)]  ### new table that merge the technical replicates
names(htseqAllCount_BL_tech) = c("GeneID","YFP_NEG_B1","YFP_NEG_B2","YFP_NEG_B3","YFP_NEG_B4","YFP_POS_B1","YFP_POS_B2","YFP_POS_B3")

##########################################################################


###get protein coding genes###
#htseqAllCount_BL <- merge(htseqAllCount_BL,gene_anno,by=c(1))
#htseqAllCount_BL <- htseqAllCount_BL[htseqAllCount_BL$V2=="protein_coding_gene",]
#htseqAllCount_BL = htseqAllCount_BL[,-16]
####
str(htseqAllCount_BL)
names(htseqAllCount_BL)
          
# make table of only 6 numeric values for downstream analysis
cm <- htseqAllCount_BL_tech[,-1]
rownames(cm) <- htseqAllCount_BL_tech[,1]

# build DGEList
group <- c(1,1,1,1,2,2,2)
y <- DGEList(counts = cm, group=c(1,1,1,1,2,2,2))
str(y)
dim(y)

# paramenter for filtering low expressed genes
min.cpm <- 2
n.min.cpm <- 3
keep <- rowSums(cpm(y)>min.cpm) >= n.min.cpm
table(keep)
y <- y[keep,]
dim(y)

y$samples$lib.size <- colSums(y$counts)



# TMM normalization
y <- calcNormFactors(y, method="TMM")
y$samples

# prepare for edgeR glm
PRT  <- factor(c("YFP_NEG", "YFP_NEG", "YFP_NEG","YFP_NEG","YFP_POS", "YFP_POS", "YFP_POS"), levels= c("YFP_NEG", "YFP_POS"))
sample.names <- c("YFP_NEG_0", "YFP_NEG_1", "YFP_NEG_2", "YFP_NEG_3", "YFP_POS_0", "YFP_POS_1", "YFP_POS_2")
targets <- as.data.frame(cbind(sample.names,PRT))
design <- model.matrix(~0+PRT)

my.contrasts <- makeContrasts(
        C1 = (PRTYFP_POS-PRTYFP_NEG), levels=design  
)
interesting.contrasts <- c("C1")

# variance estimate
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

# QC plots
#plotMDS(y) # MDS plot
#plotMeanVar(y) # Mean vs Variance
#plotBCV(y) #

# EdgeR analysis & report
fit <- glmFit(y,design)
fit

# fdr threshold
fdr.t <- 0.05  # get all the expressed genes, then u can filter based on FDR for DEGs
project.name <- "RNAseq_CMM_Blanes_unique"
cat("experiment: ", project.name, "\n")
cat("thresholds: ", min.cpm, " ", n.min.cpm,"\n")
for (my.contrast in interesting.contrasts) {
        lrt <- glmLRT(fit, contrast=my.contrasts[,my.contrast])
        etable <- topTags(lrt, n=nrow(lrt$table), adjust.method="BH")
        etable <- etable$table[etable$table$FDR<fdr.t,]
        etable <- etable[ order(etable$FDR), ]  
        cat(my.contrast," ", dim(etable)[1], " genes\n")
        # write result
        # write.table( etable[,-c(2,3,4)], file=paste(project.name,my.contrast,"cpm", min.cpm,"nsamples",n.min.cpm,sep="."), row.names=FALSE)
}



################################################
#RPKM CALCULTING (Qiwen TAIR10.gene.length) #### 
################################################

# calculate gene lenght 
len = read.table("tair10.whole.genelength.txt",sep="\t")
names(len) = c("ID","length")
htseqAllCount_BL_tech = merge(htseqAllCount_BL_tech,len,by=c(1))
x = htseqAllCount_BL_tech[,c(2:9)] 
rownames(x) = htseqAllCount_BL_tech[,1]

rpkm = rpkm(x,htseqAllCount_BL_tech[,9],normalized.lib.sizes=TRUE)
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

str(rpkm)


####RPKM for all expressed genes####
rpkm_exp = rpkm[rpkm$name %in% rownames(y$count),] 

#write.table(x = rpkm_exp,"YFPs_RPKM_17Kgenes")

### Merge Table ### 

# merge p-values/FDR table with RPKM table
etable$name = rownames(etable)
edgeR_table = merge(rpkm,etable,by=c("name"))
str(table)

#write.table(edgeR_table,file="edgeR_table")

#######################################################################################
# Get non-coding RNAs (lncRNAs)
#######################################################################################

names(gene_anno) = c("name","anno")
all_rpkm = merge(edgeR_table,gene_anno,by=c("name")) #8858 genes DEGs+lncRNA
prot.cod = all_rpkm[all_rpkm$anno=="protein_coding_gene",1] #character 8716 genes only

#etable=list with 8858 (rowanmes 1st col names last col)
lncRNAs_plus_prot=rownames(etable) #character with ATG only 8858 genes

lncRNAs_only=setdiff(lncRNAs_plus_prot,prot.cod) #83 
lnc_rpkm = all_rpkm[all_rpkm$name %in% lncRNAs_only,]
lnc_rpkm

write.table(x=lnc_rpkm, "YFPs.Sig.FDR<0.05.lnc_rpkm")
