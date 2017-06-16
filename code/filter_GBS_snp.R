##########################################
### Script to filter GBS hapmap files ####
##########################################
# date: 2017-06-11
# singhdj2@ksu.edu
# *Input HapMap file processed: Tassel Raw SNP 'hmp.txt' >> join parts with bed (Liang Gao's python script) 
# Filter criteria: MAF>0.05; >5% Het; missing >80% per marker ; line >75% missing 

require(data.table)
require(dplyr)
require(tidyr)
#require(RMySQL)
#require(rrBLUP)

## Execute Liang's python scripts to:
##     i. join chrom position parts with assembly bed file
##     ii. convert alleles to genotype calls: AA,AB,BB,NA; adds a summary stats at the end
# Example usage on Beocat environment:
system("python2 ~/scripts/format.hmp.id.parts.to.chr.py <yourFILE>.hmp.txt > <yourFile>.chr.hmp.txt")
system("python2 ~/scripts/GBS_filter_TASSEL5_hmp.py <yourFile>.chr.hmp.txt > <yourFile>.chr.hmp.filt.MSTmap.txt")

hmp <- fread("gbs/projects/BISA14-17/BISA14-17.chr.hmp.filt.MSTmap.txt",sep = "\t",header = T)
class(hmp)
hmp <- as.data.frame(hmp)
snp <- hmp
snp[1:5,1:20]
# convert NA coding from "-" to "N"
snp[snp=='-'] <- "N"
length(unique(substr(snp$`rs#`,2,3))) #distinct chromosomes

## remove some undesired columns ##
snp <- snp %>% select(-chrom, -pos, -center, -strand, -contains("blank"), -contains("test"))
snp[1,1:15]         


## calculate minor allele frequencey and percent het ##
MAF = apply(cbind(snp$cntAA, snp$cntBB), 1, min)/ apply(cbind(snp$cntAA, snp$cntBB, snp$cntAB), 1, sum)
HET = snp$cntAB / apply(cbind(snp$cntAA, snp$cntBB, snp$cntAB), 1, sum)
snpNA = snp$cntNA / apply(cbind(snp$cntAA, snp$cntBB, snp$cntAB, snp$cntNA), 1, sum) # i.e. missing genos per marker
snp = cbind(snp[,c(1:4)], MAF, HET, snpNA, snp[,c(5:ncol(snp))])

hist(MAF)
hist(snp$HET)
hist(snp$cntAA)
hist(snp$cntNA)
snp[1:3,1:50]

##########################################
### filter with MAF, HET and NA per marker
snp <- snp[snp$snpNA < 0.20,]  # NA filter
snp <- snp[snp$MAF > 0.05,]  # MAF filter
snp <- snp[snp$HET < 0.10,]  # HET filter


## summary: data per line 
lineData = apply(snp!="N", 2, sum)[-c(1:14)] #data present per line
hist(lineData/nrow(snp), main="GS 2014-17 GID's", xlab="# SNPs", ylab="# lines", breaks=60)
snp1 <- snp[,-c(1:14)] #remove columns
na25 = lineData[lineData/nrow(snp1)>0.75]  ## get rid of lines with more than 75% missing data
snp1 <- snp1[,colnames(snp1) %in% names(na25)] # apply line filter
dim(snp1); dim(snp)
snp <- cbind(snp[,8],snp1) #combine with snp column with marker positions (number 8 here)
dim(snp)
snp[1:2,1:20]

colnames(snp)[1] <- "markerName"
## write to file
outfile = "~/bisa14-17_snps_20170611.txt"
write.table(snp, file = outfile, quote=FALSE)

#############################################
##### turn into -1,0,1 matrix for rrBLUP ####
hap = read.table(file= "~/bisa14-17_snps_20170611.txt", header=TRUE)
hap[1:5,1:5]
rownames(hap) <- hap$markerName
hap = hap[,-1]
#hap = read.table(file= "", header=TRUE)
hap[1:10,1:20]
hap01 = hap
hap01[,1:ncol(hap01)]=NA
# convert codes to numerical
hap01[hap=="A"]= 1
hap01[hap=="B"]= -1
hap01[hap=="X"]= 0
hap01[hap=="N"]= NA

hap01[1:3,1:3]

#write to file
write.table(hap01, file= "~/bisa14-17_snps_numeric_20170611.txt" , col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

# remove the extra data objects from workspace
rm(snp1)
rm(lineData)
rm(snpNA)
rm(snpData)
rm(het)
rm(HET)
rm(MAF)
rm(hap)



