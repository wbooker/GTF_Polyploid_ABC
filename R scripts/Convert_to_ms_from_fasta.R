library(seqinr)
library(ape)
setwd("/home/wbooker/Dropbox/Documents/Research/ChrysoscelisVersicolor/fastaFiles_fulldataset/") #### Folder where fasta files are
#oldLoci <- read.table("EC_AE_locusInfo.txt")[,1] ##########uncomment this and two lines down if you want locus reference
files_list <- list.files(pattern = "*._L") #### Change pattern if in a different format
#files_list <- files_list[which(files_list %in% oldLoci)]
files <- files_list
files <- sample(files_list,244) ############ shuffle loci to randomize 
popsList <- read.csv("/home/wbooker/Dropbox/Documents/Research/ChrysoscelisVersicolor/SnpDocs/pops_all_MinMax.csv", header = F) #### csv file with information on populations and which to include
popsList <- popsList[which(popsList[,3]==1),] #### third column of the above csv file is set to 1 to include those individuals
popA <- c("EC")  ### diploid population name from popslist column 2
popB <- c("AE") ### tetraploid population name from popslist column 2
outgroup = ("AVI") ### Outgroup population name, one individual will be chosen randomly to be the reference at each locus
popBMin<- paste(c(popB, "_MIN"), collapse = "") #### subgenomes were split for our data into MIN and MAX sequences and had those suffixes for each individual, change if yours are different
popBMax<-  paste(c(popB, "_MAX"), collapse = "")
nucleotides <- c("a","c","t","g")
clockRate_avg <- 0.000000000876
clockRate_sd <- 0.0000000001
n0 = 100000
numLoci <- 50 #### Number of loci to include 
popSNP <- NULL
lociInfoMat <- c()
bpFile <- NULL
msMat <- NULL
D3 <- NULL
msMat <- rbind("msnsam tbs 20 -t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -n 1 tbs -n 2 tbs -n 3 tbs -g 2 tbs -g 3 tbs -ej tbs 3 2 -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs", "3579 27011 59243", "")
includedLoci = 0
for(locus in files){
  if(includedLoci < numLoci){
    alignment <- read.fasta(file=locus, as.string = FALSE, set.attributes = FALSE)
    popAList <- which(sapply(strsplit(names(alignment), "_"), '[',5) %in% sapply(strsplit(as.character(popsList[popsList[,2]%in%popA,1]) , "_"), '[', 5))
    popB_Min_List <- which(sapply(strsplit(names(alignment), "_"), '[',5) %in% sapply(strsplit(as.character(popsList[popsList[,2]%in%popBMin,1]) , "_"), '[', 5))
    popB_Min_List <- popB_Min_List[which(sapply(strsplit(names(alignment[popB_Min_List]), "_"), '[',7)=="min")]
    popB_Max_List <- which(sapply(strsplit(names(alignment), "_"), '[',5) %in% sapply(strsplit(as.character(popsList[popsList[,2]%in%popBMax,1]) , "_"), '[', 5))
    popB_Max_List <- popB_Max_List[which(sapply(strsplit(names(alignment[popB_Max_List]), "_"), '[',7)=="max")]
    popBList <- c(popB_Min_List, popB_Max_List)
    outInd <- sample(which(names(alignment) %in% popsList[popsList[,2]%in%outgroup,1]), 1)

    locusLength <- as.numeric(length(as.vector(alignment[[1]])))
    lociInfoMat <- rbind(lociInfoMat, c(locus, locusLength, length(popAList), length(popBList)/2))
    ind_Gmin <- NULL
    ind_Gmax <- NULL
    ind_minDivxy <- NULL
    ind_maxDivxy <- NULL
    popSNPtmp <- NULL
    popAseq <- NULL
    popBseq <- NULL
    popB_min_seq <- NULL
    popB_max_seq <- NULL
    
    for(i in popAList){
      alignment[i][[1]][which(alignment[i][[1]]=="-")] <- "N"
      popA_ind_seq <- paste(alignment[i][[1]], collapse = "")
      names(popA_ind_seq) <- names(alignment[i])
      popAseq <- c(popAseq,popA_ind_seq)
      #popSNPtmp <- rbind(popSNPtmp, alignment[i][[1]])
    }
    for(j in popBList){
      alignment[j][[1]][which(alignment[j][[1]]=="-")] <- "N"
      popB_ind_seq <- paste(alignment[j][[1]], collapse = "")
      names(popB_ind_seq) <- names(alignment[j])
      popBseq <- c(popBseq,popB_ind_seq)
    }
    for(j in popB_Min_List){
      alignment[j][[1]][which(alignment[j][[1]]=="-")] <- "N"
      popB_min_ind_seq <- paste(alignment[j][[1]], collapse = "")
      names(popB_min_ind_seq) <- names(alignment[j])
      popB_min_seq <- c(popB_min_seq,popB_min_ind_seq)
    }
    for(j in popB_Max_List){
      alignment[j][[1]][which(alignment[j][[1]]=="-")] <- "N"
      popB_max_ind_seq <- paste(alignment[j][[1]], collapse = "")
      names(popB_max_ind_seq) <- names(alignment[j])
      popB_max_seq <- c(popB_max_seq,popB_max_ind_seq)
    }
    
    
    pre_alignment <- as.DNAbin.DNAStringSet(c(popAseq,popB_min_seq,popB_max_seq))
    snps <- seg.sites(pre_alignment)
    if(length(snps) > 0){
      print(locus)
      popAseq <- NULL
      popBseq <- NULL
      popB_min_seq <- NULL
      popB_max_seq <- NULL
      popAseq_ms <- NULL
      popBseq_ms <- NULL
      popB_min_seq_ms <- NULL
      popB_max_seq_ms <- NULL
      
      for(i in popAList){
        alignment[i][[1]][which(alignment[i][[1]]=="-")] <- "N"
        popA_ind_seq <- paste(alignment[i][[1]][snps], collapse = "")
        names(popA_ind_seq) <- names(alignment[i])
        popAseq <- c(popAseq,popA_ind_seq)
        popA_ind_seq_ms <- alignment[i][[1]][snps]
        popAseq_ms <- rbind(popAseq_ms,popA_ind_seq_ms)
      }
      for(j in popBList){
        alignment[j][[1]][which(alignment[j][[1]]=="-")] <- "N"
        popB_ind_seq <- paste(alignment[j][[1]][snps], collapse = "")
        names(popB_ind_seq) <- names(alignment[j])
        popBseq <- c(popBseq,popB_ind_seq)
        popB_ind_seq_ms <- alignment[j][[1]][snps]
        popBseq_ms <- rbind(popBseq_ms,popB_ind_seq_ms)
      }
      for(j in popB_Min_List){
        alignment[j][[1]][which(alignment[j][[1]]=="-")] <- "N"
        popB_min_ind_seq <- paste(alignment[j][[1]][snps], collapse = "")
        names(popB_min_ind_seq) <- names(alignment[j])
        popB_min_seq <- c(popB_min_seq,popB_min_ind_seq)
        popB_min_ind_seq_ms <- alignment[j][[1]][snps]
        popB_min_seq_ms <- rbind(popB_min_seq_ms,popB_min_ind_seq_ms)
      }
      for(j in popB_Max_List){
        alignment[j][[1]][which(alignment[j][[1]]=="-")] <- "N"
        popB_max_ind_seq <- paste(alignment[j][[1]][snps], collapse = "")
        names(popB_max_ind_seq) <- names(alignment[j])
        popB_max_seq <- c(popB_max_seq,popB_max_ind_seq)
        popB_max_ind_seq_ms <- alignment[j][[1]][snps]
        popB_max_seq_ms <- rbind(popB_max_seq_ms,popB_max_ind_seq_ms)
      }
      
      alignment[outInd][[1]][which(alignment[outInd][[1]]=="-")] <- "N"
      outInd_seq_ms <- alignment[outInd][[1]][snps]
      
      ms_alignment <- rbind(popAseq_ms,popB_max_seq_ms,popB_min_seq_ms)
    
      remCol <- NULL
      for(i in 2:nrow(ms_alignment)){
        for(j in 1:ncol(ms_alignment)){
          if(as.character(ms_alignment[i,j]) == "N" | as.character(ms_alignment[i,j]) == "-"){
            remCol <- c(remCol,j)
            }
        }
      }
      if(length(remCol) != 0){
        ms_alignment <- ms_alignment[,-remCol]
        outInd_seq_ms <- outInd_seq_ms[-remCol]
      }
      if(!is.null(nrow(ms_alignment)) | !is.null(ncol(ms_alignment)) & locusLength >= 500){
        includedLoci = includedLoci + 1 
        clock <- strtrim(as.character(4*n0*locusLength*rnorm(1,clockRate_avg, clockRate_sd)),5)
        bpFile <- cbind(bpFile, c(locusLength,length(popAList), length(popBList)/2,clock,clock))
        ref_alignment <- outInd_seq_ms
        msMat <- rbind(msMat, "// 	54	1.17074	1.17074	1737	28	13	13	0.06531	8.47518	0.06531	8.47518	0.00000	0.00000	0.13154	1.63965	1.63965	60.63307	60.63307	0.18660	0.18660	1.49007	0.34644	0.34644	2.0751", paste(c("segsites: ",length(ref_alignment)), collapse = ""),"positions: ")
        for(i in 1:nrow(ms_alignment)){
          biallelic_alignment <- NULL
          for(j in 1:ncol(ms_alignment)){
            if((as.character(ms_alignment[i,j]) != "N") & (as.character(ms_alignment[i,j]) != as.character(ref_alignment[j]))){
              biallelic_alignment <- c(biallelic_alignment,1)
            }
            else{
              biallelic_alignment <- c(biallelic_alignment,0)
            }
          }
          biallelic_alignment <- paste(biallelic_alignment, collapse = "")  
          msMat <- rbind(msMat, biallelic_alignment)
        }
        msMat <- rbind(msMat, "")
      }
    }
  }
}
bpFile[1:3,] <- as.character(as.integer(bpFile[1:3,]))
bpFile <- rbind(c("#Sp1=chrys Sp2=vers", rep("",ncol(bpFile)-1)),bpFile)
write.table(msMat, file = paste(c(popA,"_",popB,"_msMat.txt"), collapse = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(lociInfoMat, file = paste(c(popA,"_",popB,"_locusInfo.txt"), collapse = ""), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(bpFile, file =paste(c("bpFile_",popA,"_",popB,".txt"), collapse = ""), col.names = F, row.names = F, quote = F, sep = "\t")
