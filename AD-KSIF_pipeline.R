##--------------------load the function files--------------------##
setwd("<PRINT_YOUR_PATH_TO_FUNCTIONS_HERE>")
source('PRE_trt.R')
source('TRE_crt.R')
source('PARA_cal.R')

##--------------------data pre-treatment--------------------##
library('stringr')
setwd("<PRINT_YOUR_PATH_TO_INPUT_FILE_HERE>")
#read the sequence file in FASTA format
otu_seq <- read.table("<PRINT_NAME_OF_INPUT_FASTA_FILE_HERE>", header = F, sep = "\n")
#sequence data pre-treatment
otuseq<-fasta_treatment(otu_seq)

##--------------------phylogenetic tree construction--------------------##
library(muscle)
otu_seq2<-readDNAStringSet("<PRINT_NAME_OF_INPUT_FASTA_FILE_HERE>", format="fasta",
                           nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
tree<-Tree_Construction(otu_seq2)
write.tree(tree,"otu_phylo.tre")

##--------------------binning--------------------##
#read the OTU abundance table
otu_abund <- read.table("<PRINT_NAME_OF_OTU_TABLE_FILE_HERE>", header = T, sep = ",")
#read the phylogenetic tree file in text format
otu_tre <- read.table("otu_phylo.tre", header = F, sep = ";")

#binning
binning.result<-binning(0.2,24,otu_tre)
binning.info<-allocation(binning.result,otuseq,otu_abund)
write.csv(binning.info,"binning.result.csv",row.names = F)

#output the sequence, abundance and phylogenetic tree of each bin
setwd("<PRINT_YOUR_PATH_TO_BIN_FILE_HERE>")
bin.num=length(unique(binning.info$Bin))#total number of the bins

for(i in 1:bin.num){
  #OTU abundance extraction of each bin
  bin<-binning.info[(binning.info$Bin==i),]
  abund<-subset(bin, select = -c(Bin, Seq))
  write.table(abund, paste("bin",as.character(i),".txt",sep=""), quote = F,row.names=FALSE,sep="\t")
  #OTU sequence extraction of each bin
  seq.ori<-subset(bin, select = c(OTU.ID, Seq))
  for(j in 1:nrow(seq.ori)){
    seq.ori$OTU.ID[j]<-paste('>',seq.ori$OTU.ID[j])
  }
  k=1
  j=1
  seq<-data.frame(matrix(NA,nrow = 2*nrow(seq.ori),1))
  while(k<=nrow(seq)){
    seq$V1[k]<-seq.ori$OTU.ID[j]
    k<-k+1
    seq$V1[k]<-seq.ori$Seq[j]
    k<-k+1
    j<-j+1
  }
  write.table(seq$V1,paste("seq_bin",as.character(i),".fasta",sep=""), quote = F,
              row.names=FALSE, col.name=FALSE)
  #phylogenetic tree construction of each bin
  seq2<-readDNAStringSet(paste("seq_bin",as.character(i),".fasta",sep=""), format="fasta",
                             nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
  tree.bin<-Tree_Construction(seq2)
  write.tree(tree.bin,paste("bin",as.character(i),"_align.phylip.tre",sep=""))
}

##--------------------index calculation--------------------##
library(ape)
library(dplyr)
#set number of iterations
n=999
#create matrix for index storation
EDindexes<-matrix(0,ncol=2,nrow=bin.num)
StaFunc<-matrix(0,ncol=2,nrow=bin.num)

#Loop through all bins files and calculate parameters
for(i in 1:bin.num){
  #read the OTU table
  otu = read.table(paste("bin",as.character(i),".txt",sep=""), header=T , row.names=1) 
  #read the phylogenetic tree
  phylo=read.tree(paste("bin",as.character(i),"_align.phylip.tre",sep=""))
  #extract the target samples by ID
  otu2<-otu[,grepl("<TARGET_SAMPLES_ID>", colnames(otu))]
  rownames(otu2)<-rownames(otu)
  #sample filtering
  otu3<-matrix(rep(1:nrow(otu2)),nrow=nrow(otu2),ncol=1)
  c="code"
  for (k in 1:ncol(otu2)){
    if (colSums(otu2)[k]>5){
      c<-cbind(c,colnames(otu2[k]))
    }
  }
  c<-c[-1]
  otu2<-subset(otu2,select=c)
  if (ncol(otu2)<2){
    next
  }
  #beta-NTI calculation
  bNTI=beta_NTI(otu2,phylo,n)
  write.csv(bNTI,paste("<PRINT_YOUR_PATH_TO_INDEX_FILE_HERE>/bNTI_",as.character(i),".csv"),quote=F)
  #RC-Bray-Curtis calculation
  RC=RC_BC(otu2,n)
  write.csv(RC,paste("<PRINT_YOUR_PATH_TO_INDEX_FILE_HERE>/RC_",as.character(i),".csv"),quote=F)
  bNTI[is.na(bNTI)]<-0
  RC[is.na(RC)]<-0
  #evolutionary & dispersal index calculation
  XY=X_Y(otu2,bNTI,RC)
  EDindexes[i,]<-XY
  #stability & category index calculation
  XY[is.na(XY)]<-0
  rhotheta=Rho_Theta(XY)
  StaFunc[i,]<-rhotheta
}

##--------------------result output--------------------##
result<-cbind(EDindexes,StaFunc)
colnames(result)<-c("X","Y","Rho","Theta")
write.csv(result,"result.csv",quote=F)