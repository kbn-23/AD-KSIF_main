Tree_Construction<-function(df){
  library(ape)
  #multi-sequence alignment
  aln <- muscle::muscle(df)
  #sequence trimming
  auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
  #alignment file output
  DNAStr = as(auto, "DNAStringSet")
  writeXStringSet(DNAStr, file="alignment.fasta")
  #integration of alignment results
  otu_seq3 <- read.table("alignment.fasta", header = F, sep = "\n")
  otu_seq3$V1[2]
  #obtain the fragments
  for (i in 1:nrow(otu_seq3)){
    if(grepl('>',otu_seq3$V1[i])){
      i<-i+1
      while(!grepl('>',otu_seq3$V1[i])){
        i<-i+1
      }
      cycle<-i-1
      break  
    }
  }
  #sequence assembly
  if(cycle > 2){
    for (i in 1:nrow(otu_seq3)){
      if(i %% cycle==2)
        blank<-i
      if((i %% cycle!=1) & (i %% cycle!=2)){
        otu_seq3$V1[blank]<-paste(otu_seq3$V1[blank],otu_seq3$V1[i],sep="")
        otu_seq3$V1[i]<-0
      }
    }
    otu_seq4 <- otu_seq3$V1[!(otu_seq3$V1==0)]
  }else{
    otu_seq4 <- otu_seq3$V1
  }
  write.table(otu_seq4,"otu_seq_alignment.fasta",row.names=FALSE,col.names=FALSE,quote=F)
  #phylogenetic tree construction by NJ method
  df<-read.dna("otu_seq_alignment.fasta", format = "fasta")
  df.distance <- dist.dna(df, model = "K80")
  tree <- nj(df.distance)
  return(tree)
}