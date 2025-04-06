Tree_Construction<-function(df){
  library(ape)
  #多序列比对
  aln <- muscle::muscle(df)
  #序列修剪
  auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
  #比对文件输出
  DNAStr = as(auto, "DNAStringSet")
  writeXStringSet(DNAStr, file="alignment.fasta")
  #比对文件序列整理
  otu_seq3 <- read.table("alignment.fasta", header = F, sep = "\n")
  otu_seq3$V1[2]
  #获得分段数
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
  #拼接完整序列
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
  #NJ法建立系统发育树并输出
  df<-read.dna("otu_seq_alignment.fasta", format = "fasta")
  df.distance <- dist.dna(df, model = "TS")
  tree <- nj(df.distance)
  return(tree)
}
