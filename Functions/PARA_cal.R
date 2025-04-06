beta_NTI <- function(otu,phylo,n){
  library(picante)
  #ensure the consistence of the ID in phylogenetic tree and the OTU table
  match.phylo.otu = match.phylo.data(phylo, otu)  
  #observed betaMNTD calculation, where abundance.weighted=T is for considering the OTU abundance (default = F).
  beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T))  
  #ID check
  #identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)) 
  #identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted))
  
  ##null model betaMNTD calculation
  #set repeats of randomization 
  beta.reps = n
  #null model generation
  rand.weighted.bMNTD.comp=array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps)); 
  for (rep in 1: beta.reps) { 
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T ,exclude.conspecifics = F));  
    print(c(date(),rep)); 
  } 
  ##βNTI calculation
  weighted.bNTI = matrix(0,nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data)); 
  for (columns in 1:(ncol(match.phylo.otu$data)-1)) { 
    for (rows in (columns+1): ncol(match.phylo.otu$data)) {
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,]; 
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals); 
      rm("rand.vals"); 
    }; 
  }; 
  rownames(weighted.bNTI) = colnames(match.phylo.otu$data); 
  colnames(weighted.bNTI) = colnames(match.phylo.otu$data); 
  return(weighted.bNTI)
}





RC_BC<-function(x,n){
  x = t(x)
  gamma<-ncol(x)
  #obtain 0-1 matirx for OTU existence
  x.inc<-ceiling(x / max(x))
  RC<-matrix(rep(0,nrow(x)*nrow(x)),nrow(x),nrow(x))
  dimnames(RC)<-list(rownames(x),rownames(x))
  #number of OTUs
  occur <- apply(x.inc, MARGIN = 2, FUN = sum)
  #distrubution of OTU abundance in sample
  abundance <- apply(x, MARGIN = 2, FUN = sum)
  library(foreach)
  library(doParallel)
  for (null.one in 1:(nrow(x) - 1)) {
    for (null.two in (null.one + 1):nrow(x)) {
      #calculate expectation of Bray-Curtis distance between null models
      null_bray_curtis <- NULL
      null_bray_curtis <- foreach (i = 1:n,.combine = "c") %dopar%{
        library(vegan)
        #generate null com 1
        com1<-rep(0,gamma)
        #existence condition
        com1[sample(1:gamma,
                    sum(x.inc[null.one, ]),
                    replace = FALSE,
                    prob = occur)] <- 1
        #assign abundance value
        com1.samp.sp = sample(which(com1 > 0),
                              (sum(x[null.one, ]) - sum(com1)),
                              replace = TRUE,
                              prob = abundance[which(com1 > 0)])
        com1.samp.sp = cbind(com1.samp.sp, 1)
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[, 2], com1.samp.sp[, 1], FUN = sum))
        colnames(com1.sp.counts) = 'counts'
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts))
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts
        #generate com2 by the same algorithm
        com2<-rep(0,gamma)
        com2[sample(1:gamma,
                    sum(x.inc[null.two, ]),
                    replace = FALSE,
                    prob = occur)] <- 1
        com2.samp.sp = sample(which(com2 > 0),
                              (sum(x[null.two, ]) - sum(com2)),
                              replace = TRUE,
                              prob = abundance[which(com2 > 0)])
        com2.samp.sp = cbind(com2.samp.sp, 1)
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[, 2], com2.samp.sp[, 1], FUN =
                                                sum))
        colnames(com2.sp.counts) = 'counts'
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts))
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts
        com2
        null.x<-rbind(com1,com2)
        null_bray_curtis = vegdist(null.x, method = 'bray')
      }
      #calculate observed Bray-Curtis distance between samples
      obs.bray = vegdist(x[c(null.one, null.two), ], method = 'bray')
      #extract the number of the observed results lower than null model expectation
      num_less_than_in_null = sum(null_bray_curtis < obs.bray)
      num_exact_matching_in_null = sum(null_bray_curtis == obs.bray)
      #RC bray-curtis calculation
      RC[null.two,null.one]=((
        num_less_than_in_null + (num_exact_matching_in_null) / 2
      ) / n)
      #normalization
      RC[null.two,null.one] = (RC[null.two,null.one] - .5) * 2
      print(c(date(),null.two))
    }
    print(c(date(),null.one))
    }
  return(RC)
}





X_Y<-function(otu,bNTI,RC){
  #relative abundance in bins
  relativeab.samp<-matrix(0,1,ncol(otu))
  relativeab.samp[1,]=colSums(otu)/sum(otu)
  abweighted_bNTI<-0
  total_weight<-0
  #evolution index calculation
  for (i in 1:(ncol(otu)-1)){
    for (j in (i+1):ncol(otu)){
      if(bNTI[j,i]!=0){
        abweighted_bNTI<-abweighted_bNTI+bNTI[j,i]*(relativeab.samp[1,i]+relativeab.samp[1,j])/2
        total_weight<-total_weight+(relativeab.samp[1,i]+relativeab.samp[1,j])/2
      }
    }
  }
  X<-abweighted_bNTI/total_weight
  #dispersal index calculation
  abweighted_RC<-0
  total_weight<-0
  for (i in 1:(ncol(otu)-1)){
    for (j in (i+1):ncol(otu)){
      if(RC[j,i]!=0){
        abweighted_RC<-abweighted_RC+RC[j,i]*(relativeab.samp[1,i]+relativeab.samp[1,j])/2
        total_weight<-total_weight+(relativeab.samp[1,i]+relativeab.samp[1,j])/2
      }
    }
  }
  Y<-abweighted_RC/total_weight
  # normalization
  Standard_X<-X/1.96
  Standard_Y<-Y/0.95
  X_Y<-rep(0,2)
  X_Y[1]=Standard_X
  X_Y[2]=Standard_Y
  return(X_Y)
}





Rho_Theta<-function(XY){
  #coordinate transformation
  X=XY[1]
  Y=XY[2]
  rho_theta<-rep(0,2)
  rho=sqrt(X^2+Y^2)
  if(X>=0&Y>=0){
    theta=atan(Y/X)/pi
  }else if(X>=0&Y<0){
    theta=atan(Y/X)/pi+2
  }else{
    theta=atan(Y/X)/pi+1
  }
  rho_theta[1]=rho
  rho_theta[2]=theta
  return(rho_theta)
}
