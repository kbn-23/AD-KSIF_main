beta_NTI <- function(otu,phylo,n){
  #加载这个数据包 
  library(picante)
  #确保系统发育树上的名称与 otu表中的名称排序一致 
  match.phylo.otu = match.phylo.data(phylo, otu)  
  #计算 betaMNTD，这里 abundance.weighted=T表示考虑物种丰度，如果不设定该参数，则默认 abundance.weighted=F，不考虑物种丰度 
  beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T))  
  #检查 beta.mntd.weighted矩阵中的列名是否和系统发育树以及 OTU表中一致，结果应该为 T 
  identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)) 
  #同上，只是一个检查，结果应该为 T 
  identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted))  
  ##计算 betaMNTD随机值 
  #随机化次数 
  beta.reps = n
  rand.weighted.bMNTD.comp=array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps)); 
  for (rep in 1: beta.reps) { 
    #这一步执行的随机化是基于对系统发育树中 OTU标签的随机重排进行的，如果沿用 1中的方法对每 一个群落进行随机构建后再计算 betaMNTD,则可以选择使用rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(randomizeMatrix(comun, null.model = c("frequency", "richness","independentswap", "trialswap"), iterations = 1000),cophenetic(match.phylo.comun$phy),abundance.weighted=T ,exclude.conspecifics = F))命令完成该步运算 
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T ,exclude.conspecifics = F));  
    print(c(date(),rep)); 
  } 
  ##计算 βNTI 
  weighted.bNTI = matrix(0,nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data)); 
  for (columns in 1:(ncol(match.phylo.otu$data)-1)) { 
    for (rows in (columns+1): ncol(match.phylo.otu$data)) { 
      #按照 βNTI计算公式计算每一对样本间的 βNTI值 
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
  #获得OTU存在性矩阵
  x.inc<-ceiling(x / max(x))
  RC<-matrix(rep(0,nrow(x)*nrow(x)),nrow(x),nrow(x))
  dimnames(RC)<-list(rownames(x),rownames(x))
  #获得样本物种数量
  occur <- apply(x.inc, MARGIN = 2, FUN = sum)
  #获得样本丰度分布
  abundance <- apply(x, MARGIN = 2, FUN = sum)
  library(foreach)
  library(doParallel)
  for (null.one in 1:(nrow(x) - 1)) {
    for (null.two in (null.one + 1):nrow(x)) {
      #获得样本间空模型BC距离期望
      null_bray_curtis <- NULL
      null_bray_curtis <- foreach (i = 1:n,.combine = "c") %dopar%{
        library(vegan)
        com1<-rep(0,gamma)
        #确定存在性条件
        com1[sample(1:gamma,
                    sum(x.inc[null.one, ]),
                    replace = FALSE,
                    prob = occur)] <- 1
        #赋予丰度值
        com1.samp.sp = sample(which(com1 > 0),
                              (sum(x[null.one, ]) - sum(com1)),
                              replace = TRUE,
                              prob = abundance[which(com1 > 0)])
        com1.samp.sp = cbind(com1.samp.sp, 1)
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[, 2], com1.samp.sp[, 1], FUN = sum))
        colnames(com1.sp.counts) = 'counts'
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts))
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts
        #按同样方法生成com2
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
      #获得实际群落样本间BC距离
      obs.bray = vegdist(x[c(null.one, null.two), ], method = 'bray')
      #预期值低于实际值的数量
      num_less_than_in_null = sum(null_bray_curtis < obs.bray)
      #令计算结果包含等于预期值的数量
      num_exact_matching_in_null = sum(null_bray_curtis == obs.bray)
      RC[null.two,null.one]=((
        num_less_than_in_null + (num_exact_matching_in_null) / 2
      ) / n)
      #标准化到[-1,1]之间
      RC[null.two,null.one] = (RC[null.two,null.one] - .5) * 2
      print(c(date(),null.two))
    }
    print(c(date(),null.one))
    }
  return(RC)
}





X_Y<-function(otu,bNTI,RC){
  #替换缺省值

  #子组内各样本的相对丰度
  relativeab.samp<-matrix(0,1,ncol(otu))
  relativeab.samp[1,]=colSums(otu)/sum(otu)
  abweighted_bNTI<-0
  total_weight<-0
  #计算进化指数
  for (i in 1:(ncol(otu)-1)){
    for (j in (i+1):ncol(otu)){
      if(bNTI[j,i]!=0){
        abweighted_bNTI<-abweighted_bNTI+bNTI[j,i]*(relativeab.samp[1,i]+relativeab.samp[1,j])/2
        total_weight<-total_weight+(relativeab.samp[1,i]+relativeab.samp[1,j])/2
      }
    }
  }
  X<-abweighted_bNTI/total_weight
  #计算扩散指数
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
  #标准化
  Standard_X<-X/1.96
  Standard_Y<-Y/0.95
  X_Y<-rep(0,2)
  X_Y[1]=Standard_X
  X_Y[2]=Standard_Y
  return(X_Y)
}





Rho_Theta<-function(XY){
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
