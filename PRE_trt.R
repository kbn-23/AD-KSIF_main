fasta_treatment<-function(df){
  #divide the ID and the sequence in the FASTA file in to two rows
  otuseq <- data.frame(OTU = c(), Seq = c())
  tmp <-data.frame(OTU = c("OTU"), Seq = c("SEQ"))
  for (i in 1:nrow(df)) {
    if (i %% 2 == 1) {
      ch <- substr(df$V1[i], 2, str_length(df$V1[i]))
      tmp$OTU[1] <- ch
    } 
    else {
      tmp$Seq[1] <- df$V1[i]
      otuseq <- rbind(otuseq, tmp)
    }
  }
  return(otuseq)
}





binning<-function(d.max,n.min,input){
  tre <- input$V1
  
  df <- data.frame(Name = c("1"), #node
                   Dist = c(0), #The distance from the current node to the parent node
                   Lth = c(0.00), #The distance from the current node to the root node
                   Prt = c(0), #The ID of the parent node of the current node
                   nNode = c(0), #The number of leaf nodes owned by the current node
                   nChd = c(0), #The number of child nodes of the current node
                   Chd1 = c(0), #The ID of the first child node of the current node
                   Chd2 = c(0), #The number of the second child node of the current node
                   Chd3 = c(0) #The number of the third child node of the current node, if any
  )
  
  #Convert a one-dimensional phylogenetic tree to a data frame
  i <- 2 #Current string character position
  pst <- 1 #The number of the current presenting node
  n <- 1 #Number of nodes in the tree
  
  while(i < str_length(tre)) { #Process the one-dimensional phylogenetic tree tre character by character
    
    ch <- substr(tre, i, i) #Extract the i-th character
    
    if (ch == '(') { #Construct new non-leaf nodes (all leaf nodes are OTUs)
      n = n+1 #Total number of nodes +1
      tdf <- data.frame(Name = c(as.character(n)), 
                        Dist = c(0), 
                        Lth = c(0.00),
                        Prt = c(pst),
                        nNode = c(0),
                        nChd = c(0), 
                        Chd1 = c(0), 
                        Chd2 = c(0), 
                        Chd3 = c(0)
      ) #Create a new non-leaf node tdf, the parent node of tdf is the current node pst
      df <- rbind(df,tdf) #Add tdf to the df node list
      df$nChd[pst] <- df$nChd[pst] + 1 #Total number of child nodes of pst + 1
      if (df$nChd[pst] == 1){ 
        df$Chd1[pst] = n
      } else if (df$nChd[pst] == 2){
        df$Chd2[pst] = n
      } else {
        df$Chd3[pst] = n
      } #Add the position of tdf in df to the child node position corresponding to pst
      pst <- n #Current node update
      i <- i+1
    } else if (ch == ",") { #Returns the parent node
      pst <- df$Prt[pst]
      i <- i+1
    } else if (ch == ")") { #Remove the boostrap value and return to the parent node
      repeat {
        i <- i+1
        if (substr(tre, i, i) == ":"){
          break
        }
      } #Skip boostrap value
      pst <- df$Prt[pst] #Returns the parent node
    } else if (ch == ":") { #Processing distance data
      pt <- i+1 #Record the starting position of the distance data in the string
      while (substr(tre, i, i) != "," && substr(tre, i, i) != ")") {
        i <- i+1
      } #Find the end position of the distance data in the string backward
      df$Dist[pst] <- as.numeric(substr(tre, pt, i-1))
    } else { #Create a new leaf node
      pt <- i #Record the starting position of the leaf node name in the string
      while (substr(tre, i, i) != ":"){ #Look backwards for the end position of the leaf node name in the string
        i <- i+1
      }
      n = n+1 #Total number of nodes +1
      tdf <- data.frame(Name = c(substr(tre, pt, i-1)), 
                        Dist = c(0), 
                        Lth = c(0.00),
                        Prt = c(pst),
                        nNode = c(1),
                        nChd = c(0), 
                        Chd1 = c(0), 
                        Chd2 = c(0), 
                        Chd3 = c(0)
      ) #Create a new leaf node tdf, the parent node of tdf is the current node pst
      df <- rbind(df,tdf) #Add tdf to the df node list
      df$nChd[pst] <- df$nChd[pst] + 1 #Total number of child nodes of pst + 1
      if (df$nChd[pst] == 1){ 
        df$Chd1[pst] = n
      } else if (df$nChd[pst] == 2){
        df$Chd2[pst] = n
      } else {
        df$Chd3[pst] = n
      } #Add the position of tdf in df to the child node position corresponding to pst
      pst <- n #Current node update
    }
  }
  
  #Create a traversal stack S
  S <- c(1) #S is used as a stack to store the tree structure of the breadth-first search
  p <- 1 #Points to the current node in the stack
  q <- 2 #Points to the empty space at the top of the stack
  while (p < q){ #Until there are no untraversed nodes
    if (df$Chd1[S[p]] != 0){ #If there are child nodes, they are stored in the stack
      S[q] <- df$Chd1[S[p]]
      q <- q+1
    }
    if (df$Chd2[S[p]] != 0){ #If there are child nodes, they are stored in the stack
      S[q] <- df$Chd2[S[p]]
      q <- q+1
    }
    if (df$Chd3[S[p]] != 0){ #If there are child nodes, they are stored in the stack
      S[q] <- df$Chd3[S[p]]
      q <- q+1
    }
    p <- p+1 #The current node has been processed, accessing unprocessed nodes
  }
  
  #Count the number of OTUs contained in each node
  while (q > 1) { #Count the number of leaf nodes contained in each node, and q points to the empty site/processed site of the stack
    q <- q-1
    if (df$Chd1[S[q]] != 0){ #If there is a leaf node in the child node, it is counted in the total number of leaf nodes
      df$nNode[S[q]] <- df$nNode[S[q]]+df$nNode[df$Chd1[S[q]]]
    }
    if (df$Chd2[S[q]] != 0){ #If there is a leaf node in the child node, it is counted in the total number of leaf nodes
      df$nNode[S[q]] <- df$nNode[S[q]]+df$nNode[df$Chd2[S[q]]]
    }
    if (df$Chd3[S[q]] != 0){ #If there is a leaf node in the child node, it is counted in the total number of leaf nodes
      df$nNode[S[q]] <- df$nNode[S[q]]+df$nNode[df$Chd3[S[q]]]
    }
  }
  
  #Calculate the distance from each node to the root node
  for (i in 2:nrow(df)) { #Calculate the distance from each node to the root node
    df$Lth[i] <- df$Lth[df$Prt[i]] + df$Dist[i]
  }
  
  #Calculate the maximum value of each node derived from the leaf node to the root node
  q <- n+1 #Points to the empty space at the top of the stack
  lf2rt <- rep(0,n) #The maximum value from the leaf node to the root node of the i-th node
  while (q > 1) { #Count the maximum distance from the leaf nodes to the root node contained in each node, and q points to the empty site/processed site of the stack
    q <- q-1
    if (df$nNode[S[q]] == 1) {
      lf2rt[S[q]] <- df$Lth[S[q]]
    } else {
      if (df$Chd1[S[q]] != 0){ #If there is a leaf node in the child node, it is counted in the total number of leaf nodes
        lf2rt[S[q]] <- max(lf2rt[S[q]], lf2rt[df$Chd1[S[q]]])
      }
      if (df$Chd2[S[q]] != 0){ #If there is a leaf node in the child node, it is counted in the total number of leaf nodes
        lf2rt[S[q]] <- max(lf2rt[S[q]], lf2rt[df$Chd2[S[q]]])
      }
      if (df$Chd3[S[q]] != 0){ #If there is a leaf node in the child node, it is counted in the total number of leaf nodes
        lf2rt[S[q]] <- max(lf2rt[S[q]], lf2rt[df$Chd3[S[q]]])
      }
    }
  }
  
  #Calculate the maximum distance between leaf nodes within each node
  ds <- rep(0,n) #The maximum distance between leaf nodes within the i-th node
  q <- n+1
  while (q > 1) { #Count the maximum distance from the leaf nodes to the root node contained in each node, and q points to the empty site/processed site of the stack
    q <- q-1
    if (df$nNode[S[q]] == 1) {
      ds[S[q]] <- 0
    } else {
      t <- c(0,0,0) #Stores the maximum distance from leaf to root
      ts <- c(0,0,0) #Maximum distance between leaf nodes within a child node
      if (df$Chd1[S[q]] != 0){ #If the child node exists, store the maximum distance from the leaf node to which the child node belongs to the root node
        t[1] <- lf2rt[df$Chd1[S[q]]]
        ts[1] <- ds[df$Chd1[S[q]]]
      }
      if (df$Chd2[S[q]] != 0){ #If the child node exists, store the maximum distance from the leaf node to which the child node belongs to the root node
        t[2] <- lf2rt[df$Chd2[S[q]]]
        ts[2] <- ds[df$Chd2[S[q]]]
      }
      if (df$Chd3[S[q]] != 0){ #If the child node exists, store the maximum distance from the leaf node to which the child node belongs to the root node
        t[3] <- lf2rt[df$Chd3[S[q]]]
        ts[3] <- ds[df$Chd3[S[q]]]
      }
      ds[S[q]] <- max(sum(t)-min(t)-2*df$Lth[S[q]], ts)
    }
  }
  
  #Determine the phylogenetic distance L, ensuring that after L is cut off, the distance between all connected species is less than the threshold
  dL <- 0
  for (i in 1:n) {
    if (ds[i] > 0.2 && dL < df$Lth[i]){
      dL <- df$Lth[i]
    }
  }
  
  #Binning
  mark_bin <- function(x, bp, p){ #The bin number of the xth node is bp, and the representative node is p
    if(x != 0){
      bin[x] <<- bp
      hBin[x] <<- p
      mark_bin(df$Chd1[x], bp, p)
      mark_bin(df$Chd2[x], bp, p)
      mark_bin(df$Chd3[x], bp, p)
    }
  }
  
  bin <- rep(0,n) #The bin number of the i-th OTU
  hBin <- rep(0,n) #The root of the bin corresponding to the i-th node (hBin[i]==i)
  nBin <- rep(0,n) #The number of species in the i-th bin
  bp <- 1 #Points to unassigned bin numbers
  for (i in 1:n) {
    if (df$Lth[i] > dL && bin[i] == 0){
      nBin[bp] <- df$nNode[i]
      hBin[i] <- i
      mark_bin(i,bp,i)
      bp <- bp+1
    } else if (df$nNode[i] == 1 && bin[i] == 0) {
      nBin[bp] <- df$nNode[i]
      hBin[i] <- i
      bin[i] <- bp
      bp <- bp+1
    }
  }
  
  #merging
  hash <- matrix(FALSE,n,n) #hash[i,j]==TRUE means i is the ancestor node of j
  for (i in 1:n) {
    p <- i
    while (p != 0) {
      hash[p,i] <- TRUE
      p <- df$Prt[p]
    }
  }
  
  NNdist <- function(i,j){
    p <- j #p represents an ancestor node of j
    while (!hash[p,i]) { #If node p is not an ancestor of node i
      p <- df$Prt[p] #p is updated to the ancestor node of p
    }
    return(df$Lth[i]+df$Lth[j]-2*df$Lth[p]) #Returns the distance between i and j
  }
  
  ds <- 24 #Merge until each bin contains at least ds leaf nodes
  flag <- TRUE
  while (flag) { #Flag determines whether there is still a bin with less leaf nodes than ds
    flag <- FALSE
    subbin <- unique(hBin[which(hBin != 0)]) #Extract the representative node of the current bin (hBin[i]==i)
    #  print(length(subbin)) #Debugging in loop
    for(i in rev(subbin)) {
      if (hBin[i] == i && nBin[bin[i]] < ds) { #If the i node is still the representative node and the number of OTUs in the bin is less than ds
        flag <- TRUE #This round of loop performs the merge operation
        dL <-lf2rt[1] * 2 #dL is the shortest distance from each representative node to the i node
        tp <- 0 #The representative node with the shortest distance to the i node
        for(j in subbin) { 
          if (hBin[j] == j && i != j) { #If the j node is still the representative node and is not the i node
            tL <- NNdist(i, j) #tL is the distance between i and j
          }
          if (tL < dL) { #If a representative node j with a shorter distance to i is found
            dL <- tL #dL replaced by tL
            tp <- j #tp is denoted by j
          }
        }
        nBin[bin[tp]] <- nBin[bin[tp]] + nBin[bin[i]] #The total number of OTUs in the bin where tp is located plus the number of OTUs in the bin where i is located
        nBin[bin[i]] <- 0 #The original bin of i is empty
        hBin[which(hBin == i)] <- tp #The nodes belonging to the i bin are changed to belong to the tp bin.
        bin[which(bin == bin[i])] <- bin[tp]
        #print(paste(i," bin to ",tp,sep = "")) #Output binning details
      }
    }
  }
  #obtain all the nodes and the associated root nodes in the raw data
  result <- data.frame(OTU = df$Name, Bin = hBin)
  #delete non-leaf nodes
  for (i in nrow(result):1) {
    if (!str_detect(result$OTU[i],"OTU")){
      result <- result[-i,]
    }
  }
  #obtain the code of binning root nodes
  outbin <- unique(result$Bin)
  #replace the binning root nodes as:
  for (i in 1:length(outbin)) {
    result$Bin[which(result$Bin == outbin[i])] <- i
  }
  #reorder the results
  result <- result[order(result$Bin),]
}





allocation<-function(otubin,otuseq,otuabund){
  #allocate each OTU into bins
  result <- cbind(Bin = rep(0,nrow(otuabund)), Seq = rep("",nrow(otuabund)), otuabund)
  
  for (i in 1:nrow(result)) {
    result$Bin[i] <- otubin$Bin[which(otubin$OTU == result$OTU.ID[i])]
    result$Seq[i] <- otuseq$Seq[which(otuseq$OTU == result$OTU.ID[i])]
  }
  
  result <- result[order(result$Bin),]
  return(result)
}