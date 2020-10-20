get_overrepresented_groups<-function(M,groups,fold=NULL){
  # M: Matrix or data.frame with the data to analyze. Groups will based on rows, analyzing the data on columns
  # groups: Vector of length nrow(M), assigning a group to each row (sample)
  # fold: If used, it gets as positive those cases when data is over a specific fold-change on each group.
  #       It should be higher than 1.5, so data is at least 50% over-represented over the other groups.
  G<-unique(groups)
  U<-vector("numeric",length=ncol(M))
  names(U)<-colnames(M)
  OR<-vector("numeric",length=ncol(M))
  names(OR)<-colnames(M)
  
  for(g in 1:length(G)){
    gM<-M[groups==G[g],] # Data matrix for the specific group
    rM<-M[!groups==G[g],] # Data matrix for the rest of data
    for(c in 1:ncol(M)){
      gS<-sum(gM[,c])
      rS<-sum(rM[,c])
      if(gS>0 & rS==0){ 
        U[colnames(M)[c]]<-g 
      }else{
        if(!is.null(fold)){
          
          if(fold<=1.5){stop("fold should be >= 1.5, when it is used")}
          
          if(gS/rS>=fold){
            OR[colnames(M)[c]]<-g
          }
        }
      }
    }
  }
  if(is.null(fold)){
    return(U)
  }else{
    UOR<-rbind(U,OR)
    or<-paste("Over",fold,"fold",sep=" ")
    rownames(UOR)<-c("Unique",or)
    return(UOR)
  }
}

