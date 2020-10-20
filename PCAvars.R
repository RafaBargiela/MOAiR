        
PCAvars<-function(pca,n.vars=nrow(pca$rotation),xlim=c(-1,1),ylim=c(-1,1),arrows=TRUE,
                  vars.lwd=par("lwd"),vars.cols="black",vars.labs.col="black",lab.x.corr=rep(0,nrow(pca$rotation)),lab.y.corr=rep(0,nrow(pca$rotation)),
                  get_pca_res=TRUE,labels=TRUE,plot=FALSE,add=FALSE, ...){
  # pca : object of class "prcomp"
  # n.vars: numeric, number of variables to draw in order of contribution.
  # vars.cols: character value or vector defining colors for variables arrows
  # vars.lwd: numeric, width for variable arrows line.
  # get_pca_res: if TRUE, returns pca with additional attributs for variables
  # Labels: Logical, if TRUE variables labels are printed
  
  if(grepl("prcomp",paste(class(pca),collapse=" "))==FALSE){
      stop(pca, "is not of class prcomp")
  }
  
  message("#####################################################################
          NOTE: If you are adding variables arrows to previous plot set add=TRUE
          ######################################################################")
  
  # Analysing the variables data
  pca$vars.contribution<-(pca$rotation^2)*100 # Percentage of contribution of each variable 
  get_coords_func<-function(loadings,sdev){loadings*sdev} # Function to help getting coords when processing loadings (rotation) matrix by column
  pca$vars.coords<-t(apply(pca$rotation,1,get_coords_func,pca$sdev)) # Scaled coordinates for variables (usaully PC1 and PC2 are between -1 and 1)
  pca$vars.cos2<-pca$vars.coords^2
  
  # Getting total contribution for first two components
  Cont2D<-vector("numeric")
  Cos2_2D<-vector("numeric")
  for(v in 1:nrow(pca$vars.contribution)){
    Cont2D[rownames(pca$vars.contribution)[v]]<-(sum(pca$vars.contribution[v,1:2]))*100/200
    Cos2_2D[rownames(pca$vars.contribution)[v]]<-sum(pca$vars.cos2[v,1:2])
  }
  Cont2D<-sort(Cont2D,decreasing=TRUE)
  pca$rotation<-pca$rotation[names(Cont2D),] # Sorting vars coords by contribution on first 2 PCs
  pca$vars.coords<-pca$vars.coords[names(Cont2D),] # Sorting vars coords by contribution on first 2 PCs
  pca$vars.cos2<-pca$vars.cos2[names(Cont2D),] # Sorting vars cos2 by contribution on first 2 PCs
  pca$vars.contribution<-pca$vars.contribution[names(Cont2D),] # Sorting vars contribution by contribution on first 2 PCs
  pca$vars.2D.contribution<-Cont2D # Relative contribution of variables to first 2 components
  pca$vars.2D.cos2<-Cos2_2D[names(Cont2D)] # Cos2 relative to the 2D representation
  
  if(plot==TRUE){
    if(add==TRUE){
      par(new=TRUE)
    }
    plot(NA,NA,axes=FALSE,ylab="",xlab="",xlim=xlim,ylim=ylim)
    if(n.vars==0){
      n<-nrow(pca$var.coords)
    }else{
      n<-n.vars
    }
    for(v in 1:n){
      if(arrows==TRUE){
        arrows(0,0,pca$vars.coords[v,1],pca$vars.coords[v,2],code=2,lwd=vars.lwd,xpd=TRUE,angle=30,length=0.05,
               col=ifelse(length(vars.cols)>1,vars.cols[v],vars.cols))
      }
      if(labels==TRUE){
        dx<-pca$vars.coords[v,1]
        dy<-pca$vars.coords[v,2]
        a<-atan2(dy,dx)*(180/pi)
        xl<-(pca$vars.coords[v,1]+0.075*cos(as_radians(a)))+lab.x.corr[v]
        yl<-(pca$vars.coords[v,2]+0.075*sin(as_radians(a)))+lab.y.corr[v]
        shadowtext(xl,yl,labels=rownames(pca$vars.coords)[v],xpd=TRUE,adj=c(0.5,0.5),
                   col=ifelse(length(vars.labs.col)>1,vars.labs.col[v],vars.labs.col),bg="white",r=0.1,...)
      }
    }
  }
  if(get_pca_res==TRUE){
    return(pca)
  }
}
