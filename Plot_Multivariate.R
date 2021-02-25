plot_Multivariate<-function(multivA,TB,xlim=c(-1,1),ylim=c(-1,1),pch=19,bg="white",ellipse.groups=NULL, points=TRUE,sites.labs=FALSE,ell.cols="grey50",ell.radius=1,ell.center.pch=19,plot=TRUE,
                            xlab="",ylab="",axes=TRUE,x.ticks=seq(xlim[1],xlim[2],1),y.ticks=seq(ylim[1],ylim[2],1), ...){
  # multivA : object class "prcomp" or "metaMDS", usually produced by prcomp or metaMDS funcitons from stats and Vegan packages, respectively.
  # TB: matrix or data.frame with the data used to produce the ordination analysis (PCA or NMDS) 
  #         If groups TRUE, it should contains a column assigning the each case to a group, which number is assigned to groups.col
  # ellipse.groups: vector containing the group assignation to each sample (row of the data matrix) used to draw ellipses around them, based on group variance
  # x.ticks and y.ticks: vector defining which and where axis labels are printed.
  # plot: Logical, if TRUE axis, axis labels and figure box are drawn. Default TRUE
  # points: Logical, if TRUE Sites dots are drawn. Default TRUE
  # sites.labs: Logical, if TRUE sites labels are displayed. Default FALSE
  require(car)
  if(plot==TRUE){
    plot(NA,NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,axes=FALSE, ... )
    if(axes==TRUE){
      axis(1,x.ticks,labels=TRUE,tick=TRUE,lwd.ticks=par("lwd"),font.axis=2,xpd=TRUE,las=1, ... )
      axis(2,y.ticks,labels=TRUE,tick=TRUE,lwd.ticks=par("lwd"),font.axis=2,xpd=TRUE,las=1, ... )
      box(lwd=2)
    }
    arrows(0,ylim[1],0,ylim[2],lty=3,lwd=1,col="grey50",code=0)
    arrows(xlim[1],0,xlim[2],0,lty=3,lwd=1,col="grey50",code=0)
  }
  # Analysing the individuals/cases/species data
  if(grepl("prcomp",paste(class(multivA),collapse=" "))==TRUE){
    sam<-multivA$x[,1:2]      
  }else{
    if(grepl("metaMDS",paste(class(multivA),collapse=" "))==TRUE){
      sam<-multivA$points[,1:2]        
    }else{
      stop("multivA object is not of class prcomp or meta|monoMDS")
    }
  }
  
  
  if(!is.null(ellipse.groups)){
    G<-unique(ellipse.groups)
    if(length(cols)==1){
      cols<-brewer.pal(length(G),"Set2") # No more than 8 groups
    }else{
      cols<-cols
    }
    names(cols)<-G
    
    for(g in 1:length(G)){
      Gnames<-rownames(TB)[ellipse.groups==G[g]]
      Gcenter<-c(mean(sam[Gnames,1]),mean(sam[Gnames,2]))
      Gsigma<-var(sam[Gnames,1:2])
      ellipse(Gcenter, Gsigma, add = TRUE,center.pch=ell.center.pch,
              radius=ifelse(length(ell.radius)>1,ell.radius[g],ell.radius),
              col=ifelse(length(ell.cols)>1,ell.cols[g],ell.cols),
              fill=TRUE,fill.alpha=0.05) # car package
    }
  }
  
  if(sites.labs==TRUE){
    shadowtext(sam[,1],sam[,2],labels=rownames(sam),xpd=TRUE, ... )  
  }
  
  if(points==TRUE){
    points(sam[,1],sam[,2],xpd=TRUE, pch=pch,bg=bg, ... )
  }
  
  
  
  
}

