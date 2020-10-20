add_Abundance_Polygons<-function(M,multivA,xlim=par("usr")[1:2],ylim=par("usr")[3:4],PerRange=c(0,0.1,1,2,5,10,15,20),PerCex=seq(0.02,0.1,length.out=length(PerRange)),
                                 xcorr=rep(0,length.out=ncol(M)),ycorr=rep(0,length.out=ncol(M)),labels=NULL,forced_labels=NA,MinLabPer=1.5,expansion=NULL,
                                 col=rgb(t(col2rgb("grey60")),alpha=100,maxColorValue=255),col.labs="grey55", ... ){
  # Analysing the individuals/cases/species data
  if(grepl("prcomp",paste(class(multivA),collapse=" "))==TRUE){
    SpeciesVars<-multivA$rotation[,1:2]      
  }else{
    if(grepl("metaMDS",paste(class(multivA),collapse=" "))==TRUE){
      SpeciesVars<-multivA$species[,1:2]        
    }else{
      stop("multivA object is not of class prcomp or meta|monoMDS")
    }
  }
  
  # Adding new plotting space where abundances are drawn based on a expansion factor over previous plot axis.
  if(!is.null(expansion)){ 
    xlim<-xlim*expansion
    ylim<-ylim*expansion
    par(new=TRUE)
    plot(NA,NA,xlim=xlim,ylim=ylim,axes=FALSE,xlab="",ylab="", ... )
    scale<-1/expansion
    ylimfactor<-0.01*(abs(par("usr")[3])+abs(par("usr")[4]))
    text(par("usr")[2],par("usr")[4]+ylimfactor,labels=paste("Expansion scale: 1/",scale,sep=""),font=2,col="grey60",xpd=TRUE,cex=1.2,adj=c(1,0))
  }
  
  # Getting Global relative abundance
  TaxonSum<-apply(t(M),1,sum)
  TaxPer<-(TaxonSum*100)/sum(TaxonSum)
  TaxCex<-vector("numeric")
  for(s in 1:length(TaxPer)){
    TaxCex[s]<-PerCex[sum(PerRange<=TaxPer[s])]
  }
  names(TaxCex)<-names(TaxPer)
  # print(TaxPer)
  
  # Getting constant factor to modify squares sizes according to 
  # axis limits if x.axis and y.axis don't have same length
  standarLimitLength<-abs(-1)+abs(1) # To make the range of squares even respecting to a plot x(-1,1)y(-1,1)
  # Basically is equal to 2
  xfactor<-(abs(xlim[1])+abs(xlim[2]))/standarLimitLength
  yfactor<-(abs(ylim[1])+abs(ylim[2]))/standarLimitLength
  # Calculating squares dimensions
  polyx<-cbind(SpeciesVars[,1]-(TaxCex*xfactor),SpeciesVars[,1]+(TaxCex*xfactor),SpeciesVars[,1]+(TaxCex*xfactor),SpeciesVars[,1]-(TaxCex*xfactor))
  polyy<-cbind(SpeciesVars[,2]-(TaxCex*yfactor),SpeciesVars[,2]-(TaxCex*yfactor),SpeciesVars[,2]+(TaxCex*yfactor),SpeciesVars[,2]+(TaxCex*yfactor))
  label.x.adj<-TaxCex*xfactor # To adjust labels over x.axis
  
  # Plotting Global relative abundance squares
  for(s in 1:nrow(polyx)){
    polygon(polyx[s,],polyy[s,],border=ifelse(length(col)>1,col[s],col),col=ifelse(length(col)>1,col[s],col),xpd=TRUE,lwd=1) 
  }
  if(!is.null(labels)){
    for(s in 1:nrow(SpeciesVars)){
      if(TaxPer[s]>MinLabPer && grepl("Unassigned",rownames(SpeciesVars)[s])==FALSE || grepl(paste(forced_labels,collapse="|"),rownames(SpeciesVars)[s])==TRUE){
        if(SpeciesVars[s,1]>0){
          adj<-c(0,0.5)
        }else{
          adj<-c(1,0.5)
          label.x.adj[s]<-label.x.adj[s]*(-1)
        }
        shadowtext(SpeciesVars[s,1]+label.x.adj[s]+xcorr[s],SpeciesVars[s,2]+ycorr[s],labels=labels[s],cex=1,font=2,xpd=TRUE,adj=adj,bg="white",
                   col=ifelse(length(col.labs)>1,col.labs[s],col.labs),...)
      }
    }  
  }
  return(c(xfactor,yfactor))
  
}  
