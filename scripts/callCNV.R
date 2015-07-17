
callCNV<-function(gcNorm,samples,ref="B73",limit=1.2,limitHom=6){
  
  # make a matrix of 0, representing the "no evidence of change in copy-number"
  cnv.mat<-matrix(data=2,nrow=dim(gcNorm)[1],ncol=length(samples))
  # give column names
  colnames(cnv.mat)<-paste0(samples,"_CNV")
  #convert to data.table
  cnv.mat<-data.table<-(cnv.mat)
  # combine gcNorm and cnv.mat
  gcNorm<-cbind(gcNorm,cnv.mat)
  lims<-NULL
  stddevs<-NULL
  colnames(cnv.mat)<-NULL
  # make a matrix to store a q-value style number for each gene
  q.mat<-data.table(cnv.mat)
  # cycle through the samples
  for ( s in samples ) {
    # subset the data
    sub.df<-subset(gcNorm,select=c(s))
    sub.df<-data.frame(sub.df)
    #p1<-hist(sub.df[,s], main =paste0(s, " read-depth"),
    #   xlab="normalized coverage per kb",cex.lab=1.4, breaks = c(seq(0,20,0.2),20.00001,max(sub.df[,s])), xlim=c(0,21), 
    #   col = alpha("cornflowerblue",0.5),border=alpha("cornflowerblue",0.7))
    
    # calculate the ratio between sample and reference"
    ratio<-sub.df[,s]/subset(gcNorm,select=ref)
    ratio<-data.matrix(ratio)
    ratio<-log(ratio+0.001)
    stdev<-sd(ratio)
    d<-(mean(ratio)-ratio)
    # plot the ratio as a histogram
    #p2<-hist(ratio, main=paste0(s," ratio (sample/ref)"),
    #   cex.lab=1.4, xlab="log(coverage difference (sample/ref))", col=alpha("cornflowerblue",0.4),
    #   border=alpha("cornflowerblue",0.4), breaks=200)
    # calculate the CNV "limit
    linePos<-mean(ratio)+(limit*stdev)
    linePosHom<-mean(ratio)+(limitHom*stdev)
    lims<-c(lims,min(mean(ratio)+linePos),abs(mean(ratio)-linePos))
    #abline(v=mean(ratio), col = "black", lwd=2, lty = "solid",lend=1)
    
    #turn sub.df back into a data.table
    sub.df<-data.table(sub.df,ratio)
    ratio.df<-data.table(sub.df,ratio)
    sName<-paste0(s,"_CNV")
    # call down CNV
    calls<-(gcNorm[,sName,with=FALSE])
    calls<-as.vector(as.matrix(calls))
    calls[ratio<(mean(ratio)-linePos)]<-1
    calls[ratio<(mean(ratio)-linePosHom)]<-0
    
    # call up CNV
    calls[ratio>(mean(ratio)+linePos)]<-3
    
    # attribute 
    gcNorm[,sName:=calls,with=FALSE]
  
    ##double check that the calls are in the right place
    #abline(v=ratio[calls==1], col=alpha("red",0.1), lwd=0.1, lty="solid", lend=1)
    #abline(v=ratio[calls==-1], col=alpha("blue",0.01), lwd=0.1, lty="solid", lend=1)
    #abline(v=ratio[calls==-2], col=alpha("blue",0.1), lwd=0.1, lty="solid", lend=1)
    #abline(v=mean(ratio)-abs(linePos), col = "black", lwd=1, lty = "dashed",lend=1)
    #abline(v=mean(ratio)+(linePos), col = "black", lwd=1, lty = "dashed",lend=1)
    #abline(v=mean(ratio)-abs(linePosHom), col = "black", lwd=1, lty = "dashed",lend=1)
    #abline(v=mean(ratio)+(linePosHom), col = "black", lwd=1, lty = "dashed",lend=1)
    #now the qvalue
    qVal<-pnorm(q=d,mean=0,sd=stdev)
    q.mat[,sName:=-log(qVal),with=FALSE]
    #print(p1)
    #print(p2)
  }
  #dev.off()
  output<-list(cnvs=gcNorm,q.mat=q.mat)
}# function callCNV()