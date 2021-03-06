#Actually get 200bp file

require('BSgenome.Hsapiens.UCSC.hg38')
require('rtracklayer')

#process input arguments
args <- commandArgs(trailingOnly = TRUE)


process.simple.offsets.python<-function(data,prefix="./",window){
  colnames(data)[1:6]<-c("A","B","C","D","E","F")
  data<-data[!is.na(data$D),]
  data$D=as.numeric(as.character(as.matrix(data$D)))
  
  
  halfw=window/2
  data$C=gsub("chr","",data$C)
  data$C=as.character(data$C)
  data$C[as.character(data$C)=="23"]<-"X"
  data$C[as.character(data$C)=="24"]<-"Y"
  
  #filter by chromosome length
  chrends<-seqlengths(BSgenome.Hsapiens.UCSC.hg38)[match( paste("chr",data$C,sep=""),names(seqlengths(BSgenome.Hsapiens.UCSC.hg38)))]
  data=data[ (as.numeric(as.character(data$D))-(halfw-1))>0 &  (as.numeric(as.character(data$D))+halfw)<chrends,]

  data.rd<-GRanges(ranges=IRanges(start=as.numeric(as.character(data$D))-(halfw-1),as.numeric(as.character(data$D))+halfw),seqnames=paste("chr",as.character(data$C),sep=""),ori=data$A,mut=data$B,name=data$F)
  genome(data.rd) <- 'hg38'
  #retrieve reference genome sequences
  data.wt<-as.character(as.matrix(as.data.frame(getSeq(x=Hsapiens,names=as(data.rd,"GRanges")))))
  
  #write to fasta format
  temp<-DNAStringSet(data.wt)
  names(temp)<-paste(data.rd$ori,data.rd$mut,seqnames(data.rd),start(data.rd),end(data.rd),data.rd$name,sep = "_")
  writeXStringSet(temp,filepath=paste(prefix,".wt",window-100,".fasta",sep=""))
  print(data.rd)
}

#read vcf format

if(grepl("vcf$",args[1])){
  data=read.csv(args[1],sep='\t',header=F, comment.char = "#",colClasses = c("character"))
  data[,2]=as.numeric(data[,2])
  data=data[,c(4,5,1,2,2,3)]
  colnames(data)[1:6]<-c("A","B","C","D","E","F")
  #for vcf format, we retrieve 300bp sequences (instead of 200bp) to allow for small deletions (<100bp)
  #process.simple.offsets.python(data,window = 300,prefix=args[1])
  process.simple.offsets.python(data,window = as.integer(args[2]),prefix=args[1])
}
