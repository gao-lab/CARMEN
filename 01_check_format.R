#process input arguments
require('rtracklayer')
# Do not use scientific notation in the postion of variants
options(scipen = 200)

args <- commandArgs(trailingOnly = TRUE)

#read vcf format

if(grepl("vcf$",args[1])){
  data=read.csv(args[1],sep='\t',header=F, comment.char = "#",colClasses = c("character"))
  data[,2]=as.numeric(data[,2])
  colnames(data)[1:6]<-c("A","B","C","D","E","F")
  data=data[order(data$A,data$B),]
  file_path = args[2]
  
  if(args[3] == "GRCh37"){
  #===========Convert the hg19 to hg38==================
   #"To avoid the differece between hg38Tohg19 and the hg19, we conver hg19 to hg38 and vice again."
  chain.file <- "hg19ToHg38.over.chain"
  ch <- import.chain(chain.file)
  input <- GRanges( seqnames = Rle(data$A),
                    ranges = IRanges(data$B, end=data$B),
                    strand = Rle(strand("*")) )
  out_hg38 <- as.data.frame(liftOver(input, ch))

  out_hg38_vcf <- data[out_hg38$group,]
  out_hg38_vcf$A <- as.character(out_hg38$seqnames)
  out_hg38_vcf$B <- out_hg38$start

  data <- out_hg38_vcf
  chain.file <- "hg38ToHg19.over.chain"
  ch <- import.chain(chain.file)
  input <- GRanges( seqnames = Rle(data$A),
                    ranges = IRanges(data$B, end=data$B),
                    strand = Rle(strand("*")) )
  out_hg19 <- as.data.frame(liftOver(input, ch))

  out_hg38_vcf <- data[out_hg19$group,]
  out_hg19_vcf <- out_hg38_vcf
  out_hg19_vcf$A <- as.character(out_hg19$seqnames)
  out_hg19_vcf$B <- out_hg19$start
  
  out_hg38_vcf <- out_hg38_vcf[ order(as.numeric(row.names(out_hg38_vcf))), ]
  out_hg19_vcf <- out_hg19_vcf[ order(as.numeric(row.names(out_hg19_vcf))), ]

  out_hg19_bed <- cbind(out_hg19_vcf$A,out_hg19_vcf$B-1,out_hg19_vcf$B)

  write.table(out_hg19_vcf,file=paste(file_path,"/infile.hg19.vcf",sep=""),sep = '\t',col.names = F,row.names = F,quote=F)
  write.table(out_hg38_vcf,file=paste(file_path,"/infile.hg38.vcf",sep=""),sep = '\t',col.names = F,row.names = F,quote=F)
  write.table(out_hg19_bed,file=paste(file_path,"/infile.hg19.bed",sep=""),sep = '\t',col.names = F,row.names = F,quote=F)}
  
  if(args[3] == "GRCh38"){
  #===========Convert the hg38 to hg19==================
  chain.file <- "hg38ToHg19.over.chain"
  ch <- import.chain(chain.file)
  input <- GRanges( seqnames = Rle(data$A),
                    ranges = IRanges(data$B, end=data$B),
                    strand = Rle(strand("*")) )
  out_hg19 <- as.data.frame(liftOver(input, ch))

  out_hg38_vcf <- data[out_hg19$group,]
  out_hg19_vcf <- out_hg38_vcf
  out_hg19_vcf$A <- as.character(out_hg19$seqnames)
  out_hg19_vcf$B <- out_hg19$start

  out_hg38_vcf <- out_hg38_vcf[ order(as.numeric(row.names(out_hg38_vcf))), ]
  out_hg19_vcf <- out_hg19_vcf[ order(as.numeric(row.names(out_hg19_vcf))), ]
  out_hg19_bed <- cbind(out_hg19_vcf$A,out_hg19_vcf$B-1,out_hg19_vcf$B)

  write.table(out_hg19_vcf,file=paste(file_path,"/infile.hg19.vcf",sep=""),sep = '\t',col.names = F,row.names = F,quote=F)
  write.table(out_hg38_vcf,file=paste(file_path,"/infile.hg38.vcf",sep=""),sep = '\t',col.names = F,row.names = F,quote=F)
  write.table(out_hg19_bed,file=paste(file_path,"/infile.hg19.bed",sep=""),sep = '\t',col.names = F,row.names = F,quote=F)}}

 
