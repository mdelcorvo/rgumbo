#' LiftoverVcf
#'
#' Function to lifts over a VCF file from one genome build to another, producing a properly headered, sorted VCF in one go.
#' 
#' \code{LiftoverVcf(vcf,output,from="hg38",to="hg19")}.
#'
#' @param vcf The input VCF/BCF file to be lifted over
#' @param output name of the output file
#' @param from starting assembly
#' @param to final assembly
#' @return a messagge explaining if all packages are installed.
#'

#' @examples
#' LiftoverVcf(vcf,output,from="hg38",to="hg19")
#' @name LiftoverVcf

#' @export
#' @rdname LiftoverVcf

#' @export
#' @rdname LiftoverVcf

LiftoverVcf <- function(vcf,output,from="hg38",to="hg19") {

defined <- ls()
passed <- names(as.list(match.call())[-1])

if (any(!defined %in% passed)) {
  stop(paste("missing values for", paste(setdiff(defined, passed), collapse=", ")))
}
    
chain <- get(paste(from,'To',to,sep='')) # load LiftOver annotation

header_pos<- head(grep('#',readLines(vcf)),-1) #Line numbers of header ('#')
vec1<-fread(vcf,fill=TRUE,header=F)[header_pos,]$V1

vcf<- fread(vcf,skip=grep('#CHROM',readLines(vcf))-1)
colnames(vcf)[1]<-c('Chrom')
nlines<-nrow(vcf)

header<-as.data.frame(matrix(nrow=length(vec1),ncol=ncol(vcf)))
header$V1<-vec1

# check presence of 'chr' prefix in chromosome number
if (!all(grepl('chr',vcf$Chrom))) {
  vcf$Chrom <- paste('chr',vcf$Chrom,sep='')
}

# LiftOver
g1<-GRanges(vcf$Chrom, ranges =IRanges(vcf$POS,vcf$POS))
lift<-do.call(rbind.data.frame, lapply(lapply(liftOver(g1, chain),as.character),function(x) ifelse(length(x)==0,NA,x)))
colnames(lift)[1]<-to
list1<-strsplit(lift[,1],':',fixed = TRUE)
vcf$Chrom<-unlist(lapply(list1,function(x) x[1]))
vcf$POS<- unlist(lapply(list1,function(x) x[2]))
#

colnames(vcf)[1]<-c('#CHROM')
colnames(header)<-colnames(vcf)
vcf<-rbind(colnames(vcf),as.data.frame(vcf))
vcf<-rbind(header,vcf) # restore header

fwrite(vcf,file=output, row.names=F, col.names=F,quote=F, sep='\t')

cat(nlines,'positions successfully converted from',from,'to',to,'genome build in',output,'file',  "\n")
}
