#' SnpEff
#'
#' Function to make the SnpEff-based vcf readable by extracting the functional annotation and splitting each gene transcript into single rows.
#' 
#' \code{SnpEff(vcf,modifier=F,amino=F,rm_dup=T)}.
#'
#' @param vcf SnpEff-based annotated vcf file
#' @param modifier if TRUE variants with no impact on gene will be removed
#' @param amino if TRUE variants that do not change the amino acid sequence will be removed
#' @param rm_dup if TRUE duplicated transcripts for each variant will be removed
#'

#' @examples
#' SnpEff(vcf)
#' @name SnpEff

#' @export
#' @rdname SnpEff

#' @export
#' @rdname SnpEff

SnpEff <- function(vcf,modifier=F,amino=F,rm_dup=T) {
#vcf file with SnpEff annotation
raw<- as.data.frame(fread(vcf,skip=grep('#CHROM',readLines(vcf))-1))
raw$id<-paste(raw[,1],raw[,2],sep=':')

tmp1<-strsplit(raw$INFO,'|',fixed = TRUE)
names(tmp1)<-raw$id
l1<-lapply(tmp1,get_code) # get relevant information from SnpEff annotation
d1<-do.call(rbind.data.frame, l1)
d1$Genomic_coordinates <- gsub('[.].*','',row.names(d1))
d1<-d1[,c(8,1:7)]
if (!modifier)
	d1<-d1[d1$impact!='MODIFIER',] # remove variants with no impact on gene
if (amino)
	d1<-d1[d1$amino!='',]	# remove variants that do not change the amino acid sequence 
d1<-merge(d1,raw[,c(3:7,9:ncol(raw))],by.x='Genomic_coordinates',by.y='id')
if (rm_dup) {
	# remove duplicated transcripts for each variant
	list1<-split(d1,d1$Genomic_coordinates)
	list2<-lapply(list1,function(x) x[!duplicated(x$dna),])
	d1<-do.call(rbind.data.frame, list2)
	d1$transcript_number<-unlist(lapply(list1,
	  function(x) apply(as.array(unique(x$dna)),1,
	   function(y) paste(x[x$dna==y,]$transcript_number,collapse='|'))))
}
row.names(d1)<-NULL   
return(d1)	
}

get_code<-function(x,amino=F) {
toMatch<-c("MODIFIER","LOW", "MODERATE", "HIGH") #effect impact on gene
v1<-grep(paste(toMatch,collapse="|"),x)	
impact=x[v1][x[v1]!=""]
if (!is.null(impact)) {
	gene=x[v1+1]
	gene=tail(names(sort(table(gene))),1)
	ensembl_id=x[v1+2]
	effect=x[v1-1]
	dna= x[v1+7]
	if (is.null(dna))
		dna=''
	amino= x[v1+8]
	if (is.null(amino))
		amino=''
	transcript_number<-c(1:length(impact))
	data1<- data.frame(gene,ensembl_id,transcript_number,impact,effect,dna,amino)
	return(data1)
}
}	
