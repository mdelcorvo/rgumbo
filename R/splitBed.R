#' splitBed
#'
#' Function to split bed files by chunks, by number of rows for each chromosome or by chromosome
#'
#' \code{splitBed(bed,chunks,n,chrOnly=F,prefix=T,writeBed=F,verbose=T)}.
#'
#' @param bed data.frame in BED format with chromosome number, start and end position
#' @param chunks number of chunks to split
#' @param n number of lines for each chromosome
#' @param chrOnly split only by chromosome
#' @param prefix set to TRUE to add chr prefix if absent or FALSE to remove it if present
#' @param writeBed write new bed files
#' @param verbose give verbose output
#'

#' @examples
#' splitBed(bed,n=100)
#' @name splitBed

#' @export
#' @rdname splitBed

#' @export
#' @rdname splitBed

splitBed<-function(bed,chunks=NULL,n=NULL,chrOnly=F,prefix=T,writeBed=F,verbose=F) {

df<-data.table::fread(bed)
	if (prefix) {
		if (!all(grepl('chr',df$V1)))
			df$V1 <- paste('chr',df$V1,sep='')
	} else {
		if (all(grepl('chr',df$V1)))
			df$V1 <- gsub('chr','',df$V1)
	}	
colnames(df)[1:3]<-c('V1','V2','V3')

if (!is.null(chunks) & is.numeric(chunks)) {  # Split bed into equal-sized chunks

	l3 <- split(df, factor(sort(rank(row.names(df))%%chunks)))
	names(l3)<-unlist(lapply(l3,function(x) paste(paste(x[1,1],paste(x[1,2],x[1,3],sep='-'),sep=':'),
												   paste(tail(x,1)[,1],paste(tail(x,1)[,2],tail(x,1)[,3],sep='-'),sep=':'),sep='_')))

} else if (!is.null(n) & is.numeric(n)) { # Split bed by number of rows for each chromosome

	l1<-split(df,df$V1)
	l2<- lapply(l1,function(x) split(x, (seq(nrow(x))-1) %/% n))
	l3<-unlist(l2,recursive=F)
	names(l3)<-unlist(lapply(l3,function(x) paste(unique(x$V1),paste(x$V2[1],tail(x$V3,1),sep='-'),sep=':')))

} else if (chrOnly) { # Split bed by chromosome

	l3<-split(df,df$V1)

} else {

	stop("Muste provide at least one parameter for bed splitting")

}

if (verbose) {
		res<-data.frame(N.segments=length(l3),Min.length=as.numeric(summary(unlist(lapply(l3,nrow)))[1]),Max.length=as.numeric(summary(unlist(lapply(l3,nrow)))[6]))
		print.data.frame (res)
}

if (writeBed) {
    lapply(l3,function(x) fwrite(x,file=paste0(paste(unique(x$V1),paste(x$V2[1],tail(x$V3,1),sep='-'),sep='_'),'.bed'), row.names=F, col.names=F,quote=F, sep='\t'))
	cat(bed, ' file splitted in ',length(l3),' pieces', "\n")
	} else {
	return(l3)
}

}
