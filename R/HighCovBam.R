#' HighCovBam
#'
#' Function to find high coverage regions in a bam file.
#' 
#' \code{HighCovBam(bam,cutoff=NULL,width=500)}.
#'
#' @param bam bam file to analyze
#' @param cutoff threshold for minimum coverage
#' @param width size of the windows
#'

#' @examples
#' HighCovBam(bam,100)
#' @name hgLiftOver

#' @export
#' @rdname HighCovBam

#' @export
#' @rdname HighCovBam

HighCovBam <- function(bam,cutoff=NULL,width=500) {
	rbga<-readGAlignments(bam)
	bamCov<-coverage(rbga)

	if (is.null(cutoff))
       cutoff<-mean(quantile(bamCov,.9))

	hotslices<-slice(bamCov,lower=cutoff)
	hotwindows<-resize(hotslices,width,fix="center")

	data1<-data.frame()
	for (i in c(1:23)) {
		tmp<-as.data.frame(hotwindows@listData[i]$chr@ranges)
        if (nrow(tmp)>0) {
           tmp$chr<-ifelse(i<=22,paste0('chr',i),paste0('chr','X'))
           data1<-rbind(data1,tmp)
        }
	}
	res<-data1[,c(4,1:3)]
	return(res)
}
