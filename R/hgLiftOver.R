#' hgLiftOver
#'
#' Function to converts genome coordinates between assemblies. 
#' 
#' \code{hgLiftOver(bed,from="hg38",to="hg19",chrom=NULL,start=NULL,end=NULL)}.
#'
#' @param bed data.frame in BED format with chromosome number, start and end position or chromosome number and a single genomic position respectively
#' @param from starting assembly
#' @param to final assembly
#' @param chrom Chromosome number
#' @param start Start position
#' @param end End position
#' @return a messagge explaining if all packages are installed.
#'

#' @examples
#' hgLiftOver(bed,from="hg38",to="hg19")
#' @name hgLiftOver

#' @export
#' @rdname hgLiftOver

#' @export
#' @rdname hgLiftOver

hgLiftOver <- function(bed,from="hg38",to="hg19",chrom=NULL,start=NULL,end=NULL) {
    
    chain <- get(paste(from,'To',to,sep='')) # load LiftOver annotation
    if (ncol(bed)>2) {	
			if (is.null(chrom)) {
			colnames(bed)[1:3]<-c('Chrom','Start','End')
			bed_old<-bed
				
				if (!all(grepl('chr',bed$Chr))) {
				bed$Chrom <- paste('chr',bed$Chrom,sep='')
				}
				
			g1<-GRanges(bed$Chr, ranges =IRanges(bed$Start,bed$End))
			lift<-do.call(rbind.data.frame, lapply(lapply(liftOver(g1, chain),as.character),function(x) ifelse(length(x)==0,NA,x)))
			colnames(lift)[1]<-to
							   
			list1<-strsplit(lift[,1],':',fixed = TRUE)
			list2<-strsplit(unlist(lapply(list1,function(x) x[2])),'-',fixed=T)
						  
			bed$Chrom<-unlist(lapply(list1,function(x) x[1]))
			bed$Start<- unlist(lapply(list2,function(x) x[1]))
			bed$End<-unlist(lapply(list2,function(x) x[2]))
			bed$refGenome <- to	
					       
				if (!all(grepl('chr',bed_old$Chr))) {
				bed$Chrom <- gsub('chr','',bed$Chrom)
				}	
					       
			} else if (is.null(end)) {
			bed<-as.data.frame(bed)	
			bed_old<-bed
	    
				if (!all(grepl('chr',bed[,chrom]))) {
				bed[,chrom] <- paste('chr',bed[,chrom],sep='')
				}
	    
			g1<-GRanges(bed[,chrom], ranges =IRanges(bed[,start],bed[,start]))
			lift<-do.call(rbind.data.frame, lapply(lapply(liftOver(g1, chain),as.character),function(x) ifelse(length(x)==0,NA,x)))
			colnames(lift)[1]<-to

			list1<-strsplit(lift[,1],':',fixed = TRUE)

			bed[,chrom]<-unlist(lapply(list1,function(x) x[1]))
			bed[,start]<- unlist(lapply(list1,function(x) x[2]))
			bed$refGenome <- to
				   
				if (!all(grepl('chr',bed_old[,chrom]))) {
				bed[,chrom] <- gsub('chr','',bed[,chrom])
				}
			
			} else {
			bed<-as.data.frame(bed)
			bed_old<-bed
			
				if (!all(grepl('chr',bed[,chrom]))) {
				bed[,chrom] <- paste('chr',bed[,chrom],sep='')
				}
				
			g1<-GRanges(bed[,chrom], ranges =IRanges(bed[,start],bed[,end]))
			lift<-do.call(rbind.data.frame, lapply(lapply(liftOver(g1, chain),as.character),function(x) ifelse(length(x)==0,NA,x)))
			colnames(lift)[1]<-to
							   
			list1<-strsplit(lift[,1],':',fixed = TRUE)
			list2<-strsplit(unlist(lapply(list1,function(x) x[2])),'-',fixed=T)
						  
			bed[,chrom]<-unlist(lapply(list1,function(x) x[1]))
			bed[,start]<- unlist(lapply(list2,function(x) x[1]))
			bed[,end]<-unlist(lapply(list2,function(x) x[2]))
			bed$refGenome <- to
			
				if (!all(grepl('chr',bed_old[,chrom]))) {
				bed[,chrom] <- gsub('chr','',bed[,chrom])
				}
			}
    } else {

		        colnames(bed)[1:2]<-c('Chrom','Pos')
		        bed_old<-bed
	    
				if (!all(grepl('chr',bed$Chr))) {
				bed$Chrom <- paste('chr',bed$Chrom,sep='')
				}
	    
			g1<-GRanges(bed$Chr, ranges =IRanges(bed$Pos,bed$Pos))
			lift<-do.call(rbind.data.frame, lapply(lapply(liftOver(g1, chain),as.character),function(x) ifelse(length(x)==0,NA,x)))
			colnames(lift)[1]<-to

			list1<-strsplit(lift[,1],':',fixed = TRUE)

			bed$Chrom<-unlist(lapply(list1,function(x) x[1]))
			bed$Pos<- unlist(lapply(list1,function(x) x[2]))
			bed$refGenome <- to
				   
				if (!all(grepl('chr',bed_old$Chr))) {
				bed$Chrom <- gsub('chr','',bed$Chrom)
				}
		
	}		    
    return(bed)
}

