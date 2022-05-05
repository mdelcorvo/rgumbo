#' hgLiftOver
#'
#' Function to converts genome coordinates between assemblies. 
#' 
#' \code{hgLiftOver(df,from="hg38",to="hg19")}.
#'
#' @param df data.frame in BED format with chromosome number, start and end position or chromosome number and a single genomic position respectively
#' @param from starting assembly
#' @param to final assembly
#' @return a messagge explaining if all packages are installed.
#'

#' @examples
#' hgLiftOver(df = bed_file,from="hg38",to="hg19")
#' @name hgLiftOver

#' @export
#' @rdname hgLiftOver

#' @export
#' @rdname hgLiftOver

hgLiftOver <- function(df,from="hg38",to="hg19") {
    
    chain <- get(paste(from,'To',to,sep=''))
    if (ncol(df)>2) {	
	    colnames(df)[1:3]<-c('Chrom','Start','End')
	    df_old<-df

	    if (!all(grepl('chr',df$Chr))) {
		df$Chrom <- paste('chr',df$Chrom,sep='')
	    }

	    g1<-GRanges(df$Chr, ranges =IRanges(df$Start,df$End))
	    lift<-do.call(rbind.data.frame, lapply(lapply(liftOver(g1, chain),as.character),function(x) ifelse(length(x)==0,NA,x)))
	    colnames(lift)[1]<-to
						   
	    list1<-strsplit(lift[,1],':',fixed = TRUE)
	    list2<-strsplit(unlist(lapply(list1,function(x) x[2])),'-',fixed=T)
					  
	    df$Chrom<-unlist(lapply(list1,function(x) x[1]))
	    df$Start<- unlist(lapply(list2,function(x) x[1]))
	    df$End<-unlist(lapply(list2,function(x) x[2]))
	    df$refGenome <- to
				  
	    if (!all(grepl('chr',df_old$Chr))) {
		df$Chrom <- gsub('chr','',df$Chrom)
	    }
    } else {
	    colnames(df)[1:2]<-c('Chrom','Pos')
	    df_old<-df
	    
	    if (!all(grepl('chr',df$Chr))) {
		df$Chrom <- paste('chr',df$Chrom,sep='')
	    }
	    
	    g1<-GRanges(df$Chr, ranges =IRanges(df$Pos,df$Pos))
	    lift<-do.call(rbind.data.frame, lapply(lapply(liftOver(g1, chain),as.character),function(x) ifelse(length(x)==0,NA,x)))
	    colnames(lift)[1]<-to
		
	    list1<-strsplit(lift[,1],':',fixed = TRUE)
						   
	    df$Chrom<-unlist(lapply(list1,function(x) x[1]))
	    df$Pos<- unlist(lapply(list1,function(x) x[2]))
	    df$refGenome <- to
				   
	    if (!all(grepl('chr',df_old$Chr))) {
		df$Chrom <- gsub('chr','',df$Chrom)
	    }
    }	    
    return(df)
}
