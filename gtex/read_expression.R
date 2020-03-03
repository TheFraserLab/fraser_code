#' Efficently reads the expression or read count data from gtex
#'
#' @param sample character — vector of sample names in the long GTEx format (e.g. GTEX-WXYG-2026-SM-4E3IY)
#' @param expressionFile character — path to expression or read count uncompressed file. These can be downloaded from the GTEx portal (GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz)
#'
#' @return numeric matrix with rownames being genes and colnames being samples
#'
#' @examples
#' DO_NOT_RUN{
#' exp_mat <- readExpression('GTEX-WXYG-2026-SM-4E3IY', './GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct')
#'
#'}
readExpression <- function(samples = NULL, expressionFile) {    
                                                                                    
    sep <-  "\t"                                         
    filecon <- file(expressionFile, "r")    
    commentLines <- readLines(filecon, n = 2)    
                   
    header <- unlist(strsplit(readLines(filecon, n = 1), sep))    
                                                     
    if(is.null(samples)) {                                                                                                                                                                             
        sampleIndex <- rep(T, length(header))                                                                                
        sampleIndex[2] <- F                                                                                                                                                           
    } else {                                                                 
        sampleIndex <- header %in% samples                                                                                                                                                    
    }                                                                                            
                                                                                                               
    sampleIndex[1] <- T      
                       
    colClasses <- rep("NULL", length(header))                                                                                                                               
    colClasses[sampleIndex] <- "numeric"                                        
    colClasses[1] <- "character"                                
                                  
    # Stop if not samples found                         
    if(length(colClasses) == 0)    
        stop('No samples found in RNA-seq matrix')    
         
    result <- read.table(filecon, header = F, sep = sep, colClasses = colClasses)     
        
    colnames(result) <- header[sampleIndex]    
        
    close(filecon)    
        
    result[,1] <- gsub('\\..+', '', result[,1])    
        
    rowNames <- result[,1]    
    result <- result[,-1,drop=F]    
        
    colNames <- colnames(result)    
    result <- as.matrix(result)    
        
    if(!is.numeric(result))    
        stop('Something went wrong, non-numeric values found in gene expression matrix')    
        
    rownames(result) <- rowNames    
    colnames(result) <- colNames    
    return(as.matrix(result))    
        
}
