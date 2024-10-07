#' @name define_identifier
#' @title Define identifiers for EEG channels
#' @description 
#' Define identifiers for EEG channels to represent unique pairs of channels. 
#' The identifiers are generated based on the number of EEG channels provided.
#' 
#' @param D Numeric, number of the EEG channels 
#' @return Numeric vector, representing identifiers for EEG channels
#' 
#' @examples 
#' define_identifier(4)  # Returns: 1 2 3 4 102 103 104 203 204 304
define_identifier<-function (D) 
{   
    if (D == 1) {return(1)}
    else if (D == 2) {return(c(1,2,102))}
    else{
    tmp <- cbind(1, (1:D))
    for (i in 2:D) tmp <- rbind(tmp, cbind(i, (i:D)))
    tmp <- tmp[-which(tmp[, 1] == tmp[, 2]), ]
    ids <- c(sprintf("%02d", 1:D), paste(sprintf("%02d", tmp[, 
        1]), sprintf("%02d", tmp[, 2]), sep = ""))
    ids <- as.numeric(ids)
    return(ids)}
}

