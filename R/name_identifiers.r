
#' @name name_identifiers
#' @title Name identifiers for EEG channels
#' @description 
#' Generates named identifiers for EEG channels based on unique pairs of channels.
#' 
#' @title Name identifiers for EEG channels
#' 
#' @description 
#' Generates named identifiers for EEG channels based on unique pairs of channels.
#' 
#' @param n Numeric vector, representing EEG channel numbers
#' 
#' @return Character vector, representing named identifiers for EEG channels
#' 
#' @examples 
#' name_identifiers(c(1:4)) # Returns: "1" "2" "3" "4" "1.2" "1.3" "1.4" "2.3" "2.4" "3.4"
name_identifiers<-function (n) 
{
    n <- data.frame(1:length(n), n) # Create a data frame 
    D <- nrow(n)
    n2 <- idlist <- define_identifier(D)
    n2 <- cbind(n2, "/")
    for (id in (idlist)) {
        if (id <= D) {
            n2[n2[, 1] == id, 2] <- as.character(n[n[, 1] == 
                id, 2])
        }
        else {
            if (id < 1000) {
                a <- as.numeric(substr(id, 1, 1))
                b <- as.numeric(substr(id, 2, 3))
            }
            else {
                a <- as.numeric(substr(id, 1, 2))
                b <- as.numeric(substr(id, 3, 4))
            }
            n2[n2[, 1] == id, 2] <- paste(as.character(n[n[, 
                1] == a, 2]), as.character(n[n[, 1] == b, 2]), 
                sep = ".")
        }
    }
    return(n2[, 2])
}