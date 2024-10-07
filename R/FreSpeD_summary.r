#' @name FreSpeD_summary
#' @title Summary of FreSpeD output
#' @description Various tables and plots showing change point distribution in time, frequency and for individual
#' channels and channel pairs.
#' 
#' @param cp list of numeric arrays, List of change points per channel/channel pair given by FreSpeD function
#' @param channelNames list, default empty.
#' @param windoLen numeric, default 1
#' @param plot logical, Should global summary plots be provided?
#' 
#' @return list containing
#'         1. numeric, total number of change points 
#'         2. numeric, total number of change points in autospectra 
#'         3. numeric vector, change points over time 
#'         4. numeric vector, autospectral change points over time 
#'         5. numeric vector, change points over frequency bands 
#'         6. numeric vector, autospectral change points over frequency bands 
#'         7. dataframe, containing the number and location for each channel and channel pair
#'
FreSpeD_summary<-function (cp, channelNames = c(), windowLen = 1, plot = TRUE) 
{
    if (length(cp) == 0) 
        stop("The provided change point list of spectral quantities is empty")
    if (windowLen < 1 || windowLen%%1 != 0) 
        stop("Window length has to be a positive integer")
    # The length if cp is by definition c = (D + choose(2,D) ), where D is the number of channels. So given c, we get D by the following formula     
    D <- round_right(uniroot(function(x, c = length(cp)) { 
        x^2 + x - 2 * c
    }, c(1, 1000))$root)

    if (length(channelNames) > 0 && length(channelNames) != D) 
        stop("Bad input. The length of channelNames is not equal to the number of channels given by the list of spectral quantities")
    if (length(channelNames) == 0) 
        channelNames <- 1:D
    n <- name_identifiers(channelNames) # c("1","2",...,"D", "1.2",...,"1.D","2.3",...,"D-1.D")
    nCpC1 <- nCpC2 <- list() # number of change points, indices of change points respectively
    nCpF <- nCpFA <- array(0, c(dim(cp[[1]])[1], 1)) # ...number change points frquencies, number change points frequencies auto 
    nCpT <- nCpTA <- array(0, c(dim(cp[[1]])[2], 1)) # ...number change points times, number change points time auto
    for (i in 1:length(cp)) {
        nCpC1[[i]] <- length(which(colSums(cp[[i]]) > 0)) # number of change points in each element of the list cp
        nCpC2[[i]] <- (which(colSums(cp[[i]]) > 0)) # indices of the change points 
        nCpF <- nCpF + rowSums(cp[[i]] > 0) # Sum of all the rows of the matrices of cp
        nCpT <- nCpT + (colSums(cp[[i]]) > 0) * 1 # ??... why put *1 ?sum of the columns 
        if (i <= D) { # if the quantity is related to the auto-spectrum
            nCpFA <- nCpFA + rowSums(cp[[i]] > 0)
            nCpTA <- nCpTA + (colSums(cp[[i]]) > 0) * 1
        }
    }
    nCp <- sum(nCpT) # total number of change points 
    nCpA <- sum(nCpTA) # total number of change points in autospectra 
    # naming the rows...
    rownames(nCpT) <- rownames(nCpTA) <- seq(1:ncol(cp[[1]])) * 
        windowLen
    rownames(nCpF) <- rownames(nCpFA) <- paste("Freq. band", 
        1:length(nCpF))
    # datafram containing the number of change points per channel/channel pair, and the location of each change points 
    tmp <- data.frame(Name = n, `Number change points` = unlist(nCpC1), 
        `Change-point location` = unlist(lapply(nCpC2, paste, 
            collapse = ", ")))
    
    out <- list(`Total # change points` = nCp, `Total # change points in autospectra` = nCpA, 
        `Change points over time` = nCpT, `Autospectral change points over time` = nCpTA, 
        `Change points over frequency bands` = nCpF, `Autospectral change points over frequency bands` = nCpFA, 
        `Change points per channel/channel pair` = tmp)
    if (plot) { # ...
        print("Plotting two charts")
        try(dev.new(noRStudioGD = TRUE), silent = TRUE)
        plot.new()
        plot(cumsum(nCpT), type = "l", ylab = "Cumulative # change points", 
            xlab = "Time", xaxt = "n")
        title("Change points in time")
        lines(cumsum(nCpTA), col = 2)
        lines(cumsum(nCpT) - cumsum(nCpTA), col = 4)
        axis(1, seq(0, length(nCpT), by = round_right(0.1 * 
            length(nCpT))), c(0, rownames(nCpT)[seq(round_right(0.1 * 
            length(nCpT)), length(nCpT), by = round_right(0.1 * 
            length(nCpT)))]))
        legend("topleft", c("Change points: total", "... in autospectra", 
            "... in cross-coherences"), lwd = 2, lty = c(1, 2, 
            4), col = c(1, 2, 4), bty = "o", bg = "white")
        try(dev.new(noRStudioGD = TRUE), silent = TRUE)
        plot.new()
        par(lend = "square")
        plot(nCpF, type = "h", lwd = 15, ylab = "Number of change points", 
            xlab = "Frequency band")
        title("Total change points over frequencies")
    }
    return(out)
}