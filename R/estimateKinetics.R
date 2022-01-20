#' Calculate kinetic rates for genes using multiple thread
#'
#' @param observedData a data frame with gene expression of m_u, m_total, p_u and p_total in sorted cells.
#' @param threadN an integer the number of thread for the calculation (default: 8)
#' @param numberofSplineAnchors an integer for the number of anchor points used for smooth the gene expression profiles along the cell cycle  (default: 20)
#'
#' @return a data frame with columns of gene, cc_time (cell cycle time), transcription rate (alpha), degradation rate (gamma).
#' @export
#'
#' @examples data_rates <- estimateKinetics(observedData = observed.data, threadN = 8, numberOfSplineAnchors = 20)
#'
#'
#'
#' calculate kinetic rates for multiple genes using multiple thread
estimateKinetics <-
  function(smoothedObservedData,
           threadN = 8,
           numberOfSplineAnchors = 20){

    if (!('foreach' %in% installed.packages())) {
      stop("Please install foreach")
    }
    if (!('doParallel' %in% installed.packages())) {
      stop("Please install doParallel")
    }
    library(foreach)
    library(doParallel)
    registerDoParallel(threadN)

    dataRates <-
      foreach (i = as.vector(unique(smoothedObservedData$gene)), .combine=rbind) %dopar% {
        tryCatch(smoothedObservedData %>%
                   dplyr::filter(gene == i) %>%
                   droplevels %>%
                   as.data.frame %>%
                   do(calculateKineticsForEachGene(., numberOfSplineAnchors = numberOfSplineAnchors)),
                 error=function(e) NULL)
      }
    return(dataRates)
  }





