#' Calculate predicted gene expression based on estimated kinetic rates (only alpha and gamma)
#' using a simplified model.
#'
#' @param dataRates a data frame with estimated time-dependent transcription (alpha) and degradation (gamma) rates
#' @param threadN an integer the number of thread that used for calculation (default: 8)
#'
#' @return a data frame with predicted gene expression of m_u and m_total for genes in sorted cells.
#' @export
#'
#' @examples prediction <- getPredictions(dataRates, threadN = 8)
#'
getPredictions <-
  function(dataRates,
           threadN = 8){

    if (!('foreach' %in% installed.packages())) {
      stop("Please install foreach")
    }
    if (!('doParallel' %in% installed.packages())) {
      stop("Please install doParallel")
    }
    library(foreach)
    library(doParallel)
    registerDoParallel(threadN)

    prediction <-
      foreach(i = as.vector(unique(dataRates$gene)), .combine = rbind) %dopar% {
        tryCatch(dataRates %>%
                   dplyr::filter(gene==i) %>%
                   mutate(labeling_time = 15) %>%
                   do(getPredictionSimplifiedModelForEachGene(.)),
                 error=function(e) NULL)
      }
    return(prediction)
  }

