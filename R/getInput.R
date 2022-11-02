################################################################
## Input data (combine cell_info with cc_time and counts data)
################################################################

#' Formatting input data by combine gene counts with cell cycle time
#'
#' @param rawCountsLabeledMature a data frame
#' count gene expression matrix of labeled mature RNAs, rownames are genes are colnames are cells.
#' @param rawCountsUnlabeledMature a data frame
#' count gene expression matrix of unlabeled mature RNAs, rownames are genes are colnames are cells.
#' @param rawCountsLabeledPrecursor a data frame
#' count gene expression matrix of unlabeled precursor RNAs, rownames are genes are colnames are cells.
#' @param rawCountsUnlabeledPrecursor a data frame
#' count gene expression matrix of unlabeled precursor RNAs, rownames are genes are colnames are cells.
#'
#' @param cellInfo a data frame with columns:
#'        cell_id (unique cell id)
#'        cc_Time (cell cycle time in minutes)
#'
#' @param labelingTime numeric the 4sU labeling time for the cells
#'
#' @return a data frame with columns:
#'    gene,
#'    cc_time,
#'    labeling_time,
#'    p_u,
#'    p_total,
#'    m_u,
#'    m_total,
#'    captureEfficiencyPerCell
#' @export
#'
#' @examples observedData <- getInput(rawCounts = raw_counts, cellInfo=cells_info)
#'
#'
getInput <-
  function(rawCountsLabeledMature,
           rawCountsUnlabeledMature,
           rawCountsLabeledPrecursor,
           rawCountsUnlabeledPrecursor,
           cells,
           labelingTime){
    if (!('dplyr' %in% installed.packages())) {
      stop("Please install dplyr")
      }
    if (!('tidyr' %in% installed.packages())) {
      stop("Please install tidyr")
    }    
    if (!('tibble' %in% installed.packages())) {
      stop("Please install tibble")
    }
    if (!('magrittr' %in% installed.packages())) {
      stop("Please install magrittr")
    }
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(magrittr)

    observed.data <-
      rawCountsLabeledMature %>%
      rownames_to_column("gene") %>%
      pivot_longer(cols=-gene, names_to="cell_id", values_to="m_l") %>%
      left_join(rawCountsUnlabeledMature %>%
                  rownames_to_column("gene") %>%
                  pivot_longer(cols=-gene, names_to="cell_id", values_to="m_u")) %>%
      left_join(rawCountsLabeledPrecursor %>%
                  rownames_to_column("gene") %>%
                  pivot_longer(cols=-gene, names_to="cell_id", values_to="p_l")) %>%
      left_join(rawCountsUnlabeledPrecursor %>%
                  rownames_to_column("gene") %>%
                  pivot_longer(cols=-gene, names_to="cell_id", values_to="p_u")) %>%
      left_join(cells) %>%
      arrange(gene, cc_time) %>%
      mutate(m_total = m_u + m_l,
             p_total = p_u + p_l) %>%
      mutate(labeling_time = labelingTime) %>%
      as.data.frame

    observed.data %<>%
      mutate(total = m_total + p_total) %>%
      group_by(cell_id) %>%
      summarise(total = sum(total)) %>%
      mutate(captureEfficiencyPerCell = total / 1e6) %>%
      select(-total) %>%
      left_join(observed.data) %>%
      arrange(gene, cc_time) %>%
      as.data.frame

    ## data info
    number.of.cells <- length(unique(observed.data$cc_time))
    cat("Number of cells: ", number.of.cells, "\n")

    number.of.genes <- length((unique(observed.data$gene)))
    cat("Number of genes: ", number.of.genes, "\n")

    return(observed.data)

  }


# raw.observed.data <- readRDS("../data/observed.data.nadia.15.rds")
#
# m_l <-
#   raw.observed.data %>%
#   select(cell_id, gene, m_l) %>%
#   pivot_wider(names_from = cell_id, values_from = m_l, values_fill = 0) %>%
#   as.data.frame
#
# m_u <-
#   raw.observed.data %>%
#   select(cell_id, gene, m_u) %>%
#   pivot_wider(names_from = cell_id, values_from = m_u,  values_fill = 0) %>%
#   as.data.frame
#
# p_l <-
#   raw.observed.data %>%
#   select(cell_id, gene, p_l) %>%
#   pivot_wider(names_from = cell_id, values_from = p_l, values_fill = 0) %>%
#   as.data.frame
#
# p_u <-
#   raw.observed.data %>%
#   select(cell_id, gene, p_u) %>%
#   pivot_wider(names_from = cell_id, values_from = p_u, values_fill = 0) %>%
#   as.data.frame
#
# rownames(m_l) <- m_l$gene
# rownames(m_u) <- m_u$gene
# rownames(p_l) <- p_l$gene
# rownames(p_u) <- p_u$gene
#
# m_l <- m_l[,-1]
# m_u <- m_u[,-1]
# p_l <- p_l[,-1]
# p_u <- p_u[,-1]
# saveRDS(m_l, "raw_counts_labeled_mature.rds")
# saveRDS(m_u, "raw_counts_unlabeled_mature.rds")
# saveRDS(p_l, "raw_counts_labeled_precursor.rds")
# saveRDS(p_u, "raw_counts_unlabeled_precursor.rds")
#
# rawCountsLabeledMature <- m_l
# rawCountsUnlabeledMature <- m_u
# rawCountsLabeledPrecursor <- p_l
# rawCountsUnlabeledPrecursor <- p_u
