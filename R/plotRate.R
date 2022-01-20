#' Plot the predicted gene expression, transcription (alpha) and degradation (gamma) rate profiles over the cell cycle time.
#'
#' @param df.prediction a data frame of the predicted gene expression, transcription (alpha) and degradation (gamma) with gene in sorted cells.
#' @param gene.to.plot character the name of the gene to plot (default: "PCNA")
#' @param color.alpha character the color used for transcription rate (default: "#D95F02")
#' @param color.gamma character the color used for degradation rate (default: "#1B9E77")
#' @param color.expression character the color used for predicted expression (default: "#666666")
#' @param show.phase.boundaries logical weather to plot the cell cycle phase boundaries in the x axis (default: TRUE)
#' @param line.size numeric the line size for the profile plot (default: 1)
#' @param font.size numeric the font size for the plot (default: 16)
#'
#' @return
#' @export
#'
#' @examples plotMtotalAlphaGamma(df.prediction=df_predicted, gene.to.plot="PCNA")
#'
plotMtotalAlphaGamma <-
  function(df.prediction = prediction,   ### data frame of prediction for single gene to plot
           gene.to.plot = "PCNA",
           color.alpha = "#D95F02",
           color.gamma = "#1B9E77",
           color.expression = "#666666",
           show.phase.boundaries = TRUE,
           line.size = 1,
           font.size = 16){

    if (!('tibble' %in% installed.packages())) {
      stop("Please install tibble")
    }
    if (!('dplyr' %in% installed.packages())) {
      stop("Please install dplyr")
    }
    if (!('tidyr' %in% installed.packages())) {
      stop("Please install tidyr")
    }
    if (!('ggplot2' %in% installed.packages())) {
      stop("Please install ggplot2")
    }
    if (!('cowplot' %in% installed.packages())) {
      stop("Please install cowplot")
    }
    ### libraries
    library(tibble)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(cowplot)
    ### facet labels
    type.labs <- c("m_total"="total mature [CPM]",
                   "alpha" = "transcription [mol/min]",
                   "gamma" = "degradation [1/min]")

    #### parameters
    cell.cycle.duration <- 19.33
    if (show.phase.boundaries) {
      ## phase using [0,1] scales
      phase.proportions <- c(G1 = .487, S = .392, G2 = .093, M = .028)
      phases <-
        phase.proportions %>%
        enframe(name = "phase", value = "proportion") %>%
        mutate(phase = factor(phase, levels = phase),
               phase.index = 1:n(),
               duration = proportion * cell.cycle.duration,
               #duration = proportion

        ) %>%
        mutate(begin = lag(cumsum(duration), default = 0),
               end = lead(begin, default = cell.cycle.duration),
               mid = (begin + end) / 2)

      ### for plot cell cycle phase ticks
      tick <-
        data.frame(
          tick=sort(c(0, phases$mid, phases$end)),
          label = c("", "G1", "", "S", "", "G2", "", "M", ""),
          color = c("black", "white", "black", "white", "black", "white", "black", "white", "black"))

      p <-
        df.prediction %>%
        filter(gene == gene.to.plot) %>%
        dplyr::select(gene, cc_time, alpha, gamma, m_total) %>%
        pivot_longer(cols=-c(gene, cc_time), names_to="type", values_to="value") %>%
        mutate(data="predicted") %>%
        mutate(type=factor(type, levels = c("m_total", "alpha", "gamma"))) %>%
        ggplot(aes(cc_time/max(cc_time)*cell.cycle.duration, value, color=type)) %>%
        + geom_line(size=line.size) %>%
        + scale_color_manual(values = c(color.expression, color.alpha, color.gamma), guide = FALSE) %>%
        + scale_x_continuous(breaks = tick$tick, labels = tick$label) %>%
        + scale_y_continuous(limits = c(0, NA)) %>%
        + facet_wrap(~type, ncol=1, scale="free_y",
                     labeller = labeller(type = type.labs)) %>%
        + ggtitle(gene.to.plot) %>%
        + guides(linetype = guide_legend(title=gene.to.plot, override.aes = list(size=1))) %>%
        + theme_cowplot(font_size = font.size) %>%
        + theme(axis.ticks.x = element_line(colour = as.character(tick$color)),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                #legend.position = "top",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                panel.background = element_rect(colour = "black"))
      } else{
        p <-
          df.prediction %>%
          filter(gene == gene.to.plot) %>%
          dplyr::select(gene, cc_time, alpha, gamma, m_total) %>%
          pivot_longer(cols=-c(gene, cc_time), names_to="type", values_to="value") %>%
          mutate(data="predicted") %>%
          mutate(type=factor(type, levels = c("m_total", "alpha", "gamma"))) %>%
          ggplot(aes(cc_time/max(cc_time)*cell.cycle.duration, value, color=type)) %>%
          + geom_line(size=line.size) %>%
          + scale_color_manual(values = c(color.expression, color.alpha, color.gamma), guide = FALSE) %>%
          + xlab("cell cycle time [hour]") %>%
          + scale_y_continuous(limits = c(0, NA)) %>%
          + facet_wrap(~type, ncol=1, scale="free_y",
                       labeller = labeller(type = type.labs)) %>%
          + ggtitle(gene.to.plot) %>%
          + guides(linetype = guide_legend(title=gene.to.plot, override.aes = list(size=1))) %>%
          + theme_cowplot(font_size = font.size) %>%
          + theme(axis.title.y = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_blank(),
                  panel.background = element_rect(colour = "black"))
    }
    return(p)
  }

