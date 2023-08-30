## Eskrate

### Eskrate estimates time-dependent RNA kinetic rates using:
1. single cells with assigned times in biological process
2. newly synthesized (labeled) and pre-existing (unlabeled) RNAs gene expression data for single cells (e.g. SLAM-Drop-seq data). 

To use this package, first install it in R:

```{r}
devtools::install_github("rajewsky-lab/Eskrate")  
```

Here's a simple tutorial for estiamting the cell cycle time-dependent RNA kinetic rates with the test data provided:

#### Load required libraries

Install the required library if it is not installed yet.

```{r}
library(tibble)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(Rmpfr)
library(gam)
library(mgcv)
library(foreach)
library(doParallel)
library(Eskrate)
```

#### Input data

1. cell meta data: cell_info containing columns of cell_id and cc_time.

2. raw gene expression data: matrices of raw gene expression data of four RNA types containing columns of cell_id, gene, capture_rate, p_u, p_total, m_u, m_total.

```{r}
### cell meta data
cells_info <- readRDS("/data/cells_info.rds"))

### gene expression matrix
raw_labeled_mature <- readRDS("/data/raw_counts_labeled_mature.rds"))
raw_unlabeled_mature <- readRDS("/data/raw_counts_unlabeled_mature.rds"))
raw_labeled_precursor <- readRDS( "/data/raw_counts_labeled_precursor.rds"))
raw_unlabeled_precursor <- readRDS("/data/raw_counts_unlabeled_precursor.rds"))

```


```{r}
observed.data <- 
  getInput(rawCountsLabeledMature = raw_labeled_mature,
           rawCountsUnlabeledMature = raw_unlabeled_mature,
           rawCountsLabeledPrecursor = raw_labeled_precursor,
           rawCountsUnlabeledPrecursor = raw_unlabeled_precursor,
           cells = cells_info,
           labelingTime = 15)

### select two genes (ZNF711, KIF14) for demo
observed.data %<>% 
  filter(gene %in% c("ZNF711", "KIF14")) %>%
  droplevels
```

#### Estimate RNA kinetic rates of transcription and degradation

```{r}
data_rates <- 
  estimateKinetics(observedData = observed.data, 
                   threadN = 2)
```

#### Calculate predictions of gene expression using the estimated kinetic rates

```{r}
data_predicted <- 
  getPredictions(dataRates = data_rates, 
                 threadN=2)
```

#### Plot predicted gene expression and estimated rates

For example, plot PCNA gene expression and rates.

```{r plot_profile, fig.height=8, fig.width=10}
p1 <-
  plotMtotalAlphaGamma(df.prediction = predicted.data,   
                       gene.to.plot = "ZNF711",
                       color.alpha = "#D95F02",
                       color.gamma = "#1B9E77",
                       color.expression = "#666666",
                       show.phase.boundaries = FALSE,
                       line.size = 1,
                       font.size = 16)

p2 <-
  plotMtotalAlphaGamma(df.prediction = predicted.data,   
                       gene.to.plot = "KIF14",
                       color.alpha = "#D95F02",
                       color.gamma = "#1B9E77",
                       color.expression = "#666666",
                       show.phase.boundaries = FALSE,
                       line.size = 1,
                       font.size = 16)

plot_grid(p1, p2, ncol=2)
```

### Reference
Liu H*, Arsiè R*, Schwabe D, Schilling M, Minia I, Alles J, Boltengagen A, Kocks C, Falcke M, Friedman N, Landthaler M<sup>#</sup> Rajewsky N<sup>#</sup> (2023) [**SLAM‐DROP ‐seq reveals MRNA kinetic rates throughout the cell cycle.**](https://doi.org/10.15252/msb.202211427) Mol Syst Biol: e11427



