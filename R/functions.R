#' Normalize gene expression count data to concentration
#'
#' @param rawData a data frame of raw counts of 'm_u', [m_l', 'm_total', 'P_l', p_u', 'p_total' as columns, and 'captureEfficiencyPerCell' is also provide as a column
#' @param umiPerCell numeric molecules per cell (default: 1e6)
#'
#' @return data frame with gene counts normalized to concentrations
#' @export
#'
#' @examples data.norm <- normalizeCountsToConcentration(rawData, umiPerCell = 1e6)
#'
normalizeCountsToConcentration <-
  function(rawData,
           umiPerCell = 1e6){

    if(is.null(rawData[["captureEfficiencyPerCell"]])){
      stop(paste("captureEfficiencyPerCell is not provided !!"))
    }
    cat("Normalized gene expression counts to concentrations. \n")
    rawData[["m_u"]] <- rawData[["m_u"]] / rawData[["captureEfficiencyPerCell"]] / umiPerCell
    rawData[["m_l"]] <- rawData[["m_l"]] / rawData[["captureEfficiencyPerCell"]] / umiPerCell
    rawData[["m_total"]] <- rawData[["m_total"]] / rawData[["captureEfficiencyPerCell"]] / umiPerCell
    rawData[["p_u"]] <- rawData[["p_u"]] / rawData[["captureEfficiencyPerCell"]] / umiPerCell
    rawData[["p_l"]] <- rawData[["p_l"]] / rawData[["captureEfficiencyPerCell"]] / umiPerCell
    rawData[["p_total"]] <- rawData[["p_total"]] / rawData[["captureEfficiencyPerCell"]] / umiPerCell
    return(rawData)
  }

#' Smooth gene expression by fitting penalized splines using general additive models.
#'
#' @param dataConcentration data frame of normalized gene expression for m_u,: genes on rows and sorted cells on columns
#' @param numberOfAnchorPoints integer the number of anchors used for the spline fitting (default: 20)
#' @param gamma numeric a number that will be passed to GAM model to adjust the over-fitting (defaut: 1, use 1.4 if you want to adjust overfit)
#' @param isInputConcentrations logic
#' @param ignoreOutliers logic (default: TRUE)
#' @param outlierIQRfactor numeric (default: 3)
#'
#' @return a data frame with smoothed values in columns and sorted cells order in rows
#' @export
#'
#' @examples data_spline <- smoothGeneProfileByPenalizedSpline(dataConcentration=observed.data["PCNA",], numberOfAnchorPoints = 20, gamma=1.4)
#'

smoothGeneProfileByPenalizedSpline <-
  function(dataConcentration,
           numberOfAnchorPoints = 20,
           gamma = 1.4,
           isInputConcentrations = TRUE,
           ignoreOutliers = FALSE,
           outlierIQRfactor = 3){

    cat("Smooth sorted gene expression using GAM. \n")
    library(mgcv)

    outlierTransPseudocount = 1
    concentration2CPM = function(con) con * 1e6
    outlierTrans = function(cpm) log(outlierTransPseudocount + cpm)
    ### generating spline for each RNA type
    for (i in unique(dataConcentration[["gene"]])){
      #### m_total
      ys.all <- dataConcentration[["m_total"]][dataConcentration[["gene"]]==i]
      xs.all <- seq(0, 1, length.out = length(ys.all))
      valid <- rep(TRUE, length(xs.all))
      if (ignoreOutliers) {
        cpm <- ys.all
        if (isInputConcentrations) cpm <- concentration2CPM(ys.all)
        cpm <- outlierTrans(cpm)
        cpm.median <- median(cpm)
        cpm.iqr <- quantile(cpm, .75) - quantile(cpm, .25)
        cpm.min <- cpm.median - outlierIQRfactor * cpm.iqr
        cpm.max <- cpm.median + outlierIQRfactor * cpm.iqr
        valid <- cpm >= cpm.min & cpm <= cpm.max
      }
      xs <- xs.all[valid]
      ys <- ys.all[valid]

      model <- gam(ys ~ s(xs, bs = "cc", k = numberOfAnchorPoints), gamma = gamma)
      if(ignoreOutliers){
        dataConcentration[["m_total_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
        dataConcentration[["m_total_outlier"]][dataConcentration[["gene"]]==i] <- !valid
      }
      if (!ignoreOutliers) {
        dataConcentration[["m_total_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
      }

      #### m_u
      ys.all <- dataConcentration[["m_u"]][dataConcentration[["gene"]]==i]
      xs.all <- seq(0, 1, length.out = length(ys.all))
      valid <- rep(TRUE, length(xs.all))
      if (ignoreOutliers) {
        cpm <- ys.all
        if (isInputConcentrations) cpm <- concentration2CPM(ys.all)
        cpm <- outlierTrans(cpm)
        cpm.median <- median(cpm)
        cpm.iqr <- quantile(cpm, .75) - quantile(cpm, .25)
        cpm.min <- cpm.median - outlierIQRfactor * cpm.iqr
        cpm.max <- cpm.median + outlierIQRfactor * cpm.iqr
        valid <- cpm >= cpm.min & cpm <= cpm.max
      }
      xs <- xs.all[valid]
      ys <- ys.all[valid]

      model <- gam(ys ~ s(xs, bs = "cc", k = numberOfAnchorPoints), gamma = gamma)

      if(ignoreOutliers){
        dataConcentration[["m_u_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
        dataConcentration[["m_u_outlier"]][dataConcentration[["gene"]]==i] <- !valid
      }
      if (!ignoreOutliers) {
        dataConcentration[["m_u_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
      }

      #### m_l
      ys.all <- dataConcentration[["m_total"]][dataConcentration[["gene"]]==i] - dataConcentration[["m_u"]][dataConcentration[["gene"]]==i]
      xs.all <- seq(0, 1, length.out = length(ys.all))
      valid <- rep(TRUE, length(xs.all))
      if (ignoreOutliers) {
        cpm <- ys.all
        if (isInputConcentrations) cpm <- concentration2CPM(ys.all)
        cpm <- outlierTrans(cpm)
        cpm.median <- median(cpm)
        cpm.iqr <- quantile(cpm, .75) - quantile(cpm, .25)
        cpm.min <- cpm.median - outlierIQRfactor * cpm.iqr
        cpm.max <- cpm.median + outlierIQRfactor * cpm.iqr
        valid <- cpm >= cpm.min & cpm <= cpm.max
      }
      xs <- xs.all[valid]
      ys <- ys.all[valid]

      model <- gam(ys ~ s(xs, bs = "cc", k = numberOfAnchorPoints), gamma = gamma)

      if(ignoreOutliers){
        dataConcentration[["m_l_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
        dataConcentration[["m_l_outlier"]][dataConcentration[["gene"]]==i] <- !valid
      }
      if (!ignoreOutliers) {
        dataConcentration[["m_l_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
      }

      #### p_u
      ys.all <- dataConcentration[["p_u"]][dataConcentration[["gene"]]==i]
      xs.all <- seq(0, 1, length.out = length(ys.all))
      valid <- rep(TRUE, length(xs.all))
      if (ignoreOutliers) {
        cpm <- ys.all
        if (isInputConcentrations) cpm <- concentration2CPM(ys.all)
        cpm <- outlierTrans(cpm)
        cpm.median <- median(cpm)
        cpm.iqr <- quantile(cpm, .75) - quantile(cpm, .25)
        cpm.min <- cpm.median - outlierIQRfactor * cpm.iqr
        cpm.max <- cpm.median + outlierIQRfactor * cpm.iqr
        valid <- cpm >= cpm.min & cpm <= cpm.max
      }
      xs <- xs.all[valid]
      ys <- ys.all[valid]

      model <- gam(ys ~ s(xs, bs = "cc", k = numberOfAnchorPoints), gamma = gamma)

      if(ignoreOutliers){
        dataConcentration[["p_u_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
        dataConcentration[["p_u_outlier"]][dataConcentration[["gene"]]==i] <- !valid
      }
      if (!ignoreOutliers) {
        dataConcentration[["p_u_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
      }

      #### p_total
      ys.all <- dataConcentration[["p_total"]][dataConcentration[["gene"]]==i]
      xs.all <- seq(0, 1, length.out = length(ys.all))
      valid <- rep(TRUE, length(xs.all))
      if (ignoreOutliers) {
        cpm <- ys.all
        if (isInputConcentrations) cpm <- concentration2CPM(ys.all)
        cpm <- outlierTrans(cpm)
        cpm.median <- median(cpm)
        cpm.iqr <- quantile(cpm, .75) - quantile(cpm, .25)
        cpm.min <- cpm.median - outlierIQRfactor * cpm.iqr
        cpm.max <- cpm.median + outlierIQRfactor * cpm.iqr
        valid <- cpm >= cpm.min & cpm <= cpm.max
      }
      xs <- xs.all[valid]
      ys <- ys.all[valid]

      model <- gam(ys ~ s(xs, bs = "cc", k = numberOfAnchorPoints), gamma = gamma)

      if(ignoreOutliers){
        dataConcentration[["p_total_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
        dataConcentration[["p_total_outlier"]][dataConcentration[["gene"]]==i] <- !valid
      }
      if (!ignoreOutliers) {
        dataConcentration[["p_total_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
      }
    }

    #### p_l
    ys.all <- dataConcentration[["p_l"]][dataConcentration[["gene"]]==i]
    xs.all <- seq(0, 1, length.out = length(ys.all))
    valid <- rep(TRUE, length(xs.all))
    if (ignoreOutliers) {
      cpm <- ys.all
      if (isInputConcentrations) cpm <- concentration2CPM(ys.all)
      cpm <- outlierTrans(cpm)
      cpm.median <- median(cpm)
      cpm.iqr <- quantile(cpm, .75) - quantile(cpm, .25)
      cpm.min <- cpm.median - outlierIQRfactor * cpm.iqr
      cpm.max <- cpm.median + outlierIQRfactor * cpm.iqr
      valid <- cpm >= cpm.min & cpm <= cpm.max
    }
    xs <- xs.all[valid]
    ys <- ys.all[valid]

    model <- gam(ys ~ s(xs, bs = "cc", k = numberOfAnchorPoints), gamma = gamma)

    if(ignoreOutliers){
      dataConcentration[["p_l_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
      dataConcentration[["p_l_outlier"]][dataConcentration[["gene"]]==i] <- !valid
    }
    if (!ignoreOutliers) {
      dataConcentration[["p_l_spline"]][dataConcentration[["gene"]]==i] <- predict(model, list(xs = xs.all))
    }

    return(dataConcentration)
  }

#'
#'  Calculate time-dependent alpha (transcription rate) and gamma (degradation rate) for each gene using observed gene expression of m_u, m_total
#'
#' @param observedCounts data frame
#' @param boolCalculateInMinutes logic statement (default: TRUE)
#' @param totalTimeOfOneCCInHours numeric total cell cycle time in hours (default: 19.33)
#' @param totalTimeOfOneCCInMinutes numeric total cell cycle time in minutes (default: 19.33*60)
#' @param numberOfSplineAnchors number of spline anchors used for smoothing profiles (default: 20)
#' @param boolSetBetaConstantBeforeAlphaEstimate logic statement (default: TRUE)
#' @param ignoreOutliers logic statement (default: TRUE)
#' @param boolUseCPM logic statement (default: TRUE)
#' @param smallestLabelingTime a number the smallest metabolic labeling time (default: 15)
#'
#' @return a data frame with kinetic rates of transcription (alpha) and degradation (gamma) for each gene at each cell cycle time
#' @export
#' @examples data_rate <- calculateKineticsForEachGene(observedCounts, smallestLabelingTime = 15)
#'
calculateKineticsForEachGene <-
  function(observedData,    ## smoothed observed raw counts for one single gene with capture rate for each cell
           boolCalculateInMinutes = TRUE,
           totalTimeOfOneCCInHours = 19.33,
           totalTimeOfOneCCInMinutes = 19.33*60,
           numberOfSplineAnchors = 20,  ### parse to spline function
           boolSetBetaConstantBeforeAlphaEstimate = TRUE,
           ignoreOutliers =TRUE,
           boolUseCPM = TRUE,
           smallestLabelingTime = NULL){

    smallestLabelingTime = min(observedData[['labeling_time']])

    if (boolCalculateInMinutes){
      totalTimeOfOneCC <- totalTimeOfOneCCInMinutes
    }else{
      totalTimeOfOneCC <- totalTimeOfOneCCInHours
    }
    fixedLabelingTimesHelp <- unique(observedData[,'labeling_time'])
    fixedLabelingTimesHelp <- fixedLabelingTimesHelp[!is.na(fixedLabelingTimesHelp)]
    fixedLabelingTimesHelp <- fixedLabelingTimesHelp[order(fixedLabelingTimesHelp)]
    #smallestLabelingTime <- min(fixedLabelingTimesHelp[fixedLabelingTimesHelp>0])]
    ## to choose the smallestLabelingTime
    tValuesPerOnePeriod <- observedData[,'cc_time']  ##cc_time

    numberOfPointsRequiredToReachLabelingTime <- which.min(abs(tValuesPerOnePeriod-smallestLabelingTime))-1

    #omega <- 2*pi/totalTimeOfOneCC
    omega <- 1

    #######################################
    ## upscale the spline to CPM
    #######################################

    if (boolUseCPM){
      true_umi = rep(10^6, nrow(observedData))
      observedData[["true_umi"]] <- true_umi  ### true umi is 1 M
      ###########################
      ## spline
      ###########################
      ## p_l
      splineplYValues <- (observedData[["p_total_spline"]] - observedData[["p_u_spline"]]) * true_umi
      splineplYValues <- pmax(splineplYValues,10^(-10))
      ## p_u
      splinepuYValues <- observedData[["p_u_spline"]] * true_umi
      splinepuYValues <- pmax(splinepuYValues,10^(-10))
      ## p_total
      splineptotalYValues <- observedData[["p_total_spline"]] * true_umi
      splineptotalYValues <- pmax(splineptotalYValues,10^(-10))
      ## m_l
      splinemlYValues <- observedData[["m_l_spline"]] * true_umi
      splinemlYValues <- pmax(splinemlYValues,10^(-10))
      ## m_u
      splinemuYValues <- observedData[["m_u_spline"]] * true_umi
      splinemuYValues <- pmax(splinemuYValues,10^(-10))
      ## m total
      splinemtotalYValues <- observedData[["m_total_spline"]]* true_umi
      splinemtotalYValues <- pmax(splinemtotalYValues,10^(-10))
    }

    ###############################################
    ###############################################
    # cellsCurrentPhase <- (numberOfPointsRequiredToReachLabelingTime+1):length(tValuesPerOnePeriod)
    # cellsInitialPhase <- 1:(length(tValuesPerOnePeriod)-numberOfPointsRequiredToReachLabelingTime)

    ###############################################################################################
    ## 2020-05-11 haiyue
    ## extend the phase for 15 minutes
    ###############################################################################################
    cellsInitialPhase <- 1:(length(tValuesPerOnePeriod))
    cellsCurrentPhase <- c((numberOfPointsRequiredToReachLabelingTime+1):length(tValuesPerOnePeriod), 1: numberOfPointsRequiredToReachLabelingTime)
    ###########################################
    ### [function 25] m - m_u is used for m_l
    ###########################################
    gamma <- -1/(omega*smallestLabelingTime)*log(splinemuYValues[cellsCurrentPhase]/splinemtotalYValues[cellsInitialPhase])
    alpha <- gamma * (splinemtotalYValues[cellsCurrentPhase]-splinemuYValues[cellsCurrentPhase]) * splinemtotalYValues[cellsInitialPhase] / (splinemtotalYValues[cellsInitialPhase]-splinemuYValues[cellsCurrentPhase])
    #alpha <- gamma * splinemlYValues[cellsCurrentPhase] * splinemtotalYValues[cellsInitialPhase] / (splinemtotalYValues[cellsInitialPhase]-splinemuYValues[cellsCurrentPhase])
    alpha[is.na(alpha)] <- 10^(-10)
    beta <- alpha / (splineplYValues[cellsCurrentPhase])

    ################################
    ## take mean beta as constant
    ################################
    if (boolSetBetaConstantBeforeAlphaEstimate){
      beta <- rep(mean(beta), length(beta))
    }

    # alpha <- append(alpha, rep(alpha[length(alpha)], numberOfPointsRequiredToReachLabelingTime))
    # beta <- append(beta, rep(beta[length(beta)], numberOfPointsRequiredToReachLabelingTime))
    # gamma <- append(gamma, rep(gamma[length(gamma)], numberOfPointsRequiredToReachLabelingTime))

    gamma <- pmax(gamma,10^(-10))
    alpha <- pmax(alpha,10^(-10))
    beta <- pmax(beta,10^(-10))

    ##########################################
    ## return data frame with parameters
    ##########################################
    rateData <-
      data.frame('gene' = observedData[,'gene'],
                 'cc_time' = observedData[,'cc_time'],
                 'labeling_time' = observedData[, 'labeling_time'],
                 'alpha' = alpha,
                 'beta' = beta,
                 'gamma' = gamma,
                 'm_total' = observedData[,'m_total'], ### observed raw counts
                 'm_l' = observedData[, 'm_l'],
                 'm_u' = observedData[,'m_u'],
                 'p_total' = observedData[,'p_total'],
                 'p_l' = observedData[,'p_l'],
                 'p_u' = observedData[, 'p_u'],
                 'm_total_spline' = splinemtotalYValues, ###spline of the upscale concentration [CPM]
                 'm_l_spline' = splinemlYValues,
                 'm_u_spline' = splinemuYValues,
                 'p_total_spline' = splineptotalYValues,
                 'p_l_spline' = splineplYValues,
                 'p_u_spline' = splinepuYValues,
                 'captureEfficiencyPerCell' = observedData[['captureEfficiencyPerCell']])

    return(rateData)
  }

#'
#'
#' Get predictions from estimated kinetic rates using simplified model
#'
#' @param parameterDF a data frame with transcription (alpha) and degradation (gamma) rate for each gene at each cell cycle time
#' @param boolCheckedInputCCTimeHasNoDuplicates logical statement (default: TRUE)
#' @param boolCalculateInMinutes logical statement (default: TRUE)
#' @param totalTimeOfOneCCInHours numeric total cell cycle time in hours (default: 19.33)
#' @param totalTimeOfOneCCInMinutes numeric total cell cycle time in minutes (default: 19.33*60)
#' @param lengthToSimulatePerRoundInHours numeric (default: 5)
#' @param lengthToSimulatePerRoundInMinutes numeric (default: 5*60)
#' @param toleranceMature numeric (default: 1)
#' @param maximumNumberOfPeriodsToSimulate numeric (default: 10)
#' @param bytePrecisionForLargerNumbersToUse numeric (default: 60)
#'
#' @return a data frame with predicted gene expression for both labeled and unlabeled mature RNAs at each cell cycle time.
#' @export
#'
#' @examples data_predicted <- getPredictionSimplifiedModelForEachGene(parameterDF = data_rate)
getPredictionSimplifiedModelForEachGene <-
  function(parameterDF,
           boolCheckedInputCCTimeHasNoDuplicates = FALSE,
           boolCalculateInMinutes = TRUE,
           totalTimeOfOneCCInHours = 19.33,
           totalTimeOfOneCCInMinutes = 19.33*60,
           lengthToSimulatePerRoundInHours = 5,
           lengthToSimulatePerRoundInMinutes = 5*60,
           #tolerancePrecursor = 0.1,
           toleranceMature = 1,
           maximumNumberOfPeriodsToSimulate = 10,
           bytePrecisionForLargerNumbersToUse = 60){

    if (!('Rmpfr' %in% installed.packages())) {
      stop("Please install Rmpfr")
    }

    ### package for large numbers
    library(Rmpfr)

    getDefiniteIntegralOfFunction <- function(xVector,
                                              yVector){
      diffX <- xVector[2]-xVector[1]
      trapezoidalIntegrals <- (yVector[1:(length(yVector)-1)]+yVector[2:length(yVector)])/2
      return(cumsum(append(0, trapezoidalIntegrals*diffX)))
    }

    if (boolCalculateInMinutes){
      totalTimeOfOneCC <- totalTimeOfOneCCInMinutes
      lengthToSimulatePerRound <- lengthToSimulatePerRoundInMinutes
    }else{
      totalTimeOfOneCC <- totalTimeOfOneCCInHours
      lengthToSimulatePerRound <- lengthToSimulatePerRoundInHours
    }

    if (!boolCheckedInputCCTimeHasNoDuplicates){
      if (length(unique(parameterDF[,'cc_time']))==length(parameterDF[,'cc_time'])){
        tValuesPerOnePeriod <- parameterDF[,'cc_time']
        numberOfPointsPerOnePeriod <- length(tValuesPerOnePeriod)

        alphaThisGene <- parameterDF[,'alpha']
        # betaThisGene <- parameterDF[,'beta']
        gammaThisGene <- parameterDF[,'gamma']
      }else{
        tValuesPerOnePeriod <- unique(parameterDF[,'cc_time'])
        numberOfPointsPerOnePeriod <- length(tValuesPerOnePeriod)

        alphaThisGene <- parameterDF[match(unique(parameterDF[,'cc_time']),parameterDF[,'cc_time']),'alpha']
        # betaThisGene <- parameterDF[match(unique(parameterDF[,'cc_time']),parameterDF[,'cc_time']),'beta']
        gammaThisGene <- parameterDF[match(unique(parameterDF[,'cc_time']),parameterDF[,'cc_time']),'gamma']
      }
    }else{
      tValuesPerOnePeriod <- parameterDF[,'cc_time']
      numberOfPointsPerOnePeriod <- length(tValuesPerOnePeriod)

      alphaThisGene <- parameterDF[,'alpha']
      betaThisGene <- parameterDF[,'beta']
      gammaThisGene <- parameterDF[,'gamma']
    }

    #error checks for parameters
    if (any(c(any(is.na(alphaThisGene)), any(alphaThisGene<0)), any(is.na(gammaThisGene)), any(gammaThisGene<0))){
      if (any(is.na(alphaThisGene))){
        warning('alpha has NA values')
      }
      if (any(alphaThisGene<0)){
        warning('alpha has negative values')
      }
      # if (any(is.na(betaThisGene))){
      #   warning('beta has NA values')
      # }
      # if (any(betaThisGene<0)){
      #   warning('beta has negative values')
      # }
      if (any(is.na(gammaThisGene))){
        warning('gamma has NA values')
      }
      if (any(gammaThisGene<0)){
        warning('gamma has negative values')
      }
      parameterDF <- parameterDF[!is.na(parameterDF[,'labeling_time']),]
      dfForValues <- data.frame('gene' = parameterDF[, 'gene'],
                                'cc_time' = parameterDF[,'cc_time'],
                                'labeling_time' = parameterDF[,'labeling_time'],
                                'alpha' = parameterDF[,'alpha'],
                                'beta' = rep(NA, dim(parameterDF)[1]),
                                'gamma' = parameterDF[,'gamma'],
                                #'p_u' = rep(NA, dim(parameterDF)[1]),
                                #'p_total' = rep(NA, dim(parameterDF)[1]),
                                'm_u' = rep(NA, dim(parameterDF)[1]),
                                'm_total' = rep(NA, dim(parameterDF)[1]))
      return(dfForValues)
    }


    alphat <- function(t){
      indicesToSelect <- which(tValuesPerOnePeriod%in%t)
      return(alphaThisGene[indicesToSelect])
    }
    betat <- function(t){
      indicesToSelect <- which(tValuesPerOnePeriod%in%t)
      return(betaThisGene[indicesToSelect])
    }
    gammat <- function(t){
      indicesToSelect <- which(tValuesPerOnePeriod%in%t)
      return(gammaThisGene[indicesToSelect])
    }

    ####Part 1: Simulate cell cycle periods until we obtain convergent trajectory###

    # pt <- function(p0, integralBeta, integralf1Helper){
    #   return(p0 * exp(-integralBeta) + integralf1Helper)
    # }
    mt <- function(m0, integralGamma, integralf1Helper){
      return(m0 * exp(-integralGamma) + integralf1Helper)
    }

    ###initial values
    # p_0Save <- alphat(0)/betat(0)
    m_0Save <- alphat(0)/gammat(0)
    # if (is.na(p_0Save)|(p_0Save<0)|is.na(m_0Save)|(m_0Save<0)){
    if (is.na(m_0Save)|(m_0Save<0)){
      warning('start value p0 or m0 is NA or negative')
      parameterDF <- parameterDF[!is.na(parameterDF[,'labeling_time']),]
      dfForValues <- data.frame('cc_time' = parameterDF[,'cc_time'],
                                'labeling_time' = parameterDF[,'labeling_time'],
                                'alpha' = parameterDF[,'alpha'],
                                'beta' = rep(NA, dim(parameterDF)[1]),
                                'gamma' = parameterDF[,'gamma'],
                                #'p_u' = rep(NA, dim(parameterDF)[1]),
                                #'p_total' = rep(NA, dim(parameterDF)[1]),
                                'm_u' = rep(NA, dim(parameterDF)[1]),
                                'm_total' = rep(NA, dim(parameterDF)[1]))
      return(dfForValues)
    }

    #variables for saving time courses
    # ptValue <- NULL
    mtValue <- NULL

    #simulate one cell cycle period at a time
    counterNumberOfPeriods <- 0
    boolConvergenceAchieved <- FALSE
    boolMaximumIterationsAchieved <- FALSE
    boolStopSimulationDueToUnfeasibleValues <- FALSE
    while((!boolConvergenceAchieved&!boolMaximumIterationsAchieved)|(counterNumberOfPeriods<2)){
      counterNumberOfPeriods <- counterNumberOfPeriods+1
      #if not first iteration: get values from end of previous period (see end of loop)
      # p_0 <- p_0Save
      m_0 <- m_0Save
      for (j in 1:ceiling((totalTimeOfOneCC)/lengthToSimulatePerRound)){
        if (j>1){
          startIndex <- endIndex  #min(which(tValuesPerOnePeriod>=(j-1)*lengthToSimulatePerRound))-1
        }else{
          startIndex <- 1   #min(which(tValuesPerOnePeriod>=(j-1)*lengthToSimulatePerRound))
        }
        if (j==ceiling((totalTimeOfOneCC)/lengthToSimulatePerRound)){
          endIndex <- length(tValuesPerOnePeriod)
        }else{
          endIndex <- max(which(tValuesPerOnePeriod<j*lengthToSimulatePerRound))
        }

        ##########################
        # print(j)
        # cat(startIndex, ":", endIndex, "\n")
        ############################

        #get window of time we simulate
        indicesFromPeriodFunctionsToSelectThisWindow <- startIndex:endIndex
        tValuesThisWindow <- tValuesPerOnePeriod[indicesFromPeriodFunctionsToSelectThisWindow]

        #calculate alpha, beta, gamma and their antiderivative
        # yValuesBeta <- betat(tValuesThisWindow)
        # BFunctionValues <- getDefiniteIntegralOfFunction(xVector = tValuesThisWindow,
        #                                                  yVector = yValuesBeta)
        # if (any(is.na(BFunctionValues))|any(is.infinite(BFunctionValues))){
        #   boolStopSimulationDueToUnfeasibleValues <- TRUE
        # }
        # if (any(BFunctionValues>700)){
        #   warning('Beta values are extremely large, potentially producing Inf downstream. Try decreasing parameter *lengthToSimulatePerRoundInMinutes*. This will increase runtime but reduce size of values handled by the algorithm.')
        # }

        yValuesGamma <- gammat(tValuesThisWindow)
        CFunctionValues <- getDefiniteIntegralOfFunction(xVector = tValuesThisWindow,
                                                         yVector = yValuesGamma)
        if (any(is.na(CFunctionValues))|any(is.infinite(CFunctionValues))){
          boolStopSimulationDueToUnfeasibleValues <- TRUE
        }
        if (any(CFunctionValues>700)){
          warning('Gamma values are extremely large, potentially producing Inf downstream. Try decreasing parameter *lengthToSimulatePerRoundInMinutes*. This will increase runtime but reduce size of values handled by the algorithm.')
        }

        # f1HelperFunctionValues <- alphat(tValuesThisWindow)*exp(BFunctionValues)
        # if (any(is.na(f1HelperFunctionValues))|any(is.infinite(f1HelperFunctionValues))){
        #   boolStopSimulationDueToUnfeasibleValues <- TRUE
        # }
        # F1FunctionValues <- getDefiniteIntegralOfFunction(xVector = tValuesThisWindow,
        #                                                   yVector = f1HelperFunctionValues)
        # if (any(is.na(F1FunctionValues))|any(is.infinite(F1FunctionValues))){
        #   boolStopSimulationDueToUnfeasibleValues <- TRUE
        # }
        f1HelperFunctionValues <- alphat(tValuesThisWindow)*exp(CFunctionValues)
        if (any(is.na(f1HelperFunctionValues))|any(is.infinite(f1HelperFunctionValues))){
          boolStopSimulationDueToUnfeasibleValues <- TRUE
        }
        F1FunctionValues <- getDefiniteIntegralOfFunction(xVector = tValuesThisWindow,
                                                          yVector = f1HelperFunctionValues)
        if (any(is.na(F1FunctionValues))|any(is.infinite(F1FunctionValues))){
          boolStopSimulationDueToUnfeasibleValues <- TRUE
        }

        ### haiyue: total precursor [function:12]
        # ptValueHelp <- pt(p0 = p_0,
        #                   integralBeta = (BFunctionValues-BFunctionValues[1]),
        #                   integralf1Helper = exp(-BFunctionValues)*(F1FunctionValues-F1FunctionValues[1]))
        # if (any(is.na(ptValueHelp))|any(is.infinite(ptValueHelp))){
        #   boolStopSimulationDueToUnfeasibleValues <- TRUE
        # }



        # f2HelperFunctionValues <- yValuesBeta*ptValueHelp*exp(CFunctionValues)
        # if (any(is.na(f2HelperFunctionValues))|any(is.infinite(f2HelperFunctionValues))){
        #   boolStopSimulationDueToUnfeasibleValues <- TRUE
        # }
        # F2FunctionValues <- getDefiniteIntegralOfFunction(xVector = tValuesThisWindow,
        #                                                   yVector = f2HelperFunctionValues)
        # if (any(is.na(F2FunctionValues))|any(is.infinite(F2FunctionValues))){
        #   boolStopSimulationDueToUnfeasibleValues <- TRUE
        # }
        ### haiyue: total mature [function 13]
        mtValueHelp <- mt(m0 = m_0,
                          integralGamma = (CFunctionValues-CFunctionValues[1]),
                          integralf1Helper = exp(-CFunctionValues)*(F1FunctionValues-F1FunctionValues[1]))  ###?
        if (any(is.na(mtValueHelp))|any(is.infinite(mtValueHelp))){
          boolStopSimulationDueToUnfeasibleValues <- TRUE
        }

        if (boolStopSimulationDueToUnfeasibleValues){
          warning('Unfeasible values were produced during simulation of multiple periods.')
          parameterDF <- parameterDF[!is.na(parameterDF[,'labeling_time']),]
          dfForValues <- data.frame('cc_time' = parameterDF[,'cc_time'],
                                    'labeling_time' = parameterDF[,'labeling_time'],
                                    'alpha' = parameterDF[,'alpha'],
                                    'beta' = rep(NA, dim(parameterDF)[1]),
                                    'gamma' = parameterDF[,'gamma'],
                                    # 'p_u' = rep(NA, dim(parameterDF)[1]),
                                    # 'p_total' = rep(NA, dim(parameterDF)[1]),
                                    'm_u' = rep(NA, dim(parameterDF)[1]),
                                    'm_total' = rep(NA, dim(parameterDF)[1]))
          return(dfForValues)
        }

        if ((j==1)&(counterNumberOfPeriods==1)){
          # ptValue <- append(ptValue, ptValueHelp)
          mtValue <- append(mtValue, mtValueHelp)
        }else{
          # ptValue <- append(ptValue, ptValueHelp[2:length(ptValueHelp)])
          mtValue <- append(mtValue, mtValueHelp[2:length(mtValueHelp)])
        }
        # p_0 <- ptValueHelp[length(ptValueHelp)]
        m_0 <- mtValueHelp[length(mtValueHelp)]
      }
      #### 2020-09-29 haiyue: do not assume p or m doubling at the end of the cell cycle
      # if (((abs(ptValueHelp[length(ptValueHelp)]/1-p_0Save)<tolerancePrecursor)&abs(mtValueHelp[length(mtValueHelp)]/1-m_0Save)<toleranceMature)|(counterNumberOfPeriods>=maximumNumberOfPeriodsToSimulate)){
      #   if ((abs(ptValueHelp[length(ptValueHelp)]/1-p_0Save)<tolerancePrecursor)&abs(mtValueHelp[length(mtValueHelp)]/1-m_0Save)<toleranceMature){
      if ((abs(mtValueHelp[length(mtValueHelp)]/1-m_0Save)<toleranceMature)|(counterNumberOfPeriods>=maximumNumberOfPeriodsToSimulate)){
        if (abs(mtValueHelp[length(mtValueHelp)]/1-m_0Save)<toleranceMature){
          boolConvergenceAchieved <- TRUE
        }else{
          boolMaximumIterationsAchieved <- TRUE
          warning(paste('convergence not achieved; trajectory after ', counterNumberOfPeriods, ' periods interpreted as converged', sep = ''))
        }
      }

      ## update the initial values by dividing the last values by two (cell division)
      ## 2020-05-07 no doubling at the end of the cell cycle since we use concentration and CPM (divide by 1)
      # p_0SaveBefore <- p_0Save
      m_0SaveBefore <- m_0Save
      # p_0Save <- ptValueHelp[length(ptValueHelp)]
      m_0Save <- mtValueHelp[length(mtValueHelp)]
    }

    indicesOfLastPeriod <- (length(mtValue)-length(tValuesPerOnePeriod)+1):length(mtValue)
    # totalPrecursorLastPeriod <- ptValue[indicesOfLastPeriod]
    # totalPrecursorLastPeriod[1] <- totalPrecursorLastPeriod[1]
    totalMatureLastPeriod <- mtValue[indicesOfLastPeriod]
    totalMatureLastPeriod[1] <- totalMatureLastPeriod[1]
    alphaLastPeriod <- alphat(tValuesPerOnePeriod)
    # betaLastPeriod <- betat(tValuesPerOnePeriod)
    gammaLastPeriod <- gammat(tValuesPerOnePeriod)

    ####Part 2: Get unlabeled true counts###

    dfForValues <- NULL
    fixedLabelingTimesHelp <- unique(parameterDF[,'labeling_time'])
    fixedLabelingTimesHelp <- fixedLabelingTimesHelp[!is.na(fixedLabelingTimesHelp)]
    fixedLabelingTimesHelp <- fixedLabelingTimesHelp[order(fixedLabelingTimesHelp)]
    for (j in 1:length(fixedLabelingTimesHelp)){

      fixedLabelingTimeThisLoop <- fixedLabelingTimesHelp[j]

      numberOfPointsRequiredToReachLabelingTime <- which.min(abs(tValuesPerOnePeriod-fixedLabelingTimeThisLoop))-1

      if (numberOfPointsRequiredToReachLabelingTime>=1){
        tValuesPerOnePeriodPlusLabelingTime <- c(tValuesPerOnePeriod, seq(from = tValuesPerOnePeriod[length(tValuesPerOnePeriod)], to = tValuesPerOnePeriod[length(tValuesPerOnePeriod)]+tValuesPerOnePeriod[numberOfPointsRequiredToReachLabelingTime], by = tValuesPerOnePeriod[2]-tValuesPerOnePeriod[1]))
      }else{
        tValuesPerOnePeriodPlusLabelingTime <- tValuesPerOnePeriod
      }

      cellsToSelectThisLabelingTime <- which(parameterDF[,'labeling_time']==fixedLabelingTimeThisLoop)
      if (any(is.na(cellsToSelectThisLabelingTime))){
        warning('labeling time has NA values')
        parameterDF <- parameterDF[!is.na(parameterDF[,'labeling_time']),]
        dfForValues <- data.frame('cc_time' = parameterDF[,'cc_time'],
                                  'labeling_time' = parameterDF[,'labeling_time'],
                                  'alpha' = parameterDF[,'alpha'],
                                  'beta' = rep(NA, dim(parameterDF)[1]),
                                  'gamma' = parameterDF[,'gamma'],
                                  # 'p_u' = rep(NA, dim(parameterDF)[1]),
                                  # 'p_total' = rep(NA, dim(parameterDF)[1]),
                                  'm_u' = rep(NA, dim(parameterDF)[1]),
                                  'm_total' = rep(NA, dim(parameterDF)[1]))
        return(dfForValues)
      }
      if (length(cellsToSelectThisLabelingTime)==0){
        warning(paste('the labeling time ', fixedLabelingTimeThisLoop, ' is not present in parameters', sep = ''))
      }

      associatedtValuesThisLoop <- parameterDF[cellsToSelectThisLabelingTime,'cc_time']

      # yValuesBeta <- mpfr(betat(tValuesPerOnePeriod),bytePrecisionForLargerNumbersToUse)
      # if (length(tValuesPerOnePeriodPlusLabelingTime)>length(tValuesPerOnePeriod)){
      #   yValuesBeta <- c(yValuesBeta, betat(tValuesPerOnePeriod[1:numberOfPointsRequiredToReachLabelingTime]))
      # }
      # BFunctionValues <- mpfr(getDefiniteIntegralOfFunction(xVector = tValuesPerOnePeriodPlusLabelingTime,
      #                                                       yVector = as.numeric(yValuesBeta)),bytePrecisionForLargerNumbersToUse)
      # if (any(is.na(BFunctionValues))|any(is.infinite(BFunctionValues))){
      #   boolStopSimulationDueToUnfeasibleValues <- TRUE
      # }

      yValuesGamma <- mpfr(gammat(tValuesPerOnePeriod),bytePrecisionForLargerNumbersToUse)
      if (length(tValuesPerOnePeriodPlusLabelingTime)>length(tValuesPerOnePeriod)){
        yValuesGamma <- c(yValuesGamma, gammat(tValuesPerOnePeriod[1:numberOfPointsRequiredToReachLabelingTime]))
      }
      CFunctionValues <- mpfr(getDefiniteIntegralOfFunction(xVector = tValuesPerOnePeriodPlusLabelingTime,
                                                            yVector = as.numeric(yValuesGamma)),bytePrecisionForLargerNumbersToUse)
      if (any(is.na(CFunctionValues))|any(is.infinite(CFunctionValues))){
        boolStopSimulationDueToUnfeasibleValues <- TRUE
      }
      #### haiyue: function[40]
      #    p_u <- totalPrecursorLastPeriod[1:(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime)]*exp(-betat(0)*fixedLabelingTimeThisLoop)
      # p_u <- totalPrecursorLastPeriod[1:(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime)]*exp(-(BFunctionValues[(numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]-BFunctionValues[1:(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime)]))
      # if (numberOfPointsRequiredToReachLabelingTime>=1){
      #   #      p_uPastCellDivision <- totalPrecursorLastPeriod[(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]*exp(-betat(0)*fixedLabelingTimeThisLoop)/2
      #   p_uPastCellDivision <- totalPrecursorLastPeriod[(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]*exp(-(BFunctionValues[length(totalPrecursorLastPeriod)]-BFunctionValues[(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]))/1*exp(-BFunctionValues[1:numberOfPointsRequiredToReachLabelingTime])
      # }else{
      #   p_uPastCellDivision <- NULL
      # }
      # p_u <- append(as.numeric(p_uPastCellDivision), as.numeric(p_u))
      # if (any(is.na(p_u))|any(is.infinite(p_u))){
      #   boolStopSimulationDueToUnfeasibleValues <- TRUE
      # }

      # f3HelperFunctionValues <- yValuesBeta*exp(CFunctionValues-BFunctionValues)
      # if (any(is.na(f3HelperFunctionValues))|any(is.infinite(f3HelperFunctionValues))){
      #   boolStopSimulationDueToUnfeasibleValues <- TRUE
      # }
      # F3FunctionValues <- mpfr(getDefiniteIntegralOfFunction(xVector = tValuesPerOnePeriodPlusLabelingTime,
      #                                                        yVector = as.numeric(f3HelperFunctionValues)),bytePrecisionForLargerNumbersToUse)
      # if (any(is.na(F3FunctionValues))|any(is.infinite(F3FunctionValues))){
      #   boolStopSimulationDueToUnfeasibleValues <- TRUE
      # }
      ### haiyue function[41]
      # m_u <- totalPrecursorLastPeriod[1:(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime)]*betat(0)/(gammat(0)-betat(0))*(exp(-betat(0)*fixedLabelingTimeThisLoop)-exp(-gammat(0)*fixedLabelingTimeThisLoop)) + totalMatureLastPeriod[1:(length(totalMatureLastPeriod)-numberOfPointsRequiredToReachLabelingTime)]*exp(-gammat(0)*fixedLabelingTimeThisLoop)
      # m_u <- totalPrecursorLastPeriod[1:(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime)]*exp(-CFunctionValues[(numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]+BFunctionValues[1:(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime)])*(F3FunctionValues[(numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]-F3FunctionValues[1:(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime)]) + totalMatureLastPeriod[1:(length(totalMatureLastPeriod)-numberOfPointsRequiredToReachLabelingTime)]*exp(-CFunctionValues[(numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]+CFunctionValues[1:(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime)])
      # if (numberOfPointsRequiredToReachLabelingTime>=1){
      #   #      m_uPastCellDivision <- (totalPrecursorLastPeriod[(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]*betat(0)/(gammat(0)-betat(0))*(exp(-betat(0)*fixedLabelingTimeThisLoop)-exp(-gammat(0)*fixedLabelingTimeThisLoop)) + totalMatureLastPeriod[(length(totalMatureLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalMatureLastPeriod)]*exp(-gammat(0)*fixedLabelingTimeThisLoop))/2
      #   m_uPastCellDivision <- totalPrecursorLastPeriod[(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]*exp(-(CFunctionValues[1:numberOfPointsRequiredToReachLabelingTime]+CFunctionValues[length(totalPrecursorLastPeriod)])+BFunctionValues[(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)])*(F3FunctionValues[(length(totalPrecursorLastPeriod)+1):(length(totalPrecursorLastPeriod)+numberOfPointsRequiredToReachLabelingTime)]-F3FunctionValues[(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)])/1 + totalMatureLastPeriod[(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]*exp(-(CFunctionValues[1:numberOfPointsRequiredToReachLabelingTime]+CFunctionValues[length(totalPrecursorLastPeriod)])+CFunctionValues[(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)])/1
      # }else{
      #   m_uPastCellDivision <- NULL
      # }
      # m_u <- append(as.numeric(m_uPastCellDivision), as.numeric(m_u))
      # if (any(is.na(m_u))|any(is.infinite(m_u))){
      #   boolStopSimulationDueToUnfeasibleValues <- TRUE
      # }

      m_u <- totalMatureLastPeriod[1:(length(totalMatureLastPeriod)-numberOfPointsRequiredToReachLabelingTime)]*exp(-(CFunctionValues[(numberOfPointsRequiredToReachLabelingTime+1):length(totalMatureLastPeriod)]-CFunctionValues[1:(length(totalMatureLastPeriod)-numberOfPointsRequiredToReachLabelingTime)]))
      if (numberOfPointsRequiredToReachLabelingTime>=1){
        #      p_uPastCellDivision <- totalPrecursorLastPeriod[(length(totalPrecursorLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalPrecursorLastPeriod)]*exp(-betat(0)*fixedLabelingTimeThisLoop)/2
        m_uPastCellDivision <- totalMatureLastPeriod[(length(totalMatureLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalMatureLastPeriod)]*exp(-(CFunctionValues[length(totalMatureLastPeriod)]-CFunctionValues[(length(totalMatureLastPeriod)-numberOfPointsRequiredToReachLabelingTime+1):length(totalMatureLastPeriod)]))/1*exp(-CFunctionValues[1:numberOfPointsRequiredToReachLabelingTime])
      }else{
        m_uPastCellDivision <- NULL
      }
      m_u <- append(as.numeric(m_uPastCellDivision), as.numeric(m_u))
      if (any(is.na(m_u))|any(is.infinite(m_u))){
        boolStopSimulationDueToUnfeasibleValues <- TRUE
      }

      if (boolStopSimulationDueToUnfeasibleValues){
        warning('Unfeasible values were produced during generation of p_u and m_u.')
        parameterDF <- parameterDF[!is.na(parameterDF[,'labeling_time']),]
        dfForValues <- data.frame('gene' = parameterDF[,'gene'],
                                  'cc_time' = parameterDF[,'cc_time'],
                                  'labeling_time' = parameterDF[,'labeling_time'],
                                  'alpha' = parameterDF[,'alpha'],
                                  'beta' = rep(NA, dim(parameterDF)[1]),
                                  'gamma' = parameterDF[,'gamma'],
                                  #'p_u' = rep(NA, dim(parameterDF)[1]),
                                  #'p_total' = rep(NA, dim(parameterDF)[1]),
                                  'm_u' = rep(NA, dim(parameterDF)[1]),
                                  'm_total' = rep(NA, dim(parameterDF)[1]))
        return(dfForValues)
      }

      if (length(cellsToSelectThisLabelingTime)==length(m_u)){
        dfForValuesThisLoop <- data.frame('gene' = parameterDF[,'gene'],
                                          'cc_time' = associatedtValuesThisLoop,
                                          'labeling_time' = rep(fixedLabelingTimeThisLoop, length(cellsToSelectThisLabelingTime)),
                                          'alpha' = alphat(associatedtValuesThisLoop),
                                          'beta' = rep(NA, length(associatedtValuesThisLoop)),
                                          'gamma' = gammat(associatedtValuesThisLoop),
                                          #'p_u' = rep(NA, length(associatedtValuesThisLoop)),
                                          #'p_total' = rep(NA, length(associatedtValuesThisLoop)),
                                          'm_u' = m_u,
                                          'm_total' = totalMatureLastPeriod)
      }else{
        dfForValuesThisLoop <- data.frame('gene' = parameterDF[,'gene'],
                                          'cc_time' = associatedtValuesThisLoop,
                                          'labeling_time' = rep(fixedLabelingTimeThisLoop, length(cellsToSelectThisLabelingTime)),
                                          'alpha' = alphat(associatedtValuesThisLoop),
                                          'beta' = rep(NA, length(associatedtValuesThisLoop)),
                                          'gamma' = gammat(associatedtValuesThisLoop),
                                          #'p_u' = rep(NA, length(associatedtValuesThisLoop)),
                                          #'p_total' = rep(NA, length(associatedtValuesThisLoop)),
                                          'm_u' = m_u[cellsToSelectThisLabelingTime],
                                          'm_total' = totalMatureLastPeriod[cellsToSelectThisLabelingTime])
      }
      dfForValues <- rbind(dfForValues, dfForValuesThisLoop)
    }
    dfForValues <- dfForValues[order(dfForValues[,'cc_time']),]

    return(dfForValues)
  }
