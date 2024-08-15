## ----include=FALSE, echo=FALSE------------------------------------------------
# --------------------------------------------------------------
# generate R file with code from this file
# --------------------------------------------------------------
knitr::purl(input = "trials_interim.qmd", output = "trials_interim.r")


## ----include=TRUE, echo=TRUE--------------------------------------------------
# --------------------------------------------------------------
# packages
# --------------------------------------------------------------
library(rpact)
packageVersion("rpact")

# computations below are done with rpact Version >= 4.0.0


## ----include=TRUE, echo=TRUE--------------------------------------------------
alpha <- 0.05
beta <- 0.2

# median OS
m1 <- 6 * 12
m2 <- 8 * 12
hr <- m1 / m2

# accrual and dropout
accrualTime <- c(0, 12)
accrualIntensity <- 100
dropoutRate1 <- 0.025
dropoutRate2 <- 0.025
dropoutTime <- 12
maxNumberOfSubjects <- accrualIntensity * max(accrualTime)


## ----include=TRUE, echo=TRUE--------------------------------------------------
nevent <- getSampleSizeSurvival(lambda1 = getLambdaByMedian(m2), 
                                lambda2 = getLambdaByMedian(m1), 
                                sided = 2, alpha = alpha, beta = beta)
mdd_no_interim <- nevent$criticalValuesEffectScaleLower
nevent <- ceiling(nevent$maxNumberOfEvents)

nevent
mdd_no_interim


## ----include=TRUE, echo=TRUE--------------------------------------------------
info <- c(1 / 3, 2 / 3, 1)
designEff <- getDesignGroupSequential(informationRates = info,
                                      typeOfDesign = "asOF", sided = 1, 
                                      alpha = alpha / 2)

samplesizeEFF <- getSampleSizeSurvival(design = designEff,
                                       lambda2 = log(2) / m1, hazardRatio = hr,
                                       dropoutRate1 = dropoutRate1, 
                                       dropoutRate2 = dropoutRate2,
                                       dropoutTime = dropoutTime,
                                       accrualTime = accrualTime, 
                                       accrualIntensity = accrualIntensity)

# number of events
nevents <- ceiling(as.vector(samplesizeEFF$eventsPerStage))

# add exploratory updated analysis
nevents <- c(nevents, 500)
nevents


## ----include=TRUE, echo=TRUE--------------------------------------------------
# Local significance levels (two-sided)
alphas <- designEff$stageLevels * 2
alphas

# MDD 
hrMDD <- as.vector(samplesizeEFF$criticalValuesEffectScale)
hrMDD

# analysis timepoints
time <- as.vector(round(samplesizeEFF$analysisTime, 1))

# predict analysis time of updated analysis
grid <- 76.3965
probEvent <- getEventProbabilities(time = grid, lambda1 = getLambdaByMedian(m2), 
                                   lambda2 = getLambdaByMedian(m1),
                                   dropoutRate1 = dropoutRate1, 
                                   dropoutRate2 = dropoutRate2, 
                                   dropoutTime = dropoutTime,
                                   accrualTime = accrualTime, 
                                   accrualIntensity = accrualIntensity, 
                                   maxNumberOfSubjects = maxNumberOfSubjects)
expEvent <- probEvent$overallEventProbabilities * maxNumberOfSubjects  
cbind(grid, expEvent)

time <- c(time, grid)

# power loss through the interims
# Futility HR and information at interim
futilityHR <- c(1, 0.9)
eventsInterim <- nevents[1:2]

# Calculate Z-score corresponding to futilityHR 
# Note: Minus sign is added to be consistent with directionUpper = FALSE 
# in getPowerSurvival call below.
allocationRatioPlanned <- 1

# approximate info according to Schoenfeld's formula
infoInterim <- eventsInterim * allocationRatioPlanned/(1 + allocationRatioPlanned) ^ 2 
ZscoreInterim <- -log(futilityHR) / sqrt(1 / infoInterim)

designFutility <- getDesignGroupSequential(sided = 1, alpha = alpha / 2,
                                           informationRates = info,
                                           typeOfDesign = "noEarlyEfficacy",
                                           futilityBounds = ZscoreInterim,
                                           bindingFutility = FALSE)

# Calculate power and timeline under H1 with original sample size if design includes 
# the non-binding futility (and the futility boundary is adhered to)
power <- getPowerSurvival(designFutility,
                          allocationRatioPlanned = allocationRatioPlanned,
                          lambda2 = log(2) / m1, hazardRatio = hr,
                          dropoutRate1 = dropoutRate1, dropoutRate2 = dropoutRate2, 
                          dropoutTime = dropoutTime,
                          accrualTime = accrualTime, accrualIntensity = accrualIntensity,
                          maxNumberOfEvents = nevents[3], directionUpper = FALSE)
summary(power)


## ----include=TRUE, echo=TRUE--------------------------------------------------
suppressMessages(library(lubridate))
fpi <- ymd("2020-04-23")
d <- round(30 * (time - floor(time)))
dates <- c(fpi, fpi %m+% months(floor(time)) + days(d))
dates


## ----include=TRUE, echo=TRUE--------------------------------------------------
# --------------------------------------------------------------
# functions
# --------------------------------------------------------------

# Quantile function for a Weibull with a potential cure proportion.
qWeibullCure <- function(p, p0, shape = 1, scale){
  res <- rep(NA, length(p))
  ind1 <- (p <= p0)
  ind2 <- (p > p0)
  res[ind1] <- Inf
  res[ind2] <- qweibull(1 - (p[ind2] - p0) / (1 - p0), shape = shape, scale = scale)
  return(res)  
}

# simulate Weibull distributed time-to-event data in one group, 
# where administrative censoring is applied after a pre-specified
# number of events. A cure proportion and censoring due to drop-out 
# can be specified in addition.
rWeibull1arm <- function(shape = 1, scale, cure = 0, recruit, dropout = 0, 
                         start.accrual = 0, cutoff, seed = NA){
  
  # shape             Weibull shape parameter. 
  # scale:            Weibull scale parameter.
  # cure:             Proportion of patients assumed to be cured, i.e. with an 
  #                   event at +infty.
  # recruit:          Recruitment.
  # dropout:          Drop-out rate, on same time scale as med.
  # start.accrual:    Time unit where accrual should start. Might be useful when 
  #                   simulating multi-stage trials.
  # cutoff:           Cutoff, #events the final censored data should have 
  #                   (can be a vector of multiple cutoffs).
  # seed:             If different from NA, seed used to generate random numbers.
  #
  # Kaspar Rufibach, June 2014
  
  if (is.na(seed) == FALSE){set.seed(seed)}
  
  n <- sum(recruit)
  
  # generate arrival times
  arrive <- rep(1:length(recruit), times = recruit)
  arrivetime <- NULL
  for (i in 1:n){arrivetime[i] <- runif(1, min = arrive[i] - 1, max = arrive[i])}
  arrivetime <- start.accrual + sort(arrivetime)
  
  # generate event times: Exp(lambda) = Weibull(shape = 1, scale = 1 / lambda)
  eventtime <- qWeibullCure(runif(n), p0 = cure, shape = shape, scale = scale)
  
  # Apply drop-out. Do this before applying the cutoff below, 
  # in order to correctly count necessary #events.
  dropouttime <- rep(Inf, n)
  if (dropout > 0){dropouttime <- rexp(n, rate = dropout)}
  event.dropout <- ifelse(eventtime > dropouttime, 0, 1)
  time.dropout <- ifelse(event.dropout == 1, eventtime, dropouttime)   
  
  # observed times, taking into account staggered entry
  tottime <- arrivetime + eventtime
  
  # find cutoff based on number of targeted events
  # only look among patients that are not considered dropped-out
  time <- data.frame(matrix(NA, ncol = length(cutoff), nrow = n))
  event <- time
  cutoff.time <- rep(NA, length(cutoff))
  
  for (j in 1:length(cutoff)){
    cutoff.time[j] <- sort(tottime[event.dropout == 1])[cutoff[j]]
    
    # apply administrative censoring at cutoff
    event[event.dropout == 1, j] <- ifelse(
      tottime[event.dropout == 1] > cutoff.time[j], 0, 1)
    event[event.dropout == 0, j] <- 0
    
    # define time to event, taking into account both types of censoring
    time[event.dropout == 1, j] <- ifelse(event[, j] == 1, 
                                          eventtime, cutoff.time[j] - 
                                            arrivetime)[event.dropout == 1]    
    time[event.dropout == 0, j] <- pmin(cutoff.time[j] - arrivetime, 
                                        time.dropout)[event.dropout == 0]
    
    # remove times for patients arriving after the cutoff
    rem <- (arrivetime > cutoff.time[j])
    if (TRUE %in% rem){time[rem, j] <- NA}
  }
  
  # generate output
  tab <- data.frame(cbind(1:n, arrivetime, eventtime, tottime, dropouttime, time, event))
  colnames(tab) <- c("pat", "arrivetime", "eventtime", "tottime", "dropouttime", 
                     paste("time cutoff = ", cutoff, sep = ""), 
                     paste("event cutoff = ", cutoff, sep = ""))
  
  res <- list("cutoff.time" = cutoff.time, "tab" = tab)
  return(res)
}

# simulate Weibull distributed time-to-event data in two arms, 
# where administrative censoring is applied after a pre-specified 
# number of events. A cure proportion and censoring due to drop-out 
# can be specified in addition.
rWeibull2arm <- function(shape = c(1, 1), scale, cure = c(0, 0), recruit, 
                         dropout = c(0, 0), start.accrual = c(0, 0), cutoff, seed = NA){
  
  # shape             2-d vector of Weibull shape parameter. 
  # scale             2-d vector of Weibull scale parameter.
  # cure:             2-d vector with cure proportion assumed in each arm.
  # recruit:          List with two elements, vector of 
  #                   recruitment in each arm.
  # dropout:          2-d vector with drop-out rate for each arm, 
  #                   on same time scale as med.
  # start.accrual:    2-d vector of time when accrual should start. 
  #                   Might be useful when simulating multi-stage trials.
  # cutoff:           Cutoff, #events the final censored data should have 
  #                   (can be a vector of multiple cutoffs).
  # seed:             If different from NA, seed used to generate random numbers.
  #
  # Kaspar Rufibach, June 2014
  
  if (is.na(seed) == FALSE){set.seed(seed)}
  
  dat1 <- rWeibull1arm(scale = scale[1], shape = shape[1], recruit = recruit[[1]], 
                       cutoff = 1, dropout = dropout[1], cure = cure[1], 
                       start.accrual = start.accrual[1], seed = NA)$tab
  dat2 <- rWeibull1arm(scale = scale[2], shape = shape[2], recruit = recruit[[2]], 
                       cutoff = 1, dropout = dropout[2], cure = cure[2], 
                       start.accrual = start.accrual[2], seed = NA)$tab
  
  n <- c(nrow(dat1), nrow(dat2))
  
  # treatment variable
  tmt <- factor(c(rep(0, n[1]), rep(1, n[2])), levels = 0:1, labels = c("A", "B"))
  
  arrivetime <- c(dat1[, "arrivetime"], dat2[, "arrivetime"])
  eventtime <- c(dat1[, "eventtime"], dat2[, "eventtime"])
  tottime <- c(dat1[, "tottime"], dat2[, "tottime"])
  dropouttime <- c(dat1[, "dropouttime"], dat2[, "dropouttime"])
  
  # Apply drop-out. Do this before applying the cutoff below, 
  #in order to correctly count necessary #events.
  event.dropout <- ifelse(eventtime > dropouttime, 0, 1)
  time.dropout <- ifelse(event.dropout == 1, eventtime, dropouttime)   
  
  # find cutoff based on number of targeted events
  # only look among patients that are not considered dropped-out
  time <- data.frame(matrix(NA, ncol = length(cutoff), nrow = sum(n)))
  event <- time
  cutoff.time <- rep(NA, length(cutoff))
  
  for (j in 1:length(cutoff)){
    cutoff.time[j] <- sort(tottime[event.dropout == 1])[cutoff[j]]
    
    # apply administrative censoring at cutoff
    event[event.dropout == 1, j] <- ifelse(tottime[event.dropout == 1] > 
                                             cutoff.time[j], 0, 1)
    event[event.dropout == 0, j] <- 0
    
    # define time to event, taking into account both types of censoring
    time[event.dropout == 1, j] <- ifelse(event[, j] == 1, eventtime, 
                                    cutoff.time[j] - arrivetime)[event.dropout == 1]    
    time[event.dropout == 0, j] <- pmin(cutoff.time[j] - arrivetime, 
                                    time.dropout)[event.dropout == 0]
    
    # remove times for patients arriving after the cutoff
    rem <- (arrivetime > cutoff.time[j])
    if (TRUE %in% rem){time[rem, j] <- NA}
  }
  
  # generate output
  tab <- data.frame(cbind(1:sum(n), tmt, arrivetime, eventtime, tottime, 
                          dropouttime, time, event))
  colnames(tab) <- c("pat", "tmt", "arrivetime", "eventtime", "tottime", 
                     "dropouttime", paste("time cutoff = ", cutoff, sep = ""), 
                     paste("event cutoff = ", cutoff, sep = ""))

  res <- list("cutoff.time" = cutoff.time, "tab" = tab)
  return(res)
}


## ----include=TRUE, echo=TRUE--------------------------------------------------
library(survival)

# accrual
recruit1 <- rep(accrualIntensity / 2, max(accrualTime))
recruit2 <- recruit1
recruit <- list(recruit1, recruit2)
n <- sum(unlist(recruit))
start.accrual <- c(0, 0)

# drop-out (2.5% annually)
dropout.year <- 0.025
dropout.month <- -log(1 - dropout.year) / 12
dropout <- rep(dropout.month, 2)

# cure proportion 
cure <- c(0, 0)

# --------------------------------------------------
# simulate a trial stopped for futility
# --------------------------------------------------
med <- c(m1, m1)
lambda <- log(2) / med
scale <- 1 / lambda
j <- 1
cut1 <- nevents[j] + 3

trial1 <- rWeibull2arm(scale = scale, shape = c(1, 1), 
                        cure = cure, recruit = recruit, dropout = dropout, 
                        start.accrual = start.accrual, cutoff = cut1, 
                        seed = 8)
tmt <- trial1$tab["tmt"][, 1]
time <- trial1$tab[paste("time cutoff = ", cut1, sep = "")][, 1]
event <- trial1$tab[paste("event cutoff = ", cut1, sep = "")][, 1]

cox1 <- summary(coxph(Surv(time, event) ~ tmt))
HR1 <- exp(coef(cox1)[1])
p1 <- coef(cox1)[1, "Pr(>|z|)"]

HR1
p1
cut1

# date of CCOD
trial1$cutoff.time
d <- round(30 * (trial1$cutoff.time - floor(trial1$cutoff.time)))
c(fpi, fpi %m+% months(floor(trial1$cutoff.time)) + days(d))


# --------------------------------------------------
# simulate a trial stopped for efficacy
# --------------------------------------------------
med <- m1 / c(1, hr)
lambda <- log(2) / med
scale <- 1 / lambda
j <- 2
cut2 <- nevents[j] - 2

trial2 <- rWeibull2arm(scale = scale, shape = c(1, 1), 
                        cure = cure, recruit = recruit, dropout = dropout, 
                        start.accrual = start.accrual, cutoff = cut2, 
                        seed = 7)
tmt <- trial2$tab["tmt"][, 1]
time <- trial2$tab[paste("time cutoff = ", cut2, sep = "")][, 1]
event <- trial2$tab[paste("event cutoff = ", cut2, sep = "")][, 1]

cox2 <- summary(coxph(Surv(time, event) ~ tmt))
HR2 <- exp(coef(cox2)[1])
ci2 <- exp(confint(coxph(Surv(time, event) ~ tmt)))
p2 <- coef(cox2)[1, "Pr(>|z|)"]

hrMDD[j]
HR2
ci2
p2
cut2
trial2$cutoff.time

# recalculate local significance level based on observed number of events
info2 <- c(info[1], cut2 / nevents[3], 1)
designUpdate2 <- getDesignGroupSequential(sided = 1, alpha = alpha / 2, 
                                          beta = beta,
                                          informationRates = info2, 
                                          typeOfDesign = "asOF")

samplesize2 <- getSampleSizeSurvival(design = designUpdate2,
                                     lambda2 = log(2) / m1, hazardRatio = hr,
                                     dropoutRate1 = dropoutRate1, 
                                     dropoutRate2 = dropoutRate2,
                                     dropoutTime = dropoutTime,
                                     accrualTime = accrualTime, 
                                     accrualIntensity = accrualIntensity)

alphas2 <- as.vector(samplesize2$criticalValuesPValueScale)
hrMDD2 <- as.vector(samplesize2$criticalValuesEffectScale)

# original and updated local significance levels
rbind(alphas, alphas2 * 2)
rbind(hrMDD, hrMDD2)

# median unbiased estimate, assuming the futility interim had happened 
# at the prespecified number of events with an observed HR of 0.9
results2 <- getDataset(
  overallEvents = c(nevents[1], cut2),
  overallLogRanks = c(log(0.9)/sqrt(4 / nevents[1]), coef(cox2)[1, "z"]),
  overallAllocationRatio = c(1, 1))
adj_result2 <- getAnalysisResults(design = designUpdate2,
                                  dataInput = results2,
                                  stage = 2, directionUpper = FALSE)
adj_result2$medianUnbiasedEstimates[2]
c(adj_result2$finalConfidenceIntervalLowerBounds[2], 
  adj_result2$finalConfidenceIntervalUpperBounds[2])

# date of CCOD
trial2$cutoff.time
d <- round(30 * (trial2$cutoff.time - floor(trial2$cutoff.time)))
c(fpi, fpi %m+% months(floor(trial2$cutoff.time)) + days(d))

