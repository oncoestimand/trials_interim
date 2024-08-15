## ----include=FALSE, echo=FALSE------------------------------------------------
# --------------------------------------------------------------
# generate R file with code from this file
# --------------------------------------------------------------
knitr::purl(input = "trials_interim.qmd", output = "trials_interim.r")


## ----include=TRUE, echo=TRUE--------------------------------------------------
# --------------------------------------------------------------
# packages
# --------------------------------------------------------------
packs <- "rpact"
for (i in 1:length(packs)){library(packs[i], character.only = TRUE)}


## ----include=TRUE, echo=TRUE--------------------------------------------------
library(rpact)
packageVersion("rpact")

# computations below with rpact Version >= 4.0.0


# design specifications
alpha <- 0.05
beta <- 0.2
m1 <- 6 * 12
m2 <- 8 * 12
hr <- m1 / m2

# further specifications
accrualTime <- c(0, 12)
accrualIntensity <- 100
dropoutRate1 <- 0.025
dropoutRate2 <- 0.025
dropoutTime <- 12
maxNumberOfSubjects <- accrualIntensity * max(accrualTime)

# Required events without interim
nevent <- getSampleSizeSurvival(lambda1 = getLambdaByMedian(m2), lambda2 = getLambdaByMedian(m1), 
                                sided = 2, alpha = alpha, beta = beta)
mdd_no_interim <- nevent$criticalValuesEffectScaleLower
nevent <- ceiling(nevent$maxNumberOfEvents)

nevent
mdd_no_interim

# OBF design without interims (to get futility bound at second interim as per 1) above)
# I use OBF also for the first interim, although strictly speaking not needed.
# it does not make a numerical difference though
# one could make this "perfect" by using a hand-knitted alpha-spending function in rpact, but I did not bother to do so
# can do once we put the script on github
info <- c(1 / 3, 2 / 3, 1)
designEff <- getDesignGroupSequential(informationRates = info,
                                      typeOfDesign = "asOF", sided = 1, 
                                      alpha = alpha / 2)

samplesizeEFF <- getSampleSizeSurvival(design = designEff,
                                           lambda2 = log(2) / m1, hazardRatio = hr,
                                           dropoutRate1 = dropoutRate1, dropoutRate2 = dropoutRate2, 
                                           dropoutTime = dropoutTime,
                                           accrualTime = accrualTime, accrualIntensity = accrualIntensity)

# number of events
nevents <- ceiling(as.vector(samplesizeEFF$eventsPerStage))

# add exploratory updated analysis
nevents <- c(nevents, 500)
nevents

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
probEvent <- getEventProbabilities(time = grid, lambda1 = getLambdaByMedian(m2), lambda2 = getLambdaByMedian(m1),
                                   dropoutRate1 = dropoutRate1, dropoutRate2 = dropoutRate2, 
                                   dropoutTime = dropoutTime,
                                   accrualTime = accrualTime, accrualIntensity = accrualIntensity, 
                                   maxNumberOfSubjects = maxNumberOfSubjects)
expEvent <- probEvent$overallEventProbabilities * maxNumberOfSubjects  
cbind(grid, expEvent)

time <- c(time, grid)

# power loss through the interims
# Futility HR and information at interim
futilityHR <- c(1, 0.9)
eventsInterim <- nevents[1:2]

# Calculate Z-score corresponding to futilityHR 
# Note: Minus sign is added to be consistent with directionUpper = FALSE in getPowerSurvival call below.
allocationRatioPlanned <- 1
infoInterim <- eventsInterim * allocationRatioPlanned/(1 + allocationRatioPlanned) ^ 2 # approximate info according to Schoenfeld's formula
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
                          dropoutRate1 = dropoutRate1, dropoutRate2 = dropoutRate2, dropoutTime = dropoutTime,
                          accrualTime = accrualTime, accrualIntensity = accrualIntensity,
                          maxNumberOfEvents = nevents[3], directionUpper = FALSE)
summary(power)

# dates of analyses
library(lubridate)
fpi <- ymd("2020-04-23")
d <- round(30 * (time - floor(time)))
dates <- c(fpi, fpi %m+% months(floor(time)) + days(d))
dates








