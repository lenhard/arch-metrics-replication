###############################################################
### Initial file parsing and data frame setup #################
###############################################################

# link to csv files
antFile <- "../data/ant-data.csv"
ant = read.csv(antFile, sep=",")
ant <- subset(ant, select=-c(X))

#################################################################################
####### Ant Data Analysis #######################################################
#################################################################################

# select the metrics for which fisher tests and correlations are computed
metricNames <- names(ant)
remove <- c("Class", "Source", "Target", "Total", "Avg.Depth", "Avg.Complexity", "Max.Depth")
metricNames <- metricNames[! metricNames %in% remove]

# compute correlations and store them in the lists
sourceCor <- list()
targetCor <- list()
totalCor <- list()

for(i in metricNames) {
 sourceCor[i] <- cor(ant$Source, ant[[i]], use="pairwise.complete.obs")
 targetCor[i] <- cor(ant$Target, ant[[i]], use="pairwise.complete.obs")
 totalCor[i] <- cor(ant$Total, ant[[i]], use="pairwise.complete.obs")
}

sourceCor <- sourceCor[(sourceCor > 0.3) | (sourceCor < -0.3)]
targetCor <- targetCor[(targetCor > 0.3) | (targetCor < -0.3)]
totalCor <- totalCor[(totalCor > 0.3) | (totalCor < -0.3)]

correlations <- list(sourceCor, targetCor, totalCor)

#compute fisher tests and store p-values and estimates in lists

sourceFish <- list()
sourceFishEffect <- list()
targetFish <- list()
targetFishEffect <- list()
totalFish <- list()
totalFishEffect <- list()

antCopy <- ant
antCopy$HasSources <- ( antCopy$Source > 0)
antCopy[c("HasSources")][is.na(antCopy[c("HasSources")])] <- FALSE
antCopy$HasSources[antCopy$HasSources == TRUE] <- "Sources"
antCopy$HasSources[antCopy$HasSources == FALSE] <- "NoSources"

antCopy$HasTargets <- ( antCopy$Target > 0)
antCopy[c("HasTargets")][is.na(antCopy[c("HasTargets")])] <- FALSE
antCopy$HasTargets[antCopy$HasTargets == TRUE] <- "Targets"
antCopy$HasTargets[antCopy$HasTargets == FALSE] <- "NoTargets"

antCopy$HasTotals <- ( antCopy$Total > 0)
antCopy[c("HasTotals")][is.na(antCopy[c("HasTotals")])] <- FALSE
antCopy$HasTotals[antCopy$HasTotals == TRUE] <- "Totals"
antCopy$HasTotals[antCopy$HasTotals == FALSE] <- "NoTotals"


for(i in metricNames) {
    antCopy[i] <- (antCopy[[i]] > median(antCopy[[i]], na.rm = TRUE))
    
    if(fisher.test(table(antCopy[[i]], antCopy$HasSources))$p.value  < 0.05) {
        sourceFish[i] <- fisher.test(table(antCopy[[i]], antCopy$HasSources))$p.value
        sourceFishEffect[i] <- fisher.test(table(antCopy[[i]], antCopy$HasSources))$estimate
    }
    
    if(fisher.test(table(antCopy[[i]], antCopy$HasTargets))$p.value  < 0.05) {
        targetFish[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTargets))$p.value
        targetFishEffect[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTargets))$estimate
    }
    
        if(fisher.test(table(antCopy[[i]], antCopy$HasTotals))$p.value  < 0.05) {
        totalFish[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTotals))$p.value
        totalFishEffect[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTotals))$estimate
    }
}

fishes <- list(sourceFish, sourceFishEffect, targetFish, targetFishEffect, totalFish, totalFishEffect)

# The chance of having no violations in classes with the lower metric value are ESTIMATE times the chance of having no violation in classes with a higher metric value


# Fisher analysis for the fourth quartile
sourceFish <- list()
sourceFishEffect <- list()
targetFish <- list()
targetFishEffect <- list()
totalFish <- list()
totalFishEffect <- list()

antCopy <- ant
antCopy$HasSources <- ( antCopy$Source > 0)
antCopy[c("HasSources")][is.na(antCopy[c("HasSources")])] <- FALSE
antCopy$HasSources[antCopy$HasSources == TRUE] <- "Sources"
antCopy$HasSources[antCopy$HasSources == FALSE] <- "NoSources"

antCopy$HasTargets <- ( antCopy$Target > 0)
antCopy[c("HasTargets")][is.na(antCopy[c("HasTargets")])] <- FALSE
antCopy$HasTargets[antCopy$HasTargets == TRUE] <- "Targets"
antCopy$HasTargets[antCopy$HasTargets == FALSE] <- "NoTargets"

antCopy$HasTotals <- ( antCopy$Total > 0)
antCopy[c("HasTotals")][is.na(antCopy[c("HasTotals")])] <- FALSE
antCopy$HasTotals[antCopy$HasTotals == TRUE] <- "Totals"
antCopy$HasTotals[antCopy$HasTotals == FALSE] <- "NoTotals"


for(i in metricNames) {
    statistic <- quantile(antCopy[[i]])[4]
    antCopy[i] <- (antCopy[[i]] > statistic)
    
    if(fisher.test(table(antCopy[[i]], antCopy$HasSources))$p.value  < 0.05) {
        sourceFish[i] <- fisher.test(table(antCopy[[i]], antCopy$HasSources))$p.value
        sourceFishEffect[i] <- fisher.test(table(antCopy[[i]], antCopy$HasSources))$estimate
    }
    
    if(fisher.test(table(antCopy[[i]], antCopy$HasTargets))$p.value  < 0.05) {
        targetFish[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTargets))$p.value
        targetFishEffect[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTargets))$estimate
    }
    
        if(fisher.test(table(antCopy[[i]], antCopy$HasTotals))$p.value  < 0.05) {
        totalFish[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTotals))$p.value
        totalFishEffect[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTotals))$estimate
    }
}

fishes <- list(sourceFish, sourceFishEffect, targetFish, targetFishEffect, totalFish, totalFishEffect)


# Fisher analysis for the ninth quantile
sourceFish <- list()
sourceFishEffect <- list()
targetFish <- list()
targetFishEffect <- list()
totalFish <- list()
totalFishEffect <- list()

antCopy <- ant
antCopy$HasSources <- ( antCopy$Source > 0)
antCopy[c("HasSources")][is.na(antCopy[c("HasSources")])] <- FALSE
antCopy$HasSources[antCopy$HasSources == TRUE] <- "Sources"
antCopy$HasSources[antCopy$HasSources == FALSE] <- "NoSources"

antCopy$HasTargets <- ( antCopy$Target > 0)
antCopy[c("HasTargets")][is.na(antCopy[c("HasTargets")])] <- FALSE
antCopy$HasTargets[antCopy$HasTargets == TRUE] <- "Targets"
antCopy$HasTargets[antCopy$HasTargets == FALSE] <- "NoTargets"

antCopy$HasTotals <- ( antCopy$Total > 0)
antCopy[c("HasTotals")][is.na(antCopy[c("HasTotals")])] <- FALSE
antCopy$HasTotals[antCopy$HasTotals == TRUE] <- "Totals"
antCopy$HasTotals[antCopy$HasTotals == FALSE] <- "NoTotals"


for(i in metricNames) {
    statistic <- quantile(antCopy[[i]], c(.9))[1]
    antCopy[i] <- (antCopy[[i]] > statistic)
    
    if(fisher.test(table(antCopy[[i]], antCopy$HasSources))$p.value  < 0.05) {
        sourceFish[i] <- fisher.test(table(antCopy[[i]], antCopy$HasSources))$p.value
        sourceFishEffect[i] <- fisher.test(table(antCopy[[i]], antCopy$HasSources))$estimate
    }
    
    if(fisher.test(table(antCopy[[i]], antCopy$HasTargets))$p.value  < 0.05) {
        targetFish[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTargets))$p.value
        targetFishEffect[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTargets))$estimate
    }
    
        if(fisher.test(table(antCopy[[i]], antCopy$HasTotals))$p.value  < 0.05) {
        totalFish[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTotals))$p.value
        totalFishEffect[i] <- fisher.test(table(antCopy[[i]], antCopy$HasTotals))$estimate
    }
}

fishes <- list(sourceFish, sourceFishEffect, targetFish, targetFishEffect, totalFish, totalFishEffect)



