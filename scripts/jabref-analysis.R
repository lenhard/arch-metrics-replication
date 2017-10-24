###############################################################
### Initial file parsing and data frame setup #################
###############################################################

jabRefAllFile <- "c:/workspaces/git/arch-metrics-replication/data/jabref-data.csv"
all = read.csv(jabRefAllFile, sep=";")
all <- subset(all, select=-c(X))

# select the metrics for which fisher tests and correlations are computed
metricNames <- names(all)
remove <- c("Class", "Source", "Target", "Total", "Avg.Depth", "Avg.Complexity", "Max.Depth")
metricNames <- metricNames[! metricNames %in% remove]

# compute correlations and store them in the lists
sourceCor <- list()
targetCor <- list()
totalCor <- list()

for(i in metricNames) {
 sourceCor[i] <- cor(all$Source, all[[i]], use="pairwise.complete.obs", method="spearman")
 targetCor[i] <- cor(all$Target, all[[i]], use="pairwise.complete.obs", method="spearman")
 totalCor[i] <- cor(all$Total, all[[i]], use="pairwise.complete.obs", method="spearman")
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

allCopy <- all
allCopy$HasViolations <- ( allCopy$Source > 0)
allCopy[c("HasViolations")][is.na(allCopy[c("HasViolations")])] <- FALSE
allCopy$HasViolations[allCopy$HasViolations == TRUE] <- "Violations"
allCopy$HasViolations[allCopy$HasViolations == FALSE] <- "NoViolations"

allCopy$HasTargets <- ( allCopy$Target > 0)
allCopy[c("HasTargets")][is.na(allCopy[c("HasTargets")])] <- FALSE
allCopy$HasTargets[allCopy$HasTargets == TRUE] <- "Targets"
allCopy$HasTargets[allCopy$HasTargets == FALSE] <- "NoTargets"

allCopy$HasTotals <- ( allCopy$Total > 0)
allCopy[c("HasTotals")][is.na(allCopy[c("HasTotals")])] <- FALSE
allCopy$HasTotals[allCopy$HasTotals == TRUE] <- "Totals"
allCopy$HasTotals[allCopy$HasTotals == FALSE] <- "NoTotals"


for(i in metricNames) {
    allCopy[i] <- (allCopy[[i]] > median(allCopy[[i]], na.rm = TRUE))
    
    if(fisher.test(table(allCopy[[i]], allCopy$HasViolations))$p.value  < 0.05) {
        sourceFish[i] <- fisher.test(table(allCopy[[i]], allCopy$HasViolations))$p.value
        sourceFishEffect[i] <- fisher.test(table(allCopy[[i]], allCopy$HasViolations))$estimate
    }
    
    if(fisher.test(table(allCopy[[i]], allCopy$HasTargets))$p.value  < 0.05) {
        targetFish[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTargets))$p.value
        targetFishEffect[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTargets))$estimate
    }
    
    if(fisher.test(table(allCopy[[i]], allCopy$HasTotals))$p.value  < 0.05) {
        totalFish[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTotals))$p.value
        totalFishEffect[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTotals))$estimate
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

allCopy <- all
allCopy$HasViolations <- ( allCopy$Source > 0)
allCopy[c("HasViolations")][is.na(allCopy[c("HasViolations")])] <- FALSE
allCopy$HasViolations[allCopy$HasViolations == TRUE] <- "Violations"
allCopy$HasViolations[allCopy$HasViolations == FALSE] <- "NoViolations"

allCopy$HasTargets <- ( allCopy$Target > 0)
allCopy[c("HasTargets")][is.na(allCopy[c("HasTargets")])] <- FALSE
allCopy$HasTargets[allCopy$HasTargets == TRUE] <- "Targets"
allCopy$HasTargets[allCopy$HasTargets == FALSE] <- "NoTargets"

allCopy$HasTotals <- ( allCopy$Total > 0)
allCopy[c("HasTotals")][is.na(allCopy[c("HasTotals")])] <- FALSE
allCopy$HasTotals[allCopy$HasTotals == TRUE] <- "Totals"
allCopy$HasTotals[allCopy$HasTotals == FALSE] <- "NoTotals"


for(i in metricNames) {
    # compute the fourth quartile
    statistic <- quantile(allCopy[[i]])[4]
    allCopy[i] <- (allCopy[[i]] > statistic)
    
    if(fisher.test(table(allCopy[[i]], allCopy$HasViolations))$p.value  < 0.05) {
        sourceFish[i] <- fisher.test(table(allCopy[[i]], allCopy$HasViolations))$p.value
        sourceFishEffect[i] <- fisher.test(table(allCopy[[i]], allCopy$HasViolations))$estimate
    }
    
    if(fisher.test(table(allCopy[[i]], allCopy$HasTargets))$p.value  < 0.05) {
        targetFish[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTargets))$p.value
        targetFishEffect[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTargets))$estimate
    }
    
    if(fisher.test(table(allCopy[[i]], allCopy$HasTotals))$p.value  < 0.05) {
        totalFish[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTotals))$p.value
        totalFishEffect[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTotals))$estimate
    }
    
}

fishes <- list(sourceFish, sourceFishEffect, targetFish, targetFishEffect, totalFish, totalFishEffect)


# Fisher analysis for the 9th quantile
sourceFish <- list()
sourceFishEffect <- list()
targetFish <- list()
targetFishEffect <- list()
totalFish <- list()
totalFishEffect <- list()

allCopy <- all
allCopy$HasViolations <- ( allCopy$Source > 0)
allCopy[c("HasViolations")][is.na(allCopy[c("HasViolations")])] <- FALSE
allCopy$HasViolations[allCopy$HasViolations == TRUE] <- "Violations"
allCopy$HasViolations[allCopy$HasViolations == FALSE] <- "NoViolations"

allCopy$HasTargets <- ( allCopy$Target > 0)
allCopy[c("HasTargets")][is.na(allCopy[c("HasTargets")])] <- FALSE
allCopy$HasTargets[allCopy$HasTargets == TRUE] <- "Targets"
allCopy$HasTargets[allCopy$HasTargets == FALSE] <- "NoTargets"

allCopy$HasTotals <- ( allCopy$Total > 0)
allCopy[c("HasTotals")][is.na(allCopy[c("HasTotals")])] <- FALSE
allCopy$HasTotals[allCopy$HasTotals == TRUE] <- "Totals"
allCopy$HasTotals[allCopy$HasTotals == FALSE] <- "NoTotals"


for(i in metricNames) {
    # compute the ninth quantile
    statistic <- quantile(allCopy[[i]], c(.9))[1]
    allCopy[i] <- (allCopy[[i]] > statistic)
    
    if(fisher.test(table(allCopy[[i]], allCopy$HasViolations))$p.value  < 0.05) {
        sourceFish[i] <- fisher.test(table(allCopy[[i]], allCopy$HasViolations))$p.value
        sourceFishEffect[i] <- fisher.test(table(allCopy[[i]], allCopy$HasViolations))$estimate
    }
    
    if(fisher.test(table(allCopy[[i]], allCopy$HasTargets))$p.value  < 0.05) {
        targetFish[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTargets))$p.value
        targetFishEffect[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTargets))$estimate
    }
    
    if(fisher.test(table(allCopy[[i]], allCopy$HasTotals))$p.value  < 0.05) {
        totalFish[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTotals))$p.value
        totalFishEffect[i] <- fisher.test(table(allCopy[[i]], allCopy$HasTotals))$estimate
    }
}

fishes <- list(sourceFish, sourceFishEffect, targetFish, targetFishEffect, totalFish, totalFishEffect)
