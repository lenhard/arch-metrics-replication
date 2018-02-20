#################################################################################
####### Ant Data Analysis #######################################################
#################################################################################

### Initial file parsing and data frame setup ###

# link to csv files
ant = read.csv(file=file.path("arch-metrics-replication/data", "ant-data.csv"), sep=",")
ant <- subset(ant, select=-c(X))

# select the metrics for which fisher tests and correlations are computed; predictor variables and metrics that could not be normalized (textual data) are excluded
metricNames <- names(ant)
remove <- c("Class", "Source", "Target", "Total", "Max.Depth")
metricNames <- metricNames[! metricNames %in% remove]

### compute correlations and store them in the lists
sourceCor <- list()
targetCor <- list()
totalCor <- list()

for(i in metricNames) {
 sourceCor[i] <- cor(ant$Source, ant[[i]], use="pairwise.complete.obs", method="spearman")
 targetCor[i] <- cor(ant$Target, ant[[i]], use="pairwise.complete.obs", method="spearman")
 totalCor[i] <- cor(ant$Total, ant[[i]], use="pairwise.complete.obs", method="spearman")
}

sourceCor <- sourceCor[(sourceCor > 0.3) | (sourceCor < -0.3)]
targetCor <- targetCor[(targetCor > 0.3) | (targetCor < -0.3)]
totalCor <- totalCor[(totalCor > 0.3) | (totalCor < -0.3)]

# to view the results, inspect the contents of "correlations"
# this is a list of three lists that include all correlations meeting the thresholds in the paper, sorted by the inconsistency type
correlations <- list(sourceCor, targetCor, totalCor)

### compute fisher tests and store p-values and estimates in lists

# first calculation round for cutoff at 50% threshold
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

# the results are stored in a list of six lists named "fishes". To view the results, inspect the contents of these lists
# the lists are paired:
# 	[1]: p-values for fisher tests regarding classes that are the source of inconsistencies
# 	[2]: odds ratios for fisher tests regarding classes that are the source of inconsistencies
# 	[3]: p-values for fisher tests regarding classes that are the targets of inconsistencies
# 	[4]: odds ratios for fisher tests regarding classes that are the targets of inconsistencies
# 	[5]: p-values for fisher tests regarding classes that are either targets or sources of inconsistencies
# 	[6]: odds ratios for fisher tests regarding classes that are either targets or sources of inconsistencies
# the lists only include metrics where p-values were below the significance level
# odds ratios can be interpreted in the following: The chance of having no violations in classes with the lower metric value are ESTIMATE times the chance of having no violation in classes with a higher metric value
fishes <- list(sourceFish, sourceFishEffect, targetFish, targetFishEffect, totalFish, totalFishEffect)

### Fisher analysis for the fourth quartile
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

# results viewing and interpretation works in the same way as above
fishes <- list(sourceFish, sourceFishEffect, targetFish, targetFishEffect, totalFish, totalFishEffect)


### Fisher analysis for the ninth quantile
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


### correlations from significant non-size metrics to size (in terms of ncloc and statements)
significantNonSize <- c("WMC","pmdAll", "pmdDesign", "code_smells", "complexity", "LCOM", "Ca", "RFC")
correlationsNcloc <-  list()
correlationsNclocP <-  list()
correlationsStatements <- list()
correlationsStatementsP <- list()
for(i in significantNonSize) {
  correlationsNcloc[i] <- cor(ant$ncloc, ant[[i]], use="pairwise.complete.obs", method="spearman")
  correlationsNclocP <-  cor.test(ant$ncloc, ant[[i]], use="pairwise.complete.obs", method="spearman", exact=FALSE)$p.value
  correlationsStatements[i] <- cor(ant$Statements, ant[[i]], use="pairwise.complete.obs", method="spearman")
  correlationsStatementsP[i] <- cor.test(ant$Statements, ant[[i]], use="pairwise.complete.obs", method="spearman", exact=FALSE)$p.value
}
# results can be viewed by inspecting four lists:
#	correlationsNcloc: the correlation coefficients between the metrics and the number of lines of code
#	correlationsNclocP: p-values for the correlation coefficients between the metrics and the number of lines of code
#	correlationsStatements: the correlation coefficients between the metrics and the number of statements
#	correlationsStatementsP: p-values for the correlation coefficients between the metrics and the number of statements
