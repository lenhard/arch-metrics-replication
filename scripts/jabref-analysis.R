#################################################################################
####### JabRef Data Analysis #######################################################
#################################################################################

### Initial file parsing and data frame setup ###

all = read.csv(file=file.path("arch-metrics-replication/data", "jabref-data.csv"), sep=",")
all <- subset(all, select=-c(X))

# select the metrics for which fisher tests and correlations are computed; predictor variables and metrics that could not be normalized (textual data) are excluded
metricNames <- names(all)
remove <- c("Class", "Source", "Target", "Total", "Max.Depth")
metricNames <- metricNames[! metricNames %in% remove]

### compute correlations and store them in the lists
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

# results viewing and interpretation works in the same way as above
fishes <- list(sourceFish, sourceFishEffect, targetFish, targetFishEffect, totalFish, totalFishEffect)

### Fisher analysis for the 9th quantile
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

### correlations from significant non-size metrics to size (in terms of ncloc and statements)
significantNonSize <- c("WMC","pmdAll", "pmdDesign", "code_smells", "complexity", "LCOM", "Ca", "RFC")
correlationsNcloc <-  list()
correlationsNclocP <-  list()
correlationsStatements <- list()
correlationsStatementsP <- list()
for(i in significantNonSize) {
  correlationsNcloc[i] <- cor(all$ncloc, all[[i]], use="pairwise.complete.obs", method="spearman")
  correlationsNclocP[i] <- cor.test(all$ncloc, all[[i]], use="pairwise.complete.obs", method="spearman", exact=FALSE)$p.value
  correlationsStatements[i] <- cor(all$Statements, all[[i]], use="pairwise.complete.obs", method="spearman")
  correlationsStatementsP[1] <- cor.test(all$Statements, all[[i]], use="pairwise.complete.obs", method="spearman", exact=FALSE)$p.value
}
# results can be viewed by inspecting four lists:
#	correlationsNcloc: the correlation coefficients between the metrics and the number of lines of code
#	correlationsNclocP: p-values for the correlation coefficients between the metrics and the number of lines of code
#	correlationsStatements: the correlation coefficients between the metrics and the number of statements
#	correlationsStatementsP: p-values for the correlation coefficients between the metrics and the number of statements
