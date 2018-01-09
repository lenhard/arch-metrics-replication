###############################################################
### Initial file parsing and data frame setup #################
###############################################################

# link to csv files
lucene = read.csv(file=file.path("arch-metrics-replication/data", "lucene-data.csv"), sep=";")
lucene <- subset(lucene, select=-c(X))

#################################################################################
####### Lucene Data Analysis ####################################################
#################################################################################

# select the metrics for which fisher tests and correlations are computed
metricNames <- names(lucene)
remove <- c("Class", "Source", "Target", "Total", "Avg.Depth", "Avg.Complexity", "Max.Depth")
metricNames <- metricNames[! metricNames %in% remove]

# compute correlations and store them in the lists
sourceCor <- list()
targetCor <- list()
totalCor <- list()

for(i in metricNames) {
 sourceCor[i] <- cor(lucene$Source, lucene[[i]], use="pairwise.complete.obs", method="spearman")
 targetCor[i] <- cor(lucene$Target, lucene[[i]], use="pairwise.complete.obs", method="spearman")
 totalCor[i] <- cor(lucene$Total, lucene[[i]], use="pairwise.complete.obs", method="spearman")
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

luceneCopy <- lucene
luceneCopy$HasSources <- ( luceneCopy$Source > 0)
luceneCopy[c("HasSources")][is.na(luceneCopy[c("HasSources")])] <- FALSE
luceneCopy$HasSources[luceneCopy$HasSources == TRUE] <- "Sources"
luceneCopy$HasSources[luceneCopy$HasSources == FALSE] <- "NoSources"

luceneCopy$HasTargets <- ( luceneCopy$Target > 0)
luceneCopy[c("HasTargets")][is.na(luceneCopy[c("HasTargets")])] <- FALSE
luceneCopy$HasTargets[luceneCopy$HasTargets == TRUE] <- "Targets"
luceneCopy$HasTargets[luceneCopy$HasTargets == FALSE] <- "NoTargets"

luceneCopy$HasTotals <- ( luceneCopy$Total > 0)
luceneCopy[c("HasTotals")][is.na(luceneCopy[c("HasTotals")])] <- FALSE
luceneCopy$HasTotals[luceneCopy$HasTotals == TRUE] <- "Totals"
luceneCopy$HasTotals[luceneCopy$HasTotals == FALSE] <- "NoTotals"


for(i in metricNames) {
    luceneCopy[i] <- (luceneCopy[[i]] > median(luceneCopy[[i]], na.rm = TRUE))
    
    if(fisher.test(table(luceneCopy[[i]], luceneCopy$HasSources))$p.value  < 0.05) {
        sourceFish[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasSources))$p.value
        sourceFishEffect[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasSources))$estimate
    }
    
    if(fisher.test(table(luceneCopy[[i]], luceneCopy$HasTargets))$p.value  < 0.05) {
        targetFish[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTargets))$p.value
        targetFishEffect[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTargets))$estimate
    }
    
        if(fisher.test(table(luceneCopy[[i]], luceneCopy$HasTotals))$p.value  < 0.05) {
        totalFish[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTotals))$p.value
        totalFishEffect[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTotals))$estimate
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

luceneCopy <- lucene
luceneCopy$HasSources <- ( luceneCopy$Source > 0)
luceneCopy[c("HasSources")][is.na(luceneCopy[c("HasSources")])] <- FALSE
luceneCopy$HasSources[luceneCopy$HasSources == TRUE] <- "Sources"
luceneCopy$HasSources[luceneCopy$HasSources == FALSE] <- "NoSources"

luceneCopy$HasTargets <- ( luceneCopy$Target > 0)
luceneCopy[c("HasTargets")][is.na(luceneCopy[c("HasTargets")])] <- FALSE
luceneCopy$HasTargets[luceneCopy$HasTargets == TRUE] <- "Targets"
luceneCopy$HasTargets[luceneCopy$HasTargets == FALSE] <- "NoTargets"

luceneCopy$HasTotals <- ( luceneCopy$Total > 0)
luceneCopy[c("HasTotals")][is.na(luceneCopy[c("HasTotals")])] <- FALSE
luceneCopy$HasTotals[luceneCopy$HasTotals == TRUE] <- "Totals"
luceneCopy$HasTotals[luceneCopy$HasTotals == FALSE] <- "NoTotals"


for(i in metricNames) {
    statistic <- quantile(luceneCopy[[i]])[4]
    luceneCopy[i] <- (luceneCopy[[i]] > statistic)
    
    if(fisher.test(table(luceneCopy[[i]], luceneCopy$HasSources))$p.value  < 0.05) {
        sourceFish[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasSources))$p.value
        sourceFishEffect[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasSources))$estimate
    }
    
    if(fisher.test(table(luceneCopy[[i]], luceneCopy$HasTargets))$p.value  < 0.05) {
        targetFish[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTargets))$p.value
        targetFishEffect[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTargets))$estimate
    }
    
        if(fisher.test(table(luceneCopy[[i]], luceneCopy$HasTotals))$p.value  < 0.05) {
        totalFish[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTotals))$p.value
        totalFishEffect[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTotals))$estimate
    }
}

fishes <- list(sourceFish, sourceFishEffect, targetFish, targetFishEffect, totalFish, totalFishEffect)

# Fisher analysis for the nineth quantile
sourceFish <- list()
sourceFishEffect <- list()
targetFish <- list()
targetFishEffect <- list()
totalFish <- list()
totalFishEffect <- list()

luceneCopy <- lucene
luceneCopy$HasSources <- ( luceneCopy$Source > 0)
luceneCopy[c("HasSources")][is.na(luceneCopy[c("HasSources")])] <- FALSE
luceneCopy$HasSources[luceneCopy$HasSources == TRUE] <- "Sources"
luceneCopy$HasSources[luceneCopy$HasSources == FALSE] <- "NoSources"

luceneCopy$HasTargets <- ( luceneCopy$Target > 0)
luceneCopy[c("HasTargets")][is.na(luceneCopy[c("HasTargets")])] <- FALSE
luceneCopy$HasTargets[luceneCopy$HasTargets == TRUE] <- "Targets"
luceneCopy$HasTargets[luceneCopy$HasTargets == FALSE] <- "NoTargets"

luceneCopy$HasTotals <- ( luceneCopy$Total > 0)
luceneCopy[c("HasTotals")][is.na(luceneCopy[c("HasTotals")])] <- FALSE
luceneCopy$HasTotals[luceneCopy$HasTotals == TRUE] <- "Totals"
luceneCopy$HasTotals[luceneCopy$HasTotals == FALSE] <- "NoTotals"


for(i in metricNames) {
    statistic <- quantile(luceneCopy[[i]], c(.9))[1]
    luceneCopy[i] <- (luceneCopy[[i]] > statistic)
    
    if(fisher.test(table(luceneCopy[[i]], luceneCopy$HasSources))$p.value  < 0.05) {
        sourceFish[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasSources))$p.value
        sourceFishEffect[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasSources))$estimate
    }
    
    if(fisher.test(table(luceneCopy[[i]], luceneCopy$HasTargets))$p.value  < 0.05) {
        targetFish[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTargets))$p.value
        targetFishEffect[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTargets))$estimate
    }
    
        if(fisher.test(table(luceneCopy[[i]], luceneCopy$HasTotals))$p.value  < 0.05) {
        totalFish[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTotals))$p.value
        totalFishEffect[i] <- fisher.test(table(luceneCopy[[i]], luceneCopy$HasTotals))$estimate
    }
}

fishes <- list(sourceFish, sourceFishEffect, targetFish, targetFishEffect, totalFish, totalFishEffect)



# correlations from significant non-size metrics to size (in terms of ncloc)
significantNonSize <- c("WMC","pmdAll", "pmdDesign", "code_smells", "complexity", "LCOM", "Ca", "RFC")
correlationsNcloc <-  list()
correlationsNclocP <-  list()
correlationsStatements <- list()
correlationsStatementsP <- list()
for(i in significantNonSize) {
  correlationsNcloc[i] <- cor(lucene$ncloc, lucene[[i]], use="pairwise.complete.obs", method="spearman")
  correlationsNclocP[i] <- cor.test(lucene$ncloc, lucene[[i]], use="pairwise.complete.obs", method="spearman", exact=FALSE)$p.value
  correlationsStatements[i] <- cor(lucene$Statements, lucene[[i]], use="pairwise.complete.obs", method="spearman")
  correlationsStatementsP[i] <- cor.test(lucene$Statements, lucene[[i]], use="pairwise.complete.obs", method="spearman", exact=FALSE)$p.value
}
