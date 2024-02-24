#BAMM Trait Evolution Analysis (part 2, hopefully cleaner)
#Morgan Muell
library(BAMMtools)
library(coda)
#library(svglite)


setwd("~/Documents/Masters/ch2_work/BAMM/tests/redhead/fuzzybanded/30genrun")

######
#PREP#
######
#make tree ultrametric
tree <- read.tree("20mil.treefile")
is.ultrametric(tree)
ut <- force.ultrametric(tree)
is.ultrametric(ut)
write.tree(ut, file="ultra-out2.treefile")

#fuzzing trait files
#"fuzzing up" trait files
traitfile <- read.delim(file="stripe_traits.txt", header=F)
traits <- traitfile$V2
x <- runif(length(traits), min = -.125, max = .125)
x
fuzzytraits <- traits + x
fuzzytraits
write.table(fuzzytraits, file="fuzzy_banded.txt")
#head over to Excel to finalize, unsure how to do it in R yet

#getting appropriate priors for control file
tree <- read.tree("ultra-out2.treefile")
setBAMMpriors(phy=tree, traits="fuzzy_banded.txt")
#use output to copy in priors into control file

###RUN BAMM###


##########
#ANALYSIS#
##########

#Load necessary files (tree and event data with rate parameters and 
#their topology mappings)
tree <- read.tree("ultra-out2.treefile")
edata <- getEventData(tree, type='trait', 
                      eventdata = "event_data.txt", burnin=0.1)
stripeedata <- getEventData(tree, type='trait', 
                            eventdata="../../../stripes/fuzzystripes/30genrun/event_data.txt", 
                            burnin=0.1)
spotedata <- getEventData(tree, type='trait', eventdata = "event_data.txt", burnin=0.1)
rededata <- getEventData(tree, type='trait', 
                         eventdata="../../../redhead/fuzzyredhead/30genrun/event_data.txt", 
                         burnin=0.1)
bandedata <- getEventData(tree, type='trait', 
                          eventdata="../../../banded/fuzzybanded/30genrun/event_data.txt", 
                          burnin=0.1)

#Assess convergence (do this, plus run multiple BAMM runs)
mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

#####Identifying Number of Rate Shifts#####
#View shift posterior distribution and calculate best model
shift_probs <- summary(edata) 
postfile <- "mcmc_out.txt"
bfmat <- computeBayesFactors(postfile, 
                             expectedNumberOfShifts=1, burnin=0.1)
bfmat #this returns pairwise matrix of Bayes factors


#plot whether prior or posterior sampled shift sets and how often
plotPrior(mcmcout, expectedNumberOfShifts=1)
title(main="Sampling Distribution - redhead30") #change title for each
#Compute credible set of shift configurations (account for 95% of
#observations in the data; DOES NOT identify core shifts)
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, 
                        threshold=5, set.limit = 0.95)
css$number.distinct #number of distinct shift configs in data
summary(css) ###LOOK AT THIS###
plot.credibleshiftset(css) #f is % of post samples with that shift
dev.off()

#Random Samples of Distinct Shift Configurations
#rank is a given set of core shifts; samples randomly displayed will
#always share the core shift but may vary in non-core shifts also
#present
dsc <- distinctShiftConfigurations(edata, expectedNumberOfShifts=1,
                                   threshold=5)
plot.bammshifts(dsc, edata, rank=2, legend=F)

#####Identifying Location of Rate Shifts#####
#Marginal Shift Probabilities (independent of rest of tree)
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, cex=.5) 
#in the above plot, longer branches have higher probabilities

#Marginal Odds Ratios (determines core shifts and likely locations)
branch_priors <- getBranchShiftPriors(tree, 
                                      expectedNumberOfShifts = 1)
plot(branch_priors, cex=.5)
#ratios above 1 have marginal shift probs greater than expected
#under prior alone 
mo <- marginalOddsRatioBranches(edata, expectedNumberOfShifts = 1)
values <- mo$edge.length
values[values <= 5] = NA
highprobs <- values
plot(mo, cex=.3)
edgelabels(highprobs, frame="none", cex=0.5) 
#how could I get it to only plot edge lengths above 1?

#####Identifying Weight of Rate Shifts#####
plotRateThroughTime(edata, ratetype="auto", ylim=c(0,1)) #ew
getMeanBranchLengthTree(edata, rate="trait")
rtt <- getRateThroughTimeMatrix(edata)
meanTraitRate <- colMeans(rtt$beta)
plot(meanTraitRate ~ rtt$times)

tiprates <- getTipRates(edata, statistic="mean")
tips <- tiprates$beta.avg
test <- getTipRates(edata) #every sample in the posterior
#but I want rate changes in instances where only the core shifts
#occurred for the scenarios in the 95% interval...but they vary
#I think an average would be best, since samples in the posterior
#change so much?


###-- generate overall best shift configuration
#maximum a posteriori shift config --distinct shift config with highest posterior prob
##    I think this is the one with the highest posterior prob from css summary
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
plot.bammdata(best, lwd=1.25, labels=TRUE, cex=0.57, legend=T)
addBAMMshifts(best, cex=2)
title(main="MAP -- spots")



#################
#Files To Export#
#################
#Table of tip rates
write.csv(tips, file="30genredhead-tip-rates.csv")
#Mean phylorate Plot (averaged rates over all models)
plot.bammdata(edata, lwd=2, legend=T)
title(main="Mean Phylorate - fuzzybanded30")

#highest prob in CSS - with locations
png(filename="MAP_best_banded.png", res=300, units="in", height=6, width=5)
plot.bammdata(best, lwd=1.25, labels=TRUE, cex=0.3, legend=T)
addBAMMshifts(best, cex=1)
title(main="MAP -- banded")
dev.off()

#Credible shift sets
plot.credibleshiftset(css)
dev.off()
#Marginal odds
plot(mo, cex=.5)
edgelabels(highprobs, frame="none", cex=0.5) 
title(main="Marginal Odds - Banded")

#all MAPs - every phenotype

#png(filename="MAP_best_striped.png", res=300, units="in", height=6, width=5)
svg(filename="MAP_best_striped.svg", height=6, width=5)
beststripe <- getBestShiftConfiguration(stripeedata, expectedNumberOfShifts=1, threshold=5)
plot.bammdata(beststripe, lwd=1.25, labels=TRUE, cex=0.3, legend=T)
addBAMMshifts(beststripe, cex=1)
title(main="MAP -- striped")
dev.off()

#png(filename="MAP_best_spotted.png", res=300, units="in", height=6, width=5)
svg(filename="MAP_best_spotted.svg", height=6, width=5)
bestspot <- getBestShiftConfiguration(spotedata, expectedNumberOfShifts=1, threshold=5)
plot.bammdata(bestspot, lwd=1.25, labels=TRUE, cex=0.3, legend=T)
addBAMMshifts(bestspot, cex=1)
#title(main="MAP -- spotted")
dev.off()

#png(filename="MAP_best_redhead.png", res=300, units="in", height=6, width=5)
svg(filename="MAP_best_redhead.svg", height=6, width=5)
bestred <- getBestShiftConfiguration(rededata, expectedNumberOfShifts=1, threshold=5)
plot.bammdata(bestred, lwd=1.25, labels=TRUE, cex=0.3, legend=T)
addBAMMshifts(bestred, cex=1)
title(main="MAP -- redhead")
dev.off()

#png(filename="MAP_best_banded.png", res=300, units="in", height=6, width=5)
svg(filename="MAP_best_banded.svg", height=6, width=5)
bestband <- getBestShiftConfiguration(bandedata, expectedNumberOfShifts=1, threshold=5)
plot.bammdata(bestband, lwd=1.25, labels=TRUE, cex=0.3, legend=T)
addBAMMshifts(bestband, cex=1)
title(main="MAP -- banded")
dev.off()

