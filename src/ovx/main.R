##  Author: doreper 
## Entry point for the ovariectomy paper analysis (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5289389/)
##
## The ovariectomy experiment  was looking for strain-by-ovariectomy effects.
## Experiment was performed as follows:
## Using a panel of 37 inbred mouse strains, with inbred females tested over 46 batches over 2 years,
## for every strain and batch combination,
## one subset of female mice was treated with ovariectomy,
## while the other mice were given a sham surgery.
## Afterwards, these mice were assayed using measures of open-field behavior, ( 6 phenotypes)
## and then the  forced swim test (1 phenotype, percent immobility). Some mice were only tested in OF.
##
## The analysis in this code is attempting to model strain-by-treatment effects, and in doing so, enable:
## 1) Estimation of effects with credible intervals
## 2) Critically, some measure of the significance of individual strain1-strain2 differences in treatment effect;
## i.e., the significance of contrasts between strain-by-treatment effects for different strains.
##
## Given the paucity of our data, as well as the reasonableness of the assumption that strain effects are drawn from a normal distribution,
## we wanted to model strain-by-treatment effects as random effects--a conservative way of modeling which
## should give more accurate results. However, although computing the significance of fixed effect contrasts is trivial,
## and packages even exist to compute fixed effect contrasts in the context of mixed models,
## there were no available packages to compute the significance of contrasts between levels of random variables in a frequentist model.
## As such, we instead employed a Bayesian hierarchical approach to modeling our data, which readily lent itself to
## computing credible intervals for contrasts between strain-by-treatment levels, even while modelling their effects as drwn from a normal.
##
## Essentially, using a JAGS sampler, the code is fitting the bayesian equivalent (with uninformative priors) of:
## deltaY ~ 1 + (1|strain)
##
## and then computing contrast signficances between the levels of strain.
##
## In the model above, each deltaY is the difference in phenotype between a paired treated and untreated mouse of the same strain and same batch.
## As there are multiple possible pairings of mice, the analysis belowe uses multiple imputation and effectively averages over
## the results for every randomly selected set of pairs.
##
## Accordingly, given that we are modeling deltaY, not Y, '1' corresponds to the treatment effect and '1|strain' corresponds to the
## strain-by-treatment-effect
##
###############################################################################

library(ggplot2)
library(data.table)
library(utils)

## load other needed R files
sources = function()
{
    ##Loads properties in from YAML file (as specified at the command line invoking R) into the global 'prop' variable"
    source("./utils/loadParams.R")

    source("./ovx/modelFunc.R")
    source("./utils/utils.R")
    source("./parallel/bsub.R")

}

sources()

rebuildData = prop$ovx$rebuildData

pracma::tic()
if(rebuildData)
{ 
    print("Load up input data, parameters, and randomly matched pairs, to compute strain-by-treatment effects")
    inputs  = ovx$loadAllInputs()

    print("running gibbs sampling to fit model of strain-by-treatment")
    gibbs.samples = ovx$runModels(datasets        = inputs$datasets,
                                matchings       = inputs$matchings,
                                observeInfo     = inputs$gibbsObserveInfo,
                                jagsModel       = inputs$jagsModel,
                                colsToIndex     = inputs$colsToIndex)

    print("done building output data")
    ##takes a while, so save
    save(file =fp(prop$output, "ovx/ovx.RData"), list = ls())
    pracma::toc()
    
} else {
    print("loading saved output data, to be used for reporting")
    load(file =fp(prop$output, "ovx/ovx.RData"))
    sources()
    print("done loading")
}


pracma::tic()


## parse the results of gibbs sampling, and in particular, taking the differences between strain-by-treatment effects for different strains,
## then computing a sort of p-value for each such difference based on the minimum tail probability
comparisonFrame = ovx$getComparisonFrame(gibbs.samples, "ranef.strain", shortname = T)

## STRAIN.BY.TREATMENT contrast p-values, for two sepecific phenotypes that were interesting
comparisonFrame.small = comparisonFrame[comparisonFrame$phen %in% c("PctCenter", "CtrDist"),]

## Generate STRAIN.BY.TREATMENT contrast plots
froot = fp(prop$output, "ovx", "strain.by.treatment")
dir.create(froot, showWarnings = F, recursive = T)

## Plots of all effects in a big facet plot. Written to output folder.
ovx$getPlot.compare.big(comparisonFrame,
                        effect = "ranef.strain",
                        fname = fp(froot, "compare.big.pdf"),
                        textsize = 4,
                        themargin = .5,
                        reallybig = T)

## plots of just the interesting effects. Written to output folder
ovx$getPlot.compare.big(comparisonFrame.small,
                        effect = "ranef.strain",
                        fname = fp(froot, "compare.small.pdf"),
                        textsize = NULL,
                        themargin = 2,
                        reallybig = F)
pracma::toc()
