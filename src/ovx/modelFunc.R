ovx = new.env(hash=T )
fp = file.path

source("./fitjags.R")
source("./parallel/bsub.R")


## Load all datasets, parameters, and compute matching information up front, and store it in a list.
ovx$loadAllInputs <- function()
{
    datasets = ovx$loadDatasets()

    ##build the matching pairs of animals in the same strain and batch, but having opposite treatment.
    ##This is done upfront to avoid having a different set of matching pairs for different phenotypes in the same dataset
    matchings = list()

    ##Although matchings are built upfront, they do need to be done separately per dataset, as a different set of
    ##mice were assayed in the OF test vs the Forced swim test
    for(d in 1:length(datasets))
    {
        matchings[[d]] = ovx$getDeltaMatchings(datasets[[d]]$frame, prop$ovx$numImp, prop$ovx$idcol)
    }

    ##A mapping that is needed for the guts of runJags code, describing which gibbs variables are observed, and what their name
    ## is in the raw dataframe
    gibbsObserveInfo = data.frame(variable = c("ranef.strain"),
                                  factorForIndexConversion = c("STRAIN.1") ,
                                  keep = T, stringsAsFactors = F)

    ##also needed for guts of runJags code, to know which variables are factors. Admittedly, it could probably impute this,
    ##but didn't have time to go into guts of code
    colsToIndex      = "STRAIN.1"
    
    ##file containing the jags model, which is the same for every phenotype
    jagsModel        = "./ovx/model.delta.bug"

    ## ##In this case, we will be contrasting ranef.strain. In other analyses (not shown), we've contrasted other variables.
    ## effectToContrast = "ranef.strain"
    
    return(list(datasets  = datasets,
                matchings = matchings,
                gibbsObserveInfo = gibbsObserveInfo,
                jagsModel = jagsModel,
                colsToIndex = colsToIndex))
}


## Load up data from two different datasets into a list,
## in which each element is itself a list consisting of a dataframe, and the names of the phenoytpes of
## interest in that dataframe
ovx$loadDatasets <- function()
{
    datasets = list()
    datasets[[1]] = list(frame = fread(fp(prop$data, "ovx/FINAL_OF_forR_v3.csv")),
                         phens = c("TotDist", "PctCenter", "VMovNo", "MovTime", "MrgDist", "CtrDist"))
    
    datasets[[2]] = list(frame = fread(fp(prop$data, "ovx/FINAL_FST_forR_v2.csv")),
                         phens = "PctImb")
    return(datasets)
}


## Generate a list of matching information of length numImputations.
## Each imputation corresponds to a different way to match all animals.
##
## Each element of the list of imputations is a dataframe, and each row of the dataframe corresponds to a pair of animal ids,
## that are in the same batch, and have the same strain, but have opposite treatments; SHAM surgery or OVX.
## As there is more than one way to match animals, we repeat for numImputations
##
## dataset: the original, unmatched dataset, specifying phenotypes for every animal individually
## numImputations: the number of times we want to generate a matching
## idcol: the column corresponding to animal id.
ovx$getDeltaMatchings <- function(dataset, numImputations, idcol)
{
    matchon  = c("STRAIN", "BATCH")
    matchoff = c("TREATMENT")
    dataset$TREATMENT = factor(dataset$TREATMENT, levels = c("SHAM", "OVX"))

    ## generate a collection of matched ids, repeating for as many imputations as desired
    matchings = list()
    for(r in 1:numImputations)
    {
        matchings[[r]]     = ovx$getDeltaMatching(dataset,
                                                  idcol = idcol,
                                                  matchon = matchon, 
                                                  matchoff = matchoff)
    }
    return(matchings)
}

## Generate a dataframe in which animal ids are matched, such that only animals with the same matchon condition can be matched,
## and they can only be matched if they have opposite matchoff conditions.
##
## df: the data frame with treated and untreated outcomes
## matchon: fields which have to be identical in order to match a pair
## matchoff: a field which has to be different in order to match a pair; i.e., the treatment factor
## idcol: the column in df corresponding to a sample ID
##
## within a given matchon condition, if there are more animals having one treatment factor than another,
## the excess animals in one of the treatments are dropped, with the specific dropped animals randomly chosen.
## Some next imputation will randomly end up including these animals in the matching, and a different set of animals will be dropped.
##
## Returns dataframe in which each row has ID.1 and ID.2, corresponding to a pair of matched animals, as well as the conditions
## within which that matched pair is matched.
ovx$getDeltaMatching <- function(df, matchon, matchoff, idcol)
{
    
    df = data.table(df)
    
    ##this 'by' statement groups ids into groups that have the same matchon condition (e.g., batch and strain)
    dfnew = df[,j=list(ID = paste(get(idcol), collapse=",")), by = c(matchon)]
    
    l1.ID = c()
    l2.ID = c()
    offlevels = levels(df[[matchoff]])
    for(i in 1:nrow(dfnew))
    {
        row = dfnew[i,]
        IDs = unlist(strsplit(row$ID, ","))
        matchingOffValues = df[[matchoff]][match(IDs, df[[idcol]])]
        
        level1.IDs  = IDs[which(matchingOffValues == offlevels[1])]
        level2.IDs  = IDs[which(matchingOffValues == offlevels[2])]
        
        numpairs = min(length(level1.IDs), length(level2.IDs))
        if(numpairs == 0)
        {
            next
        }

        ##randomly grab numpairs worth of ids whose treatement is on, and equal number whose treament is off
        ## the ith element of one vector is matched to the ith element of the other vector
        level1.IDs = sample(x = level1.IDs, replace = F,size = numpairs)
        level2.IDs = sample(x = level2.IDs, replace = F,size = numpairs)

        ## accumulate the matchings for this matchon condition to those for the other matchon conditions
        l1.ID = c(l1.ID, level1.IDs)
        l2.ID = c(l2.ID, level2.IDs)
    }

    ## create a data frame out of the matched ids, and fill in the matchon condition for each pair of matched ids
    matched = data.frame(ID.1 = l1.ID, ID.2 = l2.ID)
    for(cname in matchon)
    {
        for(j in 1:2)
        {
            matched[[paste0(cname,".",j)]] = df[[cname]][match(l1.ID, df[[idcol]])]
        }
    }

    return(matched)
}



## Run the gibbs sampling for our model of strain by treament effects
##
##
## datasets: a list of datasets, in which each dataset is itself a list of dataframe, and the phenotypes in that dataframe
##
## matchings: a list whose lenght is as long as there are datasets.
## Each element of matchings, corresponding to a dataset, is itself a list of data frames,
## of length equal to the number of desired matching imputations.
## For each imputation, the dataframe contains a set of matched pairs of mice, in which each pair is in the same strain and batch,
## and one animal is sham, while the other had ovariectomy.
##
## observeinfo: the variables which should be recorded from the Gibbs sampling process. Gibbs sampling runs faster and uses less memory if fewer variable are observed
##
## colsToIndex: the columns that are factors, which the gibbs sampler needs to be aware of.
##
## Returns a list of mcmc objects (gibbs samples), with each element corresponding to a different phenotype. All phenotypes
## in any dataset are included.

ovx$runModels <- function(datasets,
                          matchings,
                          observeInfo,
                          jagsModel,
                          colsToIndex)
{
    pracma::tic()
    
    gibbs.samples   = list()
    for(d in 1:length(datasets))
    {
        dataset.phenz = datasets[[d]]
        dataset       = dataset.phenz$frame
        phenz         = dataset.phenz$phens
        for(phen in phenz)
        {
            gibbs.samples[[phen]] = ovx$runModel(dataset = dataset,
                                                 phen = phen,
                                                 matchings = matchings[[d]], ##the matchings for 2 datasets will differ, as they assay different mice
                                                 modelname = jagsModel,
                                                 colsToIndex = colsToIndex,
                                                 observeInfo = observeInfo)
        }
    }
    pracma::toc()
    return(gibbs.samples)
}

##
## A helper function to run the gibbs models for a single dataset, and for a single phenotype, with the appropriate matchings
## for that dataset. See ovx$vrunModels for more detail on params.
##
##
ovx$runModel <- function(dataset,
                         phen,
                         matchings,
                         modelname,
                         colsToIndex,
                         observeInfo,
                         idcol = "ANIMAL.ID")
{
    gibbs.samples = list()
    y.samples     = list()
    commandList = list()

    ##initialize the accumulator for invocations of runjags, and queue up the calls untill runAll is called
    accum = ovx$getAccumulator()
    accum$init(fitjags$runJags)    

    y = dataset[[phen]]
    ids = dataset[[idcol]]

    ##for each imputed matching, add a call to run the model
    for (imp in 1:length(matchings))
    {
        print(paste0("imputation: ", imp))
        matching = matchings[[imp]]
        ## make a copy
        trainingData = matching

        trainingData$y1 = y[match(trainingData$ID.1, ids)]
        trainingData$y2 = y[match(trainingData$ID.2, ids)]
        trainingData$y  = trainingData$y2 - trainingData$y1
        trainingData[[phen]] = trainingData$y

        accum$addCall(list(trainingData= trainingData,
                           modelname   = modelname,
                           observeInfo = observeInfo,
                           iter        = prop$ovx$iter,
                           burninFrac  = prop$ovx$burninfrac,
                           nchains     = 1,
                           thin        = prop$ovx$thin,
                           phen        = phen,
                           colsToIndex = colsToIndex))
    }
    ##blocking call, distributes all model fitting to various nodes on the cluster (if on the cluster), otherwise tries to use multicore
    outfiles = accum$runAll()

    ##collate the results of every imputation, reading them from file if necessary, and make a giant mcmc object including
    ## all imputations
    if(accum$outputs.files)
    {
        collated = ovx$collateGibbs(outfiles)
    } else {
        collated = outfiles
    }
    collated = mcmc(do.call(rbind, collated))
    
    return(collated)
 }       


## Collect outputs from the cluser of running gibbs sampling, with a separate job per matching imputation.
##
## walk over all outfiles from running on the cluster. Each outfile corresponds to an imputation.
## load them into memory ('clusterOut' is in every file)
## and put them into a list. Skip the imputation if there was an error.
ovx$collateGibbs <- function(outfiles)
{
    upperLimit = length(outfiles)

    gibbs.samples.all = list()
    counter = 1
    for(r in 1:upperLimit)
    {
        outfile = outfiles[r]
        x = try(load(file = outfile)) ##clusterOut is loaded into memory 
        if(class(x)=="try-error")
        {
            print(r)
        } else {
            ##clusterOut is loaded into memory 
            gibbs.samples.all[[r]] = clusterOut
        }
    }
    return((gibbs.samples = gibbs.samples.all))
}

## Get the 'accumulator' environment, used for parallelizing function invocations.
## It contains functions which queue up jobs, submit them to the cluster (or multicore locally), and blocking until the submitted jobs finish. If properties indicate that this run is not on the cluster, creates a multicore running accumulator
##
## In a bit more detail, the accumulator, whether on the cluster or not, allows the user to call 1) 'init', with the function that will be used, and any global arguments, 2) 'addCall', which involves passing in a list of arguments to be invoked with the function, and 'runAll', which will invoke all the calls which have been queued up, and block untill it is finished.
##
##
##
ovx$getAccumulator <- function()
{
    if(prop$onCluster)
    {
        batchSize = 1
        bsubCommand = bsub$get.default.killdevil.bsub(numProcessPerNode = 1, memoryLimit.GB = prop$ovx$clusterMemLim,  queue=prop$ovx$clusterQueue)
        accum = bsub$get.bsub.accumulator("./ovx/modelFunc.R", bsubCommand, batchSize=batchSize)
    } else {

        accum = bsub$get.mc.accumulator(mc.cores= prop$ovx$mc.cores)
    }
    return(accum)
}





#################################
##
##Functions for parsing the gibbs samples as fit by the model of strain-by-treatment effects,
## And reporting the results, specifically through flag plots.
###########################################################################3




## For every phenotype, and for every pair of levels of effectName, (e.g., NOD and B6 are levels of strain, an effect),
## extracts two sets of samples of effect size, subtracts them from one another.
## Using this sample of effect differences, computes
## a p-value like quantity for that contrast.
##
## gibbs.samples.all: a list whose elements correspond to the gibbs samples for fitting the model to a single phenotype
## effectName: the effect that is to be compared; e.g., strain or strain-by-treatment
## shortname: whether or not to shorten the names output by gibbs sampling to ones that are easier to display, without integer indexes.
##
## Returns a dataframe with the phenotype, the two effect levels,
## the p-value of the difference between those levels, a discretized significance corresponding to the p-value,
## and the direction of the effect difference-- +1 or -1
##
ovx$getComparisonFrame <- function(gibbs.samples.all, effectName, shortname = F)
{
    ##get the column names from gibbs samples which correspond to levels of effectName
    ovx$getColsForEffect <- function(gibbs.samples, effectName)
    {
        regexEffectName = gsub(effectName, pattern = "\\.", replacement = "\\\\.")
        colstring = paste0(regexEffectName,"\\[.*\\]")
        colsForEffect = grepl(pattern = colstring, perl=T, x = colnames(gibbs.samples))
        colsForEffect = colnames(gibbs.samples)[colsForEffect]
        return(colsForEffect)
    }

    ##For a given effect type and levelName (as produced by jags output), comes up with a shortened name easier to display 
    ovx$getShortName <- function(effectName, levelName)
    {
        regexEffectName = gsub(effectName, pattern = "\\.", replacement = "\\\\.")
        out = gsub(levelName, pattern = paste0(regexEffectName,"|\\[|\\]"), replacement = "")
        return(out)
    }

    ##compute the minimum tail probability (MTP) of the gibbs sample vector Y;
    ##whatever the smaller empirical probability is, p(y<0) or p(y>=0).
    ## return the MTP and the direction of the difference in effects
    ovx$get.p.value <- function(y)
    {
        pval = sum(y<0)/length(y)
        if(pval<1-pval)
        {
            direction = 1
        } else {
            pval = min(1-pval, pval)
            direction = -1
        }
        return(data.frame(pval = pval, direction = direction))
    }

    ##looping over every phenotype, parse the gibbs samples, store in dfs
    dfs = list()
    for(phen in names(gibbs.samples.all))
    {
        gibbs.samples = gibbs.samples.all[[phen]]
        compareCols = ovx$getColsForEffect(gibbs.samples, effectName)

        c1 = paste0(effectName, ".1")
        c2 = paste0(effectName, ".2")

        ## compare every pair of effect levels
        for(i in 1:length(compareCols))
        {
            col.i = compareCols[i]
            for(j in 1:length(compareCols))
            {
                col.j = compareCols[j]
                
                colvec = unname(gibbs.samples[,col.j] - gibbs.samples[,col.i])
                if(i==j)
                {
                    df = data.frame(pval=1, direction=0) ## this may be obvious, but its helpful when plotting
                } else {
                    df = ovx$get.p.value(colvec)
                }
                df[[c1]] = ifelse(shortname, ovx$getShortName(effectName, col.i), col.i)
                df[[c2]] = ifelse(shortname, ovx$getShortName(effectName, col.j), col.j)

                df$phen = phen
                dfs = util$appendToList(dfs, df)
            }
        }
    }
    dfs = do.call(rbind, dfs)

    dfs$discreteValue = cut(dfs$pval, prop$ovx$discreteBuckets, include.lowest=T)
    return(dfs)
}


## Creates a flag plot of constrast significances for a given effect, and writes it to file.
##
## dfplot.compare: the dataframe output of getComparisonFrame above
## effectName: the effect of interest
## fname: filename to write plot to.
## textsize: axis text size. can be null for default text size
## themargin: margin size between panels
## reallybig: boolean indicating whther 1 or more rows of panels are needed.

ovx$getPlot.compare.big <- function(dfplot.compare,
                                    effectName,
                                    fname,
                                    textsize,
                                    themargin,
                                    reallybig)
{
    
    dfplot.compare$phen = factor(x = dfplot.compare$phen, ordered = T, levels=c("TotDist", "PctCenter", "VMovNo", "MovTime", "MrgDist", "CtrDist","PctImb"))

    limz = c(-3,3)
    dfplot.compare$pval[-log10(dfplot.compare$pval)>(limz[2]-.01)] = 10^(-(max(limz)-.01))
    dfplot.compare$dir.string = ""
    dfplot.compare$dir.discrete = factor(paste0(dfplot.compare$direction, dfplot.compare$discreteValue))
    dfplot.compare$p.with.direction = -log10(dfplot.compare$pval) * dfplot.compare$direction
    dfplot.compare$stringg = 32
    dfplot.compare$stringg[dfplot.compare$pval<.05] = 8

    effect1 = paste0(effectName, ".1")
    effect2 = paste0(effectName, ".2")

    aplot = ggplot(dfplot.compare, aes(x=get(effect1), y=get(effect2), fill = p.with.direction, label=stringg, shape = stringg))

    aplot = aplot + geom_tile()
    aplot = aplot + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
    aplot = aplot + scale_shape_identity()
    aplot = aplot + xlab("Strain1" ) + ylab("Strain2")
    aplot = aplot + scale_fill_gradient2(low = "red", high = "blue",
                                         name = expression('-Log'[10]*'(p)'%.%'EffectDirection'),
                                         limits =limz)
    aplot = aplot + geom_point(size = 2.5, color="white")


    if(reallybig)
    {
        aplot = aplot + facet_wrap(~phen, nrow = 2)
    } else {
        aplot = aplot + facet_wrap(~phen)
    }
    aplot = aplot + theme(aspect.ratio=1, panel.margin = unit(themargin, "lines"))

    aplot = aplot + xlab("Strain 1")
    aplot = aplot + ylab("Strain 2")
    ##browser()
    if(!is.null(textsize))
    {
        aplot = aplot + theme(axis.text=element_text(size = textsize))
    }
    print(fname)
    pdf(fname, width=15, height=9)
    print(aplot)
    dev.off()
    
}
