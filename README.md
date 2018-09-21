# ovx2016

## Installation:
Assumes R>=3.2. Assumes Jags has been installed. As of this writing all other dependencies are R packages, installed from biomaRt and CRAN.

## Running the code:
All main scripts assume the working directory is ovx2016/src

To run permutation testing:
ovxfold/src$ R CMD BATCH --no-save --no-restore ./ovx/perm.testing.main.R

To run permutation testing with a yaml property override file:  
```ovx2016/src$ R CMD BATCH --no-save --no-restore '--args OVERRIDEFILENAME.yaml'./ovx/perm.testing.main.R ```
e.g.,
ovx2016/src$ R CMD BATCH --no-save --no-restore '--args ../config/defaultBayes.yaml' ./ovx/perm.testing.main.R 

To run on the cluster, use bsub; e.g.,
ovx2016/src$ bsub -M 32 R CMD BATCH --no-save --no-restore ./ovx/perm.testing.main.R '--args ../config/defaultBayes.yaml' 


To run estimation, and generate bayesian p-values, caterpillar plots, etc:
ovx2016/src$ R CMD BATCH --no-save --no-restore ./ovx/estimation.main.R
Properties can be similarly overriden as for permutation testing, and script can be run on the cluster


## YAML propery files:
The default properties have been set assuming a single small machine, and also assuming the user expects the code to run in less than an hour. As such, the default sampling iterations, number of permutations, number of random imputations, etc, are all low, and so the results may be somewhat innacurate. If in possession of a more powerful machine, or running on a compute cluster administered by LSF, property override files may be appropriate.

Existing property override files include:
ovx2016/config/defaultBayes.yaml:   Property settings assuming a single 48 processor machine.
ovx2016/config/defaultCluster.yaml: Property settings assuming the code is run on an LSF cluster (e.g. killdevil at UNC)

The default property file, which includes description of all properties, is:
ovx2016/config/defaultParams.yaml


## Output files:
ovx2016/output/ovx/strain.by.treatment: figures describing strain-by-treatment effects, not including the mean treatment effect,in the whole population.  


## Source files:  

### OVX project specific functionality  
ovx2016/src/ovx
ovx2016/src/ovx/main.R      -- driver for ovx analyis. Source this.
ovx2016/src/ovx/ovxGibbs.R  -- functions specific to ovx analysis 
ovx2016/src/ovx/model.delta.bug     -- jags model for estimation of strain-by-treatment effects  


### Parallelizing code, on single and multiple machines:  
ovx2016/src/parallel/bsub.R  
ovx2016/src/parallel/bsubScript.R  

### Loading yaml files, properties, command line args:  
ovx2016/src/utils/loadParams.R  
ovx2016/src/utils/loadParamsFunc.R  

### Miscellaneous generally useful methods:  
ovx2016/src/utils/utils.R  

### Running Jags and processing gibbs samples  
ovx2016/src/fitjags.R  


## data files:

### Forced swim data:
ovxgold/data/ovx/FINAL_FST_forR_v2.csv

### Open field data:
ovxgold/data/ovx/FINAL_OF_forR_v3.csv
