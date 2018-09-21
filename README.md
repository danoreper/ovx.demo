# ovx.demo

## Installation:
Assumes R>=3.2. Assumes Jags has been installed. As of this writing all other dependencies are R packages, installed from biomaRt and CRAN. Known R dependencies include: coda, data.table, ggplot2, mvtnorm, parallel, pracma, stringr,

## Running the code:
All main scripts assume the working directory is ovx.demo/src

To locally run estimation, and generate bayesian p-values, strain-by-treatment contrast plot, etc
ovx.demo/src$ R CMD BATCH --no-save --no-restore ./ovx/main.R

To run on the cluster, use bsub; e.g.,
ovx2016/src$ bsub -M 32 R CMD BATCH --no-save --no-restore ./ovx/main.R '--args ../config/defaultBayes.yaml' 


## YAML propery files:
The default properties have been set assuming a single small machine, and also assuming the user expects the code to run in less than an hour. As such, the default sampling iterations,  number of random imputations, etc, are all low, and so the results may be somewhat innacurate. If in possession of a more powerful machine, or running on a compute cluster administered by LSF, property override files may be appropriate.

Existing property override files include:
ovx.demo/config/defaultBayes.yaml:   Property settings assuming a single 48 processor machine.
ovx.demo/config/defaultCluster.yaml: Property settings assuming the code is run on an LSF cluster (e.g. killdevil at UNC)

The default property file, which includes description of all properties, is:
ovx.demo/config/defaultParams.yaml


## Output files:
ovx.demo/output/ovx/strain.by.treatment: figures describing strain-by-treatment effects, not including the mean treatment effect,in the whole population.  


## Source files:  

### OVX project specific functionality  
ovx.demo/src/ovx
ovx.demo/src/ovx/main.R      -- driver for ovx analyis. Source this.
ovx.demo/src/ovx/ovxGibbs.R  -- functions specific to ovx analysis 
ovx.demo/src/ovx/model.delta.bug     -- jags model for estimation of strain-by-treatment effects  


### Parallelizing code, on single and multiple machines:  
ovx.demo/src/parallel/bsub.R  
ovx.demo/src/parallel/bsubScript.R  

### Loading yaml files, properties, command line args:  
ovx.demo/src/utils/loadParams.R  
ovx.demo/src/utils/loadParamsFunc.R  

### Miscellaneous generally useful methods:  
ovx.demo/src/utils/utils.R  

### Running Jags and processing gibbs samples  
ovx.demo/src/fitjags.R  


## data files:

### Forced swim data:
ovxgold/data/ovx/FINAL_FST_forR_v2.csv

### Open field data:
ovxgold/data/ovx/FINAL_OF_forR_v3.csv
