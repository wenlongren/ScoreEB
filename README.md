# ScoreEB
**An Efficient Score Test Integrated with Empirical Bayes for Genome-Wide Association Studies**

**Installation**\
There two ways to install the ScoreEB package, one is installed from **CRAN**(The Comprehensive R Archive Network), and another way is installed from **GitHub**.\
**1.** install.packages("ScoreEB")\
**2.** install.packages("remotes")\
library(remotes)\
remotes::install_github("wenlongren/ScoreEB")\
{In some cases, use remotes::install_github("wenlongren/ScoreEB@main")}

**Running Example**\
library(ScoreEB)\
dir_input <- "your file path"\
genofile <- your genotype file\
phenofile <- your phenotype file\
ScoreEB(genofile, phenofile, popfile = NULL, trait.num = 1, EMB.tau = 0, EMB.omega = 0, B.Moment = 20, tol.pcg = 1e-4, iter.pcg = 100, bin = 100, lod.cutoff = 3.0, seed.num = 10000, dir_out) 

**Input File Format**\
Please refer to the mrMLM v4.0.2 (https://cran.r-project.org/web/packages/mrMLM/index.html). ScoreEB uses the input file format same with mrMLM v4.0.2. 

**Explanation of Input Parameters**\
**1.** genofile and phenofile are the **required** input file, while popfile is the **optional** input file.\
**2.** trait.num stands for computing trait from the 1st to the "trait.num".\
**3.** EMB.tau and EMB.omega are two values of hyperparameters in empirical Bayes step, which are set to 0 by default.\
**4.** B.Moment is a parameter to obtain trace of *N*x*N* matrix approximately using method of moment. B.Moment is set to 20 by default.\
**5.** tol.pcg and iter.pcg are tolerance and maximum iteration number in preconditioned conjugate gradient algorithm.\
**6.** bin is to choose the maximum score within a certain range.\
**7.** lod.cutoff is the threshold to determine identified QTNs.\
**8.** seed.num is to set the seed number.\
**9.** dir_out is the file path to save the results.

**Explanation of Output Results**\
**1.** The results file "ScoreEB.Result.csv" has 8 columns, including "Trait", "Id", "Chr", "Pos", "Score", "Beta", "Lod" and "Pvalue". \
**Note:** "Pvalue" is corresponding to the "Lod" obtained by R function pchisq(Lodx4.605,1,lower.tail=FALSE).\
**2.** The time file "ScoreEB.time.csv" includes 3 rows, which are "User", "System", "Elapse" time, respectively.
