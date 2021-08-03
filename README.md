# ScoreEB
**An Efficient Score Test Integrated with Empirical Bayes for Genome-Wide Association Studies**

**Installation**\
install.packages("remotes")\
library(remotes)\
remotes::install_github("wenlongren/ScoreEB")

**Running Example**\
library(ScoreEB)\
dir_input <- "your file path"\
genofile <- paste0(dir_input, "mrMLM.Genotype.csv")\
phenofile <- paste0(dir_input, "mrMLM.SimPheno1.csv")\
ScoreEB(genofile, phenofile, popfile = NULL, trait.num = 1, B.Moment = 20, tol.pcg = 1e-6, iter.pcg = 500, bin = 100, lod.cutoff = 3.0, dir_out) 

**Input File Format**\
Please refer to the mrMLM v4.0.2 (https://cran.r-project.org/web/packages/mrMLM/index.html). ScoreEB2 uses the input file format same with mrMLM v4.0.2. 

**Explanation of Input Parameters**\
**1.** genofile and phenofile are the **required** input file, while popfile is the **optional** input file.\
**2.** trait.num stands for computing trait from the 1st to the "trait.num".\
**3.** B.Moment is a parameter to obtain trace of *N*x*N* matrix approximately using method of moment. B.Moment is set to 20 by default.\
**4.** tol.pcg and iter.pcg are tolerance and maximum iteration number in preconditioned conjugate gradient algorithm.\
**5.** bin is to choose the maximum score within a certain range.\
**6.** lod.cutoff is the threshold to determine identified QTNs.\
**7.** dir_out is the file path to save the results. 

**Explanation of Output Results**\
**1.** The results file "ScoreEB.Result.csv" has 8 columns, including "Trait", "Id", "Chr", "Pos", "Score", "Beta", "Lod" and "Pvalue". \
**Note:** "Pvalue" is corresponding to the "Lod" obtained by R function pchisq(Lodx4.605,1,lower.tail=FALSE).

**2.** The time file "ScoreEB.time.csv" includes 3 rows, which are "User", "System", "Elapse" time, respectively.
