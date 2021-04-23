# ScoreEB

**An Efficient Score Test Integrated with Empirical Bayes for Genome-Wide Association Studies**

Perform association test within linear mixed model framework using score test integrated with Empirical Bayes for genome-wide association study. Firstly, score test was conducted for each single nucleotide polymorphism (SNP) under linear mixed model framework, taking into account the genetic relatedness and population structure. And then all the potentially associated SNPs were selected with a less stringent criterion. Finally, all the selected SNPs were performed Empirical Bayes in a multi-locus model to identify the true quantitative trait nucleotide (QTN).

**Installation**\
install.packages("devtools")\
install.packages("remotes")\
remotes::install_github("wenlongren/ScoreEB")

**Running Example**\
library(ScoreEB)\
setwd("C:/Users/ThinkPad/Desktop/Lucky/ScoreEB/GitHub/")\
genfile <- "Sim_200K.gen"\
famfile <- "Sim_200K.fam"\
bimfile <- "Sim_200K.bim"\
phenofile <- "SimY1.csv"\
IniPthreshold <- 0.001\
LOD <- 3.0\
repetition <- 3\
result <- ScoreEBtest(genfile,famfile,bimfile,phenofile,IniPthreshold,LOD,repetition)\
Note: setwd("path"), path is where you place the input data. And the example data can be found in "wenlongren/ScoreEB/Data" folder.

**Explanation of Input Parameters and Output Results**
**1.	Input four files: genfile, famfile, bimfile and phenofile**
(1)	*genfile*: a  M x N dimension matrix, M is the number of markers, N is the number of individuals, there is no row names and column names. And the file format is **.gen**.
(2)	*famfile*: a N x 6 dimension matrix, and it has column names, which are “famid”, “id”, “	father”,  “mother”, “sex”, and “pheno”, respectively. And the file format is **.fam**. Note: we do not use phenotype here.
(3)	*bimfile*: a M x 6 dimension matrix, and it has column names, which are “chr”, “id”, “dist”, “pos”, “A1” and “A2”, respectively. And the file format is **.bim**. 
(4)	*phenofile*: a N x R dimension matrix, R is the number of phenotypes. The row names are individual names, and the column names are the names of phenotypes. And the file format is **.csv**.

**2.	Input three parameters: IniPthreshold, LOD, repetition**
(1)	*IniPthreshold*: the p threshold used in the first stage (single locus score test stage), which is much less than bonferroni correction p threshold. IniPthreshold is usually set as 0.01 or 0.001, and it can be modified according to different data.
(2)	*LOD*: logarithm of odds, which is the threshold used in the second stage (multiple locus empirical bayes stage), that is for the final result. LOD is usually set as 2.5 or 3.0, and it can be modified by users.
(3)	*repetition*: for real data analysis, it stands for the number of phenotypes need to be computed; for simulation data analysis, it stands for number of repetitions need to be computed.

**3.	Output Results**
The output result is a S x 11 dimension matrix, and it has column names, which are “ntime”, “chr”, “pos”, “id”, “A1”, “A2”, “freq”, “score”, “p”, “lod” and “lrtb”, respectively. Note: “ntime” stands for which phenotype or which repetition that be detected, “lrtb” stands for the effect of marker obtained by likelihood ratio test.


**Please refer:\
WL Ren, et al.(2020) An Efficient Score Test Integrated with Empirical Bayes for Genome-Wide Association Studies.**
