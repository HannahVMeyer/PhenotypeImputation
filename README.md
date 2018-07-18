# PhenotypeImputation

In biological research, individual samples are increasingly phenotyped for multiple correlated and/or uncorrelated traits.
Unevitably, in experiments with large samples sizes, some samples will not have been successfully phenotyped for each trait
(technical difficulites, measurement not passing quality control etc.). Instead of exlcuding the samples with missing 
phenotype information from downstream analysis, the missing measurments can be imputed based on information of the samples' other
measurements, measurements within the cohort and the correlation between samples and/or traits within the experiment.

Here, I collected functions from R packages for imputation [1,2,3], missing data visualisation [4] and statistical tests for missingness mechanisms [5,6]. I combined these functions to be conviniently accessible, provide plotting functions to visualise the results and demonstrate a workflow for phenotype imputation on two independent datasets [7,8,9]. 

At the moment, the functions lack detailed documentation, but as a starter, the vignettes 
[Workflow_yeastData](https://github.com/HannahVMeyer/PhenotypeImputation/blob/master/vignettes/Workflow_yeastData.pdf) and [Workflow_ratData](https://github.com/HannahVMeyer/PhenotypeImputation/blob/master/vignettes/Workflow_ratData.pdf) can be taken as guides on usage. 

# References
1. Dahl A, Iotchkova V, Baud A, Johansson Å, Gyllensten U, Soranzo N, et al. A multiple-phenotype imputation method for genetic studies. Nature Genetics. Nature Publishing Group; 2016;48: 466–472. doi:10.1038/ng.3513
2. Buuren S van, Groothuis-Oudshoorn K. mice : Multivariate Imputation by Chained Equations in R. Journal of Statistical Software. 2011;45: 1–67. doi:10.18637/jss.v045.i03
3. Little RJA, Rubin DB. Statistical analysis with missing data. 2nd ed. Balding DJ, Bloomfield P, Cressie NAC, Fisher NI, Johnstone IM, Kadane JB, et al., editors. New Jersey: John Wiley & Sons, Inc; 2002. p. 408.
4. https://CRAN.R-project.org/package=VIM
5. Little RJA. A Test of Missing Completely at Random for Multivariate Data with Missing Values. Source Journal of the American Statistical Association. 1988;83: 1198–1202.
6. https://CRAN.R-project.org/package=BaylorEdPsych
7. Bloom JS, Ehrenreich IM, Loo WT, Lite T-LV, Kruglyak L. Finding the sources of missing heritability in a yeast cross. Nature. 2013;494: 234–7. doi:10.1038/nature11867
8. Baud A, Hermsen R, Guryev V, Stridh P, Graham D, McBride MW, et al. Combined sequence-based and genetic mapping analysis of complex traits in outbred rats. Nature Genetics. Nature Publishing Group; 2013;45: 767–775. doi:10.1038/ng.2644
9. Baud A, Guryev V, Hummel O, Johannesson M, Hermsen R, Stridh P, et al. Genomes and phenomes of a population of outbred rats and its progenitors. Scientific Data. Nature Publishing Group; 2014;1: 140011. doi:10.1038/sdata.2014.11
