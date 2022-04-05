# ABC analysis of polyploid specation, chromosomal inheritance, and evolutionary histories 
Polyploidization and inheritance mode model estimation using Approximate Bayesian Computation. Developed for studying speciation in the North American gray treefrog complex *Hyla chrysoscelis/versicolor*

If used please cite the following work: https://academic.oup.com/mbe/article/39/2/msab316/6427635

This repository contains scripts and information on how to run an ABC analysis to estimate the mode of polyploid formation and the pattern of chromosomal inheritance using sequence data. As of now, this analysis only works for tetraploids. 

This analysis was modified from methods outlined in Roux and Pannel (2015)[1] and Leroy et al. (2017)[2]. This readme is also modified from examples laid out in their repositories: https://github.com/popgenomics/ABC_WGD, https://github.com/ThibaultLeroyFr/WhiteOaksABC

## Summary Statistics
The following summary statistics are calculated (modified from https://github.com/popgenomics/ABC_WGD)

| Statistics         | Description                                                                 |
|:-------------------|:----------------------------------------------------------------------------|
| __bialsites__          |  number of SNPs in the alignment                                            |
| __sf AB__               |  number of fixed differences between species A and B / locus length         |
| __sx A__                |  number of exclusively polymorphic positions in species A / locus length    |
| __sx B__                |  number of exclusively polymorphic positions in species B / locus length    |
| __ss AB__               |  number of shared biallelic positions between species A and B/ locus length |
| __successive ss__       |  maximum number of successive shared biallelic positions for a locus |
| __pi A__                |  Tajima’s Theta (pi) within species A                                            |
| __pi B__                |  Tajima’s Theta (pi) within species B                                            |
| __theta A__             |  Watterson’s Theta within species A                                          |
| __theta B__             |  Watterson’s Theta within species A                                          |
| __pearson_r_pi_AB__    |  correlation’s coefficient for pi over orthologs between A and B            |
| __pearson_r_theta_AB__ |  correlation’s coefficient for theta over orthologs between A and B         |
| __Dtaj A__              |  Tajima’s D for species A                                                   |
| __Dtaj B__              |  Tajima’s D for species B                                                   |
| __divAB__              |  raw divergence Dxy measured between A and B                                |
| __netdivAB__           |  net divergence Da measured between A and B                                 |
| __minDivAB__           |  smallest divergence measured between one individual from A and one from B  |
| __maxDivAB__           |  highest divergence measured between one individual from A and one from B   |
| __GminAB__             |  minimum divergence between one sequence from A and one from B (__minDivAB__ divided by the average __divAB__)                                                             |
| __GmaxAB__             |  maximum divergence between one sequence from A and one from B (__maxDivAB__ divided by the average __divAB__)                                                             |
| __FST__             |  FST between A and B compute as 1-(pi_A + pi_B) / (2 * pi_AB)                  |
| __D3a__             | absolute value of the D3 statistic with the diploid and two polyploid subgenomes  |
| __pearson_r_divAB_netDivAB__ |  correlation’s coefficient for divAB and netDivAB            |
| __pearson_r_divAB_FST__ |  correlation’s coefficient for divAB and FST        |
| __pearson_r_netDivAB_FST__    |  correlation’s coefficient for newDivAB and FST          |
| __ss_sf__ |  proportion of loci with both shared biallelic posistions and shared fixed differences   |
| __ss_noSf__ |  proportion of loci with shared biallelic posistions but no shared fixed differences   |
| __noSs_sf__ |  proportion of loci with no shared biallelic posistions but with shared fixed differences    |
| __noSs_noSf__ |  proportion of loci with neither shared biallelic posistions or shared fixed differences |


  
  

## Software Requirements:

Python
  - pypy (be sure to change dependency location in mscalc_wgd_v3.py)
  - numpy
  
R
  - seqinr
  - ape
  - nnet
  - abc/abcrf
  
msnsam[3]: https://github.com/rossibarra/msnsam
 
## Simulation Scripts

### run_ABC_polyploid_v3.py

This script is used to simulate sequence data for a specific model of polyploidization mode and chromosomal inheritance pattern. On the command line, running this would be as follows:

```
./run_ABC_polyploid_v3.py allo disomic migA F bpFile_EC_AE.txt 25000 1 40 F
```

The script takes 9 arguments:

1. mode of polyploid formation: allo, auto, polyP
  - (allo) allopolyploid 
  - (auto) autopolyploid from a sister lineage of the sampled diploid
  - (polyP) autopolyploid from the sampled diploid lineage 
2. pattern of chromosomal inheritance: disomic, heterosomic, tetrasomic
3. migration pattern: noMig, migA, migAB
  - (noMig) no migration between diploids and tetraploids
  - (migA) migration between the diploid and the A subgenome of the tetraploid. When the polyploidization mode is allopolyploid, this means migration between the diploid and the tetraploid subgenome inhereted from the sampled diploid
  - (migAB) migration between the diploid and both tetraploid subgenomes
4. One-way migration from diploid to tetraploid only: T/F
5. The bpfile, which contains information on the number of sampled loci, loci length, sampled individuals, and the clock-rate for the chosen loci
6. The number of multilocus simulations to run per repeat.
7. Beginning repeat number (Typically begins at 1, but on the chance the full number of simulations fail, it might be useful to start at a higher number so as not to start the simulation over)
8. End repeat number. Total number of repeats = {8} - {7} - 1. The total number of simulations is {6} * ({8}-{7} + 1). In the above example, the total number of multilocus simulations is 25000 * (40-1+1) = 1,000,000
9. Whether or not to keep or delete the joint site frequency spectrum file. This file is large, so if you aren't using it for a seperate analysis I recommend removing it 

This script will generate a folder titled "allo_disomic_migA" with 40 folders, each containing the following files:
 - ABCstat.txt: summary statistics for each multilocus simulation, one simulation per line
 - output.txt: prior values used to run each simulation, one simulation per line
 - ABCjsfs: the joint site frequency spectrum, one simulation per line (unless removed)
 
Within the script, it is also necessary to set the bounds for the prior distributions. You will need to set upper and lower limits for:

- n1: effective population size of the diploid
- n2: effective population size of the tetraploid
- nA: effective population size of the ancestral populations (used for both the ancestral pop. size prior to Tsplit and for the tetraploid prior to Twgd. 
- tau: bounds of the possible Tsplit and Twgd times
- TauDist: Use a uniform or exponential distribution for Tau. If exponential, the lower limit for Twgd is the first argument, and the mean of the distribution is the second argument. 
- M: scalar for the beta distribution for migration rate
- shape1: prior for the first shape parameter of the beta distribution for the number of migrants
- shape2: prior for the second shape parameter of the beta distribution for the number of migrants
- symModel: is migration symmetric or assymetric (asym/sym)
- MVariation: does migration vary across the genome? (homo/hetero)

### priorgen_wgd_geneflow_v2.py

The run_ABC_polyploid_v3.py script will pass arguments to this file to generate the values of each parameter for every simulation, this information will then get passed to msnsam to simulate the given model using TBS arguments

### mscalc_wgd_v3.py

After running each simulation, mscalc will calculate the values for each summary statistic to generate the ABCstat.txt files

## Generating the observed data

To generate statistics for the observed data, it is necessary to first convert your sequences into the ms format that can be read by mscalc_wgd_v3.py. This ensures all summary statistics are calculated in the same way as the simulated data and in the same format. Do not use summary statistics calculated by R packages or other software as their way of handling missing data or other aspects of genetic data may not be the same. Furthermore, it is just easier to convert to ms format.

To convert to the ms format, run the script convert_to_ms_from_fasta.R. The script only works for sequences in the fasta format, and it would be easiest to convert your sequence files to fasta. However, it should be easily modifiable if your data are in a different format. 

Additionally, you will need a population input file in csv format that contains 3 columns: individual names, population they came from, and whether or not to include that individual (set to 1 to include). See example file "pops_all_MinMax.csv" in the examples folder

This script will automatically generate your bpfile as well as generate the ms file to calculate the observed data along with a table detailing which loci were used.

To generate the observed sumstat file, make a folder in the directory with the bpfile and ms file titled "Obs" (in this example) and run the following code: 
```
cat EC_AE_msFile.txt | ./mscalc_wgd_v3.py Obs bpFile_EC_AE.txt
```
That will generate the ABCjsfs.txt file for your data and the ABCstat.txt file. I recommend changing the name so you can keep track of these files.

## Running the final analysis

ABC_nnet_FULL_.R will run the final analysis using the ABC neural network approach. It is designed in a way such that all your simulated and observed data for two populations of interest will be in a folder *popA_popB* (e.g. EC_AE in the example), with the bpfile as bpFile_*popA_popB*.txt, and your observed data as *popA_popB*_ObsData.txt. This way you only have to change a single line in the R script to automatically run the analysis. It is set up to run 5 seperate analysis for 5 different models(noMig, migA, migAB, migA_1W, migAB_1W). You will need to modify it if you run different models. 

You will need to set the number of replicates closest to the observed data to extract with the weighted Epanechnikov kernel (EpReplicates), the number of models for the individual analysis (nModels), and the number of models for the final analysis (nModelsFinal). The percentage of the replicate simulations closest to the observed values of the summary statistics is then (EpReplicates/nModels * the number of multilocus simulations per model). In the example EpReplicates is 4500, so for 9 models and 1,000,000 multilocus simulations this is 0.5% of the total number of replicates. This percentage is maintained for the final analysis by setting nModelsFinal. 

You will also need to set the number of replicates for the neural network analysis as well as the number of folders for each dataset. In the example, there are 20 nnet replicates and 40 folders. 

This analysis will generate a text file for all individual analyses and the final analyses with the model probabilities on each line in the order input into the nnet analysis. In the example, the order is: allo_disomic, allo_heterosomic, allo_tetrasomic, auto_disomic, auto_heterosomic, auto_tetrasomic, polyP_disomic, polyP_heterosomic, polyP_tetrasomic. summarize by loading the text file into R and taking the column means (colMeans()). 

Finally, it will also automatically generate posterior distributions for the top 3 models in individual text files for each posterior neural network replicate. to view the posterior distribution for each parameter(example Tsplit) run something like:

```
setwd("/gpfs/research/scratch/wwb15/ABC_WGD/50Loci/Exponential/EC_AE/")
postSum = NULL
nPostReps <- 10
chosenModel <- "allo_disomic_migAB_1W"
for(i in 1:nPostReps){
  postSum <- rbind(read.table(paste(c("Posterior_EC_AE",chosenModel,i), collapse = "")))
}
colnames(postSum) <- colnames(read.table(paste(c("Posterior_EC_AE",chosenModel,"_Header"), collapse = ""), header = T))
hist(postSum$Twgd)

```
To compare to the prior, import the output.txt files in each folder using a for loop for that model and generate similar histograms.

## Important information

It's important to note that prior values are set with the assumption that the effective population size multiplier is 100,000. That means, for instance, if we set our upper bound of n1 to 5, that means it is actually 500,000. To change this, you will need to edit the n0 value at the top of the priorgen python script as well as in the convert_to_ms_from_fasta.R file. Furthermore, this affects interpretation of Tsplit and Twgd times, which are multiplied by 4 * n0. So, if you set the upper limit of tau to 1, that means it is actually 400,000 generations (since our clock rate is in substitutions/yr it can be approximated as 400,000 ybp). 

The bpFile has 5 lines. The first is a comment line which you can change to anything you like. The second and third rows are the sampled number of individuals from A and B respectively. The third and fourth rows are information on theta and rho, which are 4 * n0 * clockRate

## Final Notes

I have tried to make this type of analysis as user friendly as possible given the number of different steps that need to be taken and the amount of customization each dataset requires. However, I expect individual analyses will take some tweaking to work properly. If you have any questions, feel free to email me at wbooker14@gmail.com and I will do my best to respond to you in a timely manner. 

There are also scripts included for running an ABC random forest analysis, however this was not used for our study apart from generating LDA plots. An example script for this is in the examples folder. 


## Citations:

[1] Roux, C. and Pannell, J. R. (2015). Inferring the mode of origin of polyploid species from next-generation sequence data. Molecular Ecology 24(5), 1047–1059.

[2] Leroy, T., Roux, C., Villate, L., Bodénès, C., Romiguier, J., Paiva, J. A. P., Dossat, C., Aury, J.-M., Plomion, C. and Kremer, A. (2017). Extensive recent secondary contacts between four European white oak species. The New Phytologist 214(2), 865–878.

[3] Ross-Ibarra J, Wright SI, Foxe JP, Kawabe A, DeRose-Wilson L, et al. (2008). Patterns of Polymorphism and Demographic History in Natural Populations of Arabidopsis lyrata . PLoS ONE 3(6): e2411
