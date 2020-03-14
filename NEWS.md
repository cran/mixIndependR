# mixIndependR 0.2.1 -2020.03.11
# Fix Bugs
- Fix the error that "AlleleShare_Table" cannot deal with samples with odd sample size when "replicate=F"
# Add New Function
- Add a function named as "Prop_Pvalue" which can generate a bundle of p-values for one sample
# Update Description Files
- The Edition


# mixIndependR 0.2.0 -2020.02.01
# Fix bugs
- Fix format conflict in "AlleleShare_Table"
- Solve "NA" problem in “Dis_SimuChisq”
# Update Description Files
- Add "References" to "DistAlleleShare", "ExpProAlleleShare", "GenotypeFreq", "HWE_Fisher", and "RxpHetero"
- Correct "Usage" to include more information of input.
- Add "examples" to all R functions.
# New .R files
- "RxpHetero" to calculate average heterozygosity of observed or under Hardy-Weinberg Equilibrium
- "FreqAlleleShare" to build observed distribution for No. of shared alleles
- "FreqHetero" build observed distribution for No. of heterozygous loci
- "RealProAlleleShare" to calculate real density of shared alleles on each loci.
- "ComposPare_K" to generate a dataframe including observed and expected data of No. of heterozygous loci for easily plotting.
- "ComposPare_X" to generate a dataframe including observed and expected data of No. of sharing alleles for easily plotting.
# Update Old .R files
- Add column names to "AlleleShare_Table"


* Added a `NEWS.md` file to track changes to the package.
