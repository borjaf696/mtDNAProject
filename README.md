# Predictive analysis of mtDNA variations

The initial goal of the project is first of all to verify if it exists any potential predictive relation between the CVIndex (reported in Cannon et al. 2002) and the GenBank appaearences for LSU (RNR2) and SSU (RNR1) rRNA mithocondrial sections. To do so human mithocondrial variations data for the RNR1 and RNR2 locus has been gathered from https://www.mitomap.org/foswiki/bin/view/MITOMAP/Resources, and CVIndex were extracted from Cannon et al. 2002.

The src folder contains a script to preprocess both control and MitoMap data to build and inteligible dataframe applicable in further analysis. 

## Data

### Control

Control information with the CVIndex from "Escherichia Coli" genome, and the "alignment" between LSU and SSU against this genome.

Dataframe inputs for R analysis, for every variation we store:

* CVTOT - CVIndex
* GB.Seqs - number of occurences in GenBank.
* Gb.Seqs.Corrected - number of occurences of a variation where the individual and the variation belong to the same haplogroup.
* Genomic - position in the genomic data.
* All set of haplogroups and the number of occurences. 

## Analysis

The core of the analysis uses two different regression styles:

1. Parametric:
    - Poisson
    - Zero inflated poisson
        * Geometric
        * Poisson
        * Binomial negative
    - PIG or Poisson Inverse Gaussian
    - SICHEL
    - DEL or Delaport
    - Logistic
2. Non Parametric:
    - Spline curves
        * B-splines
        * Smoothing splines
        * P-splines
    - Loess
    - Local Logistic regression
    
All the assayed models confirm the original hipotesis of "the higher number in the CVIndex the lower number of occurences in GenBank".

To extent the analysis haplogroup information (information related with the "race" of the individual being sequenced) was added, in the first place only the number of occurences were corrected removing a variation if the individual bearing it belongs to an haplogroup associated with the variation. Afterwards a more complex analysis was performed trying to look for some not yet associated haplogroups but there were not luck. 
 
All these analysis have been performed via R and can be found in r/analisis_rmd.Rmd in RMarkdown format
 
 ## Summary conclusions so far
 
 The data shows that the more CVIndex (which means that the region is highly preserved) the lower number of occurences in GenBank (which means less variations reported). These trend is kept for both SSU and LSU section. 