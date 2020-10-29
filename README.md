# Predictive analysis of mtDNA variations

The initial goal of the project is first of all to verify if it exists any potential predictive relation between the CVIndex (reported in Cannon et al. 2002) and the GenBank appaearences for LSU (RNR2) and SSU (RNR1) rRNA mithocondrial sections. To do so human mithocondrial variations data for the RNR1 and RNR2 locus has been gathered from https://www.mitomap.org/foswiki/bin/view/MITOMAP/Resources, and CVIndex were extracted from Cannon et al. 2002.

The src folder contains a script to preprocess both control and MitoMap data to build and inteligible dataframe applicable in further analysis.

## Analysis

We have focused on 2 different types of regression models:
        
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
 
 All these analysis are gathered in r/analisis_rmd.Rmd