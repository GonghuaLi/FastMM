# FastMM - an ultra-efficient toolbox for constraint-based metabolic modeling

## Description
Constraint-based metabolic modeling has become an important computational method in biomedical research, to understand metabolism related disease mechanism, to predict potential new drug targets antimetabolites, and biomarkers of complex diseases, such as cancer. The rapid accumulation of genomic and proteomic data from large scale disease studies, such as the Cancer Genome Atlas (TCGA), provides an unprecedented opportunity for personalized metabolic modeling of large number of patient samples. The currently available metabolic modeling software, such as COBRA 2.0, requires substantial computing time conducting a metabolic modeling, which limits its applications in large scale genome-wide analysis.

The FastMM was designed to provide a easy-to-use and efficient toolbox for constraint meatolic modeling.

In FastMM, most of the time-cost functions of metabolic modeling  are written by C. 
These functions including:
flux variability analysis, 
genome-wide single gene knockout analysis, 
genome-wide double gene knockout analysis, 
genome-wide single metabolite knockout analysis, 
genome-wide double metabolite knockout analysis,
and MCMC sampling.

**The wrapper of FastMM was written by matlab/octave, and can be perferctly merged into COBRA 2.0 toolbox.**

***FastMM*** is ***60~380 times faster*** than COBRA 2.0 in performing flux variability analysis and knockout analysis and returns consistent results. When applied to MCMC sampling, ***FastMM*** is also ***10~15 times faster*** than COBRA 2.0. 

## Dependancy:
1. opencobra
2. gurobi
3. cplex

## Installation
To install FastMM, users just download FastMM and type the following command in matlab/Octave:

```matlab
Install
```
And add the ./bin sub directory into your enriroment

## Test
Type the following command in matlab:

```matlab
Test
```

## How to use
### Basic use
Basic use provide a regular one-command solution for efficient metabolic modeling.
The user just configure the inputs in pars.txt and then type one commnand in matlab/Octave:

```matlab
FastMM
```

### Advanced use
All of the advanced functions of FastMM can be found the FastMM_mannual.pdf






