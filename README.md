# LC6Analysis

## LC6_Analysis.R	
This contains R script used for the anatomical analysis of LC6 and downstream neurons. Most EM figures (except connectivity matrices) were generated using this code. It's written to be executed in sequence. First load the initial data file "LC6_Analysis_data.RData", the run the script.

However, if you only want to generate a particular figure, then load the workspace file (contains already processed variables). You need to unzip them into one file first.
- LC6_Analysis_workspace.7z.001
- LC6_Analysis_workspace.7z.002
- LC6_Analysis_workspace.7z.003

After this, look for the figure of interest by searching, eg.  "Figure 6C". Then you need to run the whole section, which starts with, eg.

```
# target neuron RF wt by synapses ---------------------------------------------------------------------------------
```

## LC6_Analysis_connMatrix.R  
contains R code to generate the connectivity matrices.

"LC6_Analysis_connMatrix_data.RData"	is the data file that needs to be loaded.

## The rest are data and calibration files used by LC6_Analysis.R

