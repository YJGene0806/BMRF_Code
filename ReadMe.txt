[ReadMe.txt file]

There are three files (.R) for illustrating our proposed method "BMRF".

[BMRF_Function.R file]
The "BMRF_Function.R" file is the function for applying BMRF to estimate probabilistic network structure. 
This function can generate posterior samples for "Beta" and "Gamma" in the BMRF model, which can be applied to estimate probabilistic network and edge strength. 
The output of BMRF_Function is a "list", which is the default setting by using the "coda" package to store the posterior samples.

[Example.R file]
The "Example. R" file consists of a toy example to illustrate applying BMRF to construct network structure. 
This example shows how to use BMRF to estimate the probability of an edge in the network and the corresponding strength of partial correlation.

[Performance_Criteria.R file]
The "Performance_Criteria.R" file is the function for calculating the performance of estimating network structure by our proposed methods.  
By utilizing this function, one can generate those commonly used criteria when presenting estimation accuracy, such as F1-score, MCC score, etc.
