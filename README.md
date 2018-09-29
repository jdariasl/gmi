# gmi
Gaussian Mixture Imputation.

This is a MatLab implementation of an imputation algorithm based on a Gaussian Mixture Model. GMI introduce an additional step into the EM algorithm. 
In every iteration, GMI uses a GMM regressor to estimate the missing values, using the parameters of the model at that iteration.

The missing values in the data should be represented as NaN.

Run ScriptDemo.m to see an example.

This code depends on the NetLab toolbox. It is free available from:
http://www.aston.ac.uk/eas/research/groups/ncrg/resources/netlab/downloads/

