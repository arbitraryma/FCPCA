# FCPCA
This repository provides the code for the fuzzy MTS clustering algorithm based on common principal component analysis (CPCA), which we name FCPCA.

To use this algorithm, please first run the related functions. You can find all the functions in the 'summary_fcpca' file. 
Some examples we used in the paper are also presented.

Recently, we also coded the method in Python; users are free to choose which version they like.

All results shown in our paper are generated using our R code. 

The method can be used to select the number of clusters (k) / fuzziness parameter (m) by using the CVI we introduced in our paper, and this is described by the S_value in our code. 
Users can make a grid search on the number of clusters and the fuzziness parameter, and choose the clustering result that generates the smallest CVI. 
