## FALCO: Factor AnaLysis in COrrelated data

The functions in this package implement the methods proposed in McKennan, 2020 (https://arxiv.org/abs/2009.11134), and perform factor analysis in high dimensional biological data. Let yg be an n-vector containing the expression or methylation of genomic unit g. To use this package, you must be able to express Var(yg) (the variance of yg) as Var(yg) = v\_{g1}B\_1 + ... + v\_{gb}B\_b, where v\_{g1},...,v\_{gb} are unknown scalars and B\_1,...,B\_b are known matrices that parametrize the relationship between samples. This encompasses nearly all modern gene expression and methylation data. Some examples include:

1) Unrelated samples. In this case, b=1 and B\_1=I_n is the identity matrix.

2) Individuals related through a kinship matrix. In this case, b=2, B\_1=I\_n, and B\_2=U, where U is the kinship matrix.

3) Multi-tissue data from unrelated individuals. For data with T tissues and arbitrary correlation structure, one can express V(yg) = sum\_{i=1}^{T(T+1)/2}v\_{gi}B\_i. One can simplify the correlation structure depending on the similarity between tissues.

4) Longitudinal data. For general longitudinal data with T time points, V(yg) can be expressed exactly as it is in 3). If one assumes the marginal variance for each sample is the same, this can be simplified to V(yg) = v\_{g1} I\_n + sum\_{i=2}^{T(T-1)/2 + 1}v\_{gi}B\_i.

The two primary functions are given below. 

### FALCO 

This implements FALCO (Algorithm 1 of McKennan, 2020), which estimates latent factors and loadings. Like PCA in data with unrelated samples, the factors are orthogonal to one another and are arranged in order of decreasing importance (i.e. decreasing variance explained). As demonstrated in McKennan, 2020, these behave like principal components, and can be used to perform quality control, identify latent patterns in the data, and de-noise the expression/methylation data matrix in eQTL and meQTL studies.

### CBCV_plus

This implements CBCV+ (Algorithm 2 of McKennan, 2020), which estimates K, the number of latent factors. If K is unspecified in FALCO, it is estimated using CBCV_plus.
