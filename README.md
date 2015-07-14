## FIQT
FDR Inverse Quantile Transformation
Fast and very accurate winner's curse adjustment for genome scans summary statistics.
Method under review at PLoS Genetics, 
**A simple yet accurate correction for winner's curse can predict signals discovered in much larger genome scans** [PDF](http://www.biorxiv.org/content/early/2015/05/13/019299)

## Outline

Given the extreme value distribution of scan statistics in the upper and lower tails and different distributions elsewhere, it is unclear (or, at least, very complicated) how to properly estimate the biases of, and thereby adjust, all statistics in a genome scan. However, it is extremely simple to adjust the p-values for the genome-wide multiple testing, e.g. using a FDR (as in this paper) or Holm procedure (which is always conservative). While the FDR procedure can be anticonservative for the extreme scenario of negatively correlated variables, in genetics, Z-scores are only locally correlated and we do not expect even these local correlations to be preponderantly negative. Thus, FDR is not expected to be anticonservative and this is the reason why it is commonly used in statistical genetics. For a first step let then 〖p_i〗^*,i=1,...,k be the FDR adjusted p-values. In the second step, we estimate the expected (adjusted) mean of Z-scores, (〖X_i〗^* ) ̂ ,i=1,...,k, by transforming the p-values in an upper tail standard normal quantile while keeping the sign of the original statistics. Or, in mathematical notation, 
(〖X_i〗^* ) ̂=sign(X_i ) 〖 Φ〗^(-1) (1-〖p_i〗^*/2). 


## R code

The 7 lines R function for computing this FIQT (winner's curse) adjustment is given below,
where:

**z** - association z-scores

**min.p** - minimum p-value admitted (to avoid zero p-values/adjusted p-values which give troubles with inverse cdf)

**min.p** - very large **z**s corresponding to **min.p** (i.e **z** > 37) are not adjusted, as their bias is essentially zero) 

```R
    FIQT<-function(z=z, min.p=10^-300){
        pvals<-2*pnorm(abs(z),low=F)
        pvals[pvals<min.p]<- min.p
        adj.pvals<-p.adjust(pvals,method="fdr")
        mu.z<-sign(z)*qnorm(adj.pvals/2,low=F)
        mu.z[abs(z)>qnorm(min.p/2,low=F)]<-z[abs(z)>qnorm(min.p/2,low=F)]
       mu.z
    }
```
