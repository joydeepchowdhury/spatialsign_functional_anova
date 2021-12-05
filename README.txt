The file SSinternal.rcpp contains some internal functions required in the function SStestpvalue.
One has to make a R package from SSinternal.rcpp to use the function SStestpvalue. The package
can be made by executing the R file makepackageSSinternal.R.


The function SStestpvalue returns the p value of the SS test based on bootstrap or asymptotic
implementation. It has the following arguments:
  (1) Data: An n-by-p matrix, where each row represents a functional observation on a fixed interval
              recorded on a equally spaced grid of length p. The first n_1 rows correspond to the n_1
	      observations in the first group, the next n_2 rows correspond to the n_2 observations
              in the second group, and so on.
  (2) blocknumbers: A vector of length K giving the group sizes: (n_1, n_2, ..., n_K).
  (3) resolution: The inverse of the common distance between two consecutive grid points at which the
	      observations are recorded.
  (4) asymporbootstrap: 1 for asymptotic implementation of the test, and 2 for bootstrap implementation
	      of the test. Default value is 1.
  (5) bootstrapnum: Number of bootstrap samples for the bootstrap implementation. Default value is 1000.
