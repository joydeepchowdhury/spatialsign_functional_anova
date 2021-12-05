# The function SStestpvalue returns the p value of the SS test based on bootstrap or asymptotic
# implementation. It has the following arguments:
# (1) Data: An n-by-p matrix, where each row represents a functional observation on a fixed
# interval recorded on a equally spaced grid of length p. The first n_1 rows correspond to the n_1
# observations in the first group, the next n_2 rows correspond to the n_2 observations in the
# second group, and so on.
# (2) blocknumbers: A vector of length K giving the group sizes: (n_1, n_2, ..., n_K).
# (3) resolution: The inverse of the common distance between two consecutive grid points at which the
# observations are recorded.
# (4) asymporbootstrap: 1 for asymptotic implementation of the test, and 2 for bootstrap implementation
# of the test. Default value is 1.
# (5) bootstrapnum: Number of bootstrap samples for the bootstrap implementation. Default value is 1000.

              

SStestpvalue = function(Data, blocknumbers, resolution, asymporbootstrap = 1, bootstrapnum = 1000)
{
  require(SSinternal)
  
  if (asymporbootstrap == 1){
    SS = SSstatistic(Data, blocknumbers, resolution)
    
    SIGMA = SSCovOperator(Data, blocknumbers, resolution)
    
    highestportion = max(abs(SIGMA - t(SIGMA))) / max(abs(SIGMA))
    if (highestportion >= 1e-10){
      stop('SIGMA not symmetric!! Check Program and debug!!')
    }else{
      SIGMA = ( SIGMA + t(SIGMA) ) / 2
    }
    
    eigendecomposition = eigen(SIGMA, only.values = TRUE)
    eigenvalues = (eigendecomposition$values) / resolution
    
    replication = 1000
    K = length(blocknumbers)
    gridlength = ncol(Data)
    
    Standard.Normal.numbers = matrix(data = rnorm(K * gridlength * replication, mean = 0, sd = 1),
                                     nrow = (K * gridlength), ncol = replication)
    
    Chisquare.1.matrix = Standard.Normal.numbers^2
    simulatedvalues = 
      matrix(eigenvalues, nrow = 1, ncol = length(eigenvalues)) %*% Chisquare.1.matrix
    
    p.value.SS = sum(SS < simulatedvalues) / length(simulatedvalues)
    
    return(p.value.SS)
  }else{
    p.value.SS.bootstrap = SSbootstrappvalue(Data, blocknumbers, resolution, bootstrapnum)
    
    return(p.value.SS.bootstrap)
  }
}