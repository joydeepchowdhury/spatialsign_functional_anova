#include <Rcpp.h>
using namespace Rcpp;

// Computing Covariance Operator of SS statistic

// [[Rcpp::export]]

NumericMatrix SSCovOperator(NumericMatrix Data, NumericVector blocksizes, double resolution) {
  int n = Data.nrow(), p = Data.ncol();
  int numblocks = blocksizes.size();
  double spatialsumcoord;
  double vector[p], vectortemp[p], spatialdistvectorsum[p];
  double spatialdistarray[n][numblocks][p];
  double spatialdistav[numblocks][numblocks][p];
  
  NumericMatrix returnmatrix((numblocks * p), (numblocks * p));
  
  int i, j, k, l, currentblocksize, counter;
  double distance;
  
  for (i = 0; i < n; ++i){
    for (l = 0; l < p; ++l){
      vector[l] = Data(i, l);
    }
    
    counter = 0;
    for (j = 0; j < numblocks; ++j){
      currentblocksize = blocksizes[j];
      
      for (l = 0; l < p; ++l){
        spatialdistvectorsum[l] = 0;
      }
      for (k = 0; k < currentblocksize; ++k){
        for (l = 0; l < p; ++l){
          vectortemp[l] = Data(counter, l);
        }
        counter = counter + 1;
        
        // L2 distance calculation between vector and vectortemp: sqrt(sum(vector - vectortemp))
        distance = 0;
        for (l = 0; l < p; ++l){
          distance = distance + (pow(vector[l] - vectortemp[l], 2) / resolution);
        }
        distance = sqrt(distance); // Euclidean distance computation complete
        
        if (distance > 1e-15){
          for (l = 0; l < p; ++l){
            spatialdistvectorsum[l] = spatialdistvectorsum[l] + ((vector[l] - vectortemp[l]) / distance);
          }
        }
      }
      
      for (l = 0; l < p; ++l){
        spatialdistarray[i][j][l] = spatialdistvectorsum[l] / currentblocksize;
      }
    }
  }
  
  // spatialdistav computation
  for (j = 0; j < numblocks; ++j){
    for (k = 0; k < p; ++k){
      counter = 0;
      for (i = 0; i < numblocks; ++i){
        currentblocksize = blocksizes[i];
        
        spatialsumcoord = 0;
        for (l = 0; l < currentblocksize; ++l){
          spatialsumcoord = spatialsumcoord + spatialdistarray[counter][j][k];
          
          counter = counter + 1;
        }
        
        spatialdistav[i][j][k] = spatialsumcoord / currentblocksize;
      }
    }
  }
  
  // Covariance operator computation
  double CovMatrix[(numblocks * p)][(numblocks * p)];
  
  int k1, k2, l1,l2;
  int cumsumblocksizes[(numblocks + 1)];
  double lambda[numblocks];
  
  double firstterm[p][p], secondterm[p][p];
  double block1[p][p], block2[p][p];
  double block3[p][p], block4[p][p];
  
  for (k = 0; k <= numblocks; ++k){
    cumsumblocksizes[k] = 0;
    for (l = 0; l < k; ++l){
      cumsumblocksizes[k] = cumsumblocksizes[k] + blocksizes[l];
    }
  }
  
  for (k = 0; k < numblocks; ++k){
    lambda[k] = blocksizes[k] / n;
  }
  
  for (k1 = 0; k1 < numblocks; ++k1){
    for (k2 = 0; k2 < numblocks; ++k2){
      // Computating first composite term in the expression of sigma(k1, k2)
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          firstterm[i][j] = 0;
        }
      }
      
      for (l = 0; l < numblocks; ++l){
        // Initialising all elements of block1, block2 and block3 to 0
        for (i = 0; i < p; ++i){
          for (j = 0; j < p; ++j){
            block1[i][j] = 0;
            block2[i][j] = 0;
            block3[i][j] = 0;
          }
        }
        
        // Computing first block
        for (i = 0; i < p; ++i){
          for (j = 0; j < p; ++j){
            counter = cumsumblocksizes[l];
            
            while (counter < cumsumblocksizes[l + 1]){
              block1[i][j] = block1[i][j] + 
                (spatialdistarray[counter][k1][i] * spatialdistarray[counter][k2][j]);
              counter = counter + 1;
            }
            block1[i][j] = block1[i][j] / blocksizes[l];
            
            block1[i][j] = block1[i][j] - (spatialdistav[l][k1][i] * spatialdistav[l][k2][j]);
          }
        }
        
        // Computing second block
        for (i = 0; i < p; ++i){
          for (j = 0; j < p; ++j){
            counter = cumsumblocksizes[k1];
            
            while (counter < cumsumblocksizes[k1 + 1]){
              block2[i][j] = block2[i][j] + 
                (spatialdistarray[counter][l][i] * spatialdistarray[counter][k2][j]);
              counter = counter + 1;
            }
            block2[i][j] = block2[i][j] / blocksizes[k1];
            
            block2[i][j] = block2[i][j] - (spatialdistav[k1][l][i] * spatialdistav[k1][k2][j]);
            
            block2[i][j] = - block2[i][j];
          }
        }
        
        // Computing third block
        for (i = 0; i < p; ++i){
          for (j = 0; j < p; ++j){
            counter = cumsumblocksizes[k2];
            
            while (counter < cumsumblocksizes[k2 + 1]){
              block3[i][j] = block3[i][j] + 
                (spatialdistarray[counter][k1][i] * spatialdistarray[counter][l][j]);
              counter = counter + 1;
            }
            block3[i][j] = block3[i][j] / blocksizes[k2];
            
            block3[i][j] = block3[i][j] - (spatialdistav[k2][k1][i] * spatialdistav[k2][l][j]);
            
            block3[i][j] = - block3[i][j];
          }
        }
        
        // Adding block1, block2, block3 over l
        for (i = 0; i < p; ++i){
          for (j = 0; j < p; ++j){
            firstterm[i][j] = firstterm[i][j] + 
              ( lambda[l] * (block1[i][j] + block2[i][j] + block3[i][j]) );
          }
        }
      }
      
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          firstterm[i][j] = sqrt(lambda[k1] * lambda[k2]) * firstterm[i][j];
        }
      }
      
      // Initialising and computing second composite term in the expression of sigma(k1, k2)
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          secondterm[i][j] = 0;
        }
      }
      
      if (k1 == k2){
        for (l1 = 0; l1 < numblocks; ++l1){
          for (l2 = 0; l2 < numblocks; ++l2){
            // Initialising all elements of block4 to 0
            for (i = 0; i < p; ++i){
              for (j = 0; j < p; ++j){
                block4[i][j] = 0;
              }
            }
            
            // Computing fourth block
            for (i = 0; i < p; ++i){
              for (j = 0; j < p; ++j){
                counter = cumsumblocksizes[k1];
                
                while (counter < cumsumblocksizes[k1 + 1]){
                  block4[i][j] = block4[i][j] + 
                    (spatialdistarray[counter][l1][i] * spatialdistarray[counter][l2][j]);
                  counter = counter + 1;
                }
                block4[i][j] = block4[i][j] / blocksizes[k1];
                
                block4[i][j] = block4[i][j] - (spatialdistav[k1][l1][i] * spatialdistav[k1][l2][j]);
              }
            }
            
            // Adding block4 over l1, l2
            for (i = 0; i < p; ++i){
              for (j = 0; j < p; ++j){
                secondterm[i][j] = secondterm[i][j] + ( lambda[l1] * lambda[l2] * block4[i][j] );
              }
            }
          }
        }
      }
      
      // Filling up sigma(k1, k2) in Sigma
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          CovMatrix[((k1 * p) + i)][((k2 * p) + j)] = 
            firstterm[i][j] + secondterm[i][j];
        }
      }
    }
  }
  
  for (i = 0; i < (numblocks * p); ++i){
    for (j = 0; j < (numblocks * p); ++j){
      returnmatrix(i, j) = CovMatrix[i][j];
    }
  }
  
  return returnmatrix;
}


// Computing Covariance Operator (adjusted) of SS statistic

// [[Rcpp::export]]

NumericMatrix SSCovOperatoradjusted(NumericMatrix Data, NumericVector blocksizes, double resolution) {
  int n = Data.nrow(), p = Data.ncol();
  int numblocks = blocksizes.size();
  double spatialsumcoord;
  double vector[p], vectortemp[p], spatialdistvectorsum[p];
  double spatialdistarray[n][numblocks][p];
  double spatialdistav[numblocks][numblocks][p];
  
  NumericMatrix returnmatrix((numblocks * p), (numblocks * p));
  
  int i, j, k, l, currentblocksize, counter;
  double distance;
  
  for (i = 0; i < n; ++i){
    for (l = 0; l < p; ++l){
      vector[l] = Data(i, l);
    }
    
    counter = 0;
    for (j = 0; j < numblocks; ++j){
      currentblocksize = blocksizes[j];
      
      for (l = 0; l < p; ++l){
        spatialdistvectorsum[l] = 0;
      }
      for (k = 0; k < currentblocksize; ++k){
        for (l = 0; l < p; ++l){
          vectortemp[l] = Data(counter, l);
        }
        counter = counter + 1;
        
        // L2 distance calculation between vector and vectortemp: sqrt(sum(vector - vectortemp))
        distance = 0;
        for (l = 0; l < p; ++l){
          distance = distance + (pow(vector[l] - vectortemp[l], 2) / resolution);
        }
        distance = sqrt(distance); // Euclidean distance computation complete
        
        if (distance > 1e-15){
          for (l = 0; l < p; ++l){
            spatialdistvectorsum[l] = spatialdistvectorsum[l] + ((vector[l] - vectortemp[l]) / distance);
          }
        }
      }
      
      for (l = 0; l < p; ++l){
        // spatialdistarray[i][j][l] = spatialdistvectorsum[l] / (currentblocksize - 1);   /////////////////
        spatialdistarray[i][j][l] = spatialdistvectorsum[l] / (currentblocksize);   /////////////////
      }
    }
  }
  
  // spatialdistav computation
  for (j = 0; j < numblocks; ++j){
    for (k = 0; k < p; ++k){
      counter = 0;
      for (i = 0; i < numblocks; ++i){
        currentblocksize = blocksizes[i];
        
        spatialsumcoord = 0;
        for (l = 0; l < currentblocksize; ++l){
          spatialsumcoord = spatialsumcoord + spatialdistarray[counter][j][k];
          
          counter = counter + 1;
        }
        
        // spatialdistav[i][j][k] = spatialsumcoord / (currentblocksize - 1);              /////////////////
        spatialdistav[i][j][k] = spatialsumcoord / (currentblocksize);              /////////////////
      }
    }
  }
  
  // Covariance operator computation
  double CovMatrix[(numblocks * p)][(numblocks * p)];
  
  int k1, k2, l1,l2;
  int cumsumblocksizes[(numblocks + 1)];
  double lambda[numblocks];
  
  double firstterm[p][p], secondterm[p][p];
  double block1[p][p], block2[p][p];
  double block3[p][p], block4[p][p];
  
  for (k = 0; k <= numblocks; ++k){
    cumsumblocksizes[k] = 0;
    for (l = 0; l < k; ++l){
      cumsumblocksizes[k] = cumsumblocksizes[k] + blocksizes[l];
    }
  }
  
  for (k = 0; k < numblocks; ++k){
    lambda[k] = blocksizes[k] / n;
  }
  
  for (k1 = 0; k1 < numblocks; ++k1){
    for (k2 = 0; k2 < numblocks; ++k2){
      // Computating first composite term in the expression of sigma(k1, k2)
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          firstterm[i][j] = 0;
        }
      }
      
      for (l = 0; l < numblocks; ++l){
        // Initialising all elements of block1, block2 and block3 to 0
        for (i = 0; i < p; ++i){
          for (j = 0; j < p; ++j){
            block1[i][j] = 0;
            block2[i][j] = 0;
            block3[i][j] = 0;
          }
        }
        
        // Computing first block
        for (i = 0; i < p; ++i){
          for (j = 0; j < p; ++j){
            counter = cumsumblocksizes[l];
            
            while (counter < cumsumblocksizes[l + 1]){
              block1[i][j] = block1[i][j] + 
                (spatialdistarray[counter][k1][i] * spatialdistarray[counter][k2][j]);
              counter = counter + 1;
            }
            block1[i][j] = block1[i][j] / (blocksizes[l] - 1);                        /////////////////
            
            block1[i][j] = block1[i][j] - (spatialdistav[l][k1][i] * spatialdistav[l][k2][j]);
          }
        }
        
        // Computing second block
        for (i = 0; i < p; ++i){
          for (j = 0; j < p; ++j){
            counter = cumsumblocksizes[k1];
            
            while (counter < cumsumblocksizes[k1 + 1]){
              block2[i][j] = block2[i][j] + 
                (spatialdistarray[counter][l][i] * spatialdistarray[counter][k2][j]);
              counter = counter + 1;
            }
            block2[i][j] = block2[i][j] / (blocksizes[k1] - 1);                       /////////////////
            
            block2[i][j] = block2[i][j] - (spatialdistav[k1][l][i] * spatialdistav[k1][k2][j]);
            
            block2[i][j] = - block2[i][j];
          }
        }
        
        // Computing third block
        for (i = 0; i < p; ++i){
          for (j = 0; j < p; ++j){
            counter = cumsumblocksizes[k2];
            
            while (counter < cumsumblocksizes[k2 + 1]){
              block3[i][j] = block3[i][j] + 
                (spatialdistarray[counter][k1][i] * spatialdistarray[counter][l][j]);
              counter = counter + 1;
            }
            block3[i][j] = block3[i][j] / (blocksizes[k2] - 1);                       /////////////////
            
            block3[i][j] = block3[i][j] - (spatialdistav[k2][k1][i] * spatialdistav[k2][l][j]);
            
            block3[i][j] = - block3[i][j];
          }
        }
        
        // Adding block1, block2, block3 over l
        for (i = 0; i < p; ++i){
          for (j = 0; j < p; ++j){
            firstterm[i][j] = firstterm[i][j] + 
              ( lambda[l] * (block1[i][j] + block2[i][j] + block3[i][j]) );
          }
        }
      }
      
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          firstterm[i][j] = sqrt(lambda[k1] * lambda[k2]) * firstterm[i][j];
        }
      }
      
      // Initialising and computing second composite term in the expression of sigma(k1, k2)
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          secondterm[i][j] = 0;
        }
      }
      
      if (k1 == k2){
        for (l1 = 0; l1 < numblocks; ++l1){
          for (l2 = 0; l2 < numblocks; ++l2){
            // Initialising all elements of block4 to 0
            for (i = 0; i < p; ++i){
              for (j = 0; j < p; ++j){
                block4[i][j] = 0;
              }
            }
            
            // Computing fourth block
            for (i = 0; i < p; ++i){
              for (j = 0; j < p; ++j){
                counter = cumsumblocksizes[k1];
                
                while (counter < cumsumblocksizes[k1 + 1]){
                  block4[i][j] = block4[i][j] + 
                    (spatialdistarray[counter][l1][i] * spatialdistarray[counter][l2][j]);
                  counter = counter + 1;
                }
                block4[i][j] = block4[i][j] / (blocksizes[k1] - 1);                  /////////////////
                
                block4[i][j] = block4[i][j] - (spatialdistav[k1][l1][i] * spatialdistav[k1][l2][j]);
              }
            }
            
            // Adding block4 over l1, l2
            for (i = 0; i < p; ++i){
              for (j = 0; j < p; ++j){
                secondterm[i][j] = secondterm[i][j] + ( lambda[l1] * lambda[l2] * block4[i][j] );
              }
            }
          }
        }
      }
      
      // Filling up sigma(k1, k2) in Sigma
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          CovMatrix[((k1 * p) + i)][((k2 * p) + j)] = 
            firstterm[i][j] + secondterm[i][j];
        }
      }
    }
  }
  
  for (i = 0; i < (numblocks * p); ++i){
    for (j = 0; j < (numblocks * p); ++j){
      returnmatrix(i, j) = CovMatrix[i][j];
    }
  }
  
  return returnmatrix;
}


// Computing trace of Covariance Operator (adjusted) of SS statistic

// [[Rcpp::export]]

double SStraceCovOperatoradjusted(NumericMatrix Data, NumericVector blocksizes, double resolution) {
  int n = Data.nrow(), p = Data.ncol();
  int numblocks = blocksizes.size();
  double spatialsumcoord;
  double vector[p], vectortemp[p], spatialdistvectorsum[p];
  double spatialdistarray[n][numblocks][p];
  double spatialdistav[numblocks][numblocks][p];
  
  double traceCov = 0;
  
  int i, j, k, l, currentblocksize, counter;
  double distance;
  
  for (i = 0; i < n; ++i){
    for (l = 0; l < p; ++l){
      vector[l] = Data(i, l);
    }
    
    counter = 0;
    for (j = 0; j < numblocks; ++j){
      currentblocksize = blocksizes[j];
      
      for (l = 0; l < p; ++l){
        spatialdistvectorsum[l] = 0;
      }
      for (k = 0; k < currentblocksize; ++k){
        for (l = 0; l < p; ++l){
          vectortemp[l] = Data(counter, l);
        }
        counter = counter + 1;
        
        // L2 distance calculation between vector and vectortemp: sqrt(sum(vector - vectortemp))
        distance = 0;
        for (l = 0; l < p; ++l){
          distance = distance + (pow(vector[l] - vectortemp[l], 2) / resolution);
        }
        distance = sqrt(distance); // Euclidean distance computation complete
        
        if (distance > 1e-15){
          for (l = 0; l < p; ++l){
            spatialdistvectorsum[l] = spatialdistvectorsum[l] + ((vector[l] - vectortemp[l]) / distance);
          }
        }
      }
      
      for (l = 0; l < p; ++l){
        // spatialdistarray[i][j][l] = spatialdistvectorsum[l] / (currentblocksize - 1);   /////////////////
        spatialdistarray[i][j][l] = spatialdistvectorsum[l] / (currentblocksize);   /////////////////
      }
    }
  }
  
  // spatialdistav computation
  for (j = 0; j < numblocks; ++j){
    for (k = 0; k < p; ++k){
      counter = 0;
      for (i = 0; i < numblocks; ++i){
        currentblocksize = blocksizes[i];
        
        spatialsumcoord = 0;
        for (l = 0; l < currentblocksize; ++l){
          spatialsumcoord = spatialsumcoord + spatialdistarray[counter][j][k];
          
          counter = counter + 1;
        }
        
        // spatialdistav[i][j][k] = spatialsumcoord / (currentblocksize - 1);              /////////////////
        spatialdistav[i][j][k] = spatialsumcoord / (currentblocksize);              /////////////////
      }
    }
  }
  
  // Covariance operator trace computation
  double tempmatrix[p][p];
  
  int k1, k2, l1,l2;
  int cumsumblocksizes[(numblocks + 1)];
  double lambda[numblocks];
  
  double firstterm[p][p], secondterm[p][p];
  double block1[p][p], block2[p][p];
  double block3[p][p], block4[p][p];
  
  for (k = 0; k <= numblocks; ++k){
    cumsumblocksizes[k] = 0;
    for (l = 0; l < k; ++l){
      cumsumblocksizes[k] = cumsumblocksizes[k] + blocksizes[l];
    }
  }
  
  for (k = 0; k < numblocks; ++k){
    lambda[k] = blocksizes[k] / n;
  }
  
  for (k1 = 0; k1 < numblocks; ++k1){
    k2 = k1;
    
    // Computating first composite term in the expression of sigma(k1, k2)
    for (i = 0; i < p; ++i){
      for (j = 0; j < p; ++j){
        firstterm[i][j] = 0;
      }
    }
    
    for (l = 0; l < numblocks; ++l){
      // Initialising all elements of block1, block2 and block3 to 0
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          block1[i][j] = 0;
          block2[i][j] = 0;
          block3[i][j] = 0;
        }
      }
      
      // Computing first block
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          counter = cumsumblocksizes[l];
          
          while (counter < cumsumblocksizes[l + 1]){
            block1[i][j] = block1[i][j] + 
              (spatialdistarray[counter][k1][i] * spatialdistarray[counter][k2][j]);
            counter = counter + 1;
          }
          block1[i][j] = block1[i][j] / (blocksizes[l] - 1);                        /////////////////
          
          block1[i][j] = block1[i][j] - (spatialdistav[l][k1][i] * spatialdistav[l][k2][j]);
        }
      }
      
      // Computing second block
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          counter = cumsumblocksizes[k1];
          
          while (counter < cumsumblocksizes[k1 + 1]){
            block2[i][j] = block2[i][j] + 
              (spatialdistarray[counter][l][i] * spatialdistarray[counter][k2][j]);
            counter = counter + 1;
          }
          block2[i][j] = block2[i][j] / (blocksizes[k1] - 1);                       /////////////////
          
          block2[i][j] = block2[i][j] - (spatialdistav[k1][l][i] * spatialdistav[k1][k2][j]);
          
          block2[i][j] = - block2[i][j];
        }
      }
      
      // Computing third block
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          counter = cumsumblocksizes[k2];
          
          while (counter < cumsumblocksizes[k2 + 1]){
            block3[i][j] = block3[i][j] + 
              (spatialdistarray[counter][k1][i] * spatialdistarray[counter][l][j]);
            counter = counter + 1;
          }
          block3[i][j] = block3[i][j] / (blocksizes[k2] - 1);                       /////////////////
          
          block3[i][j] = block3[i][j] - (spatialdistav[k2][k1][i] * spatialdistav[k2][l][j]);
          
          block3[i][j] = - block3[i][j];
        }
      }
      
      // Adding block1, block2, block3 over l
      for (i = 0; i < p; ++i){
        for (j = 0; j < p; ++j){
          firstterm[i][j] = firstterm[i][j] + 
            ( lambda[l] * (block1[i][j] + block2[i][j] + block3[i][j]) );
        }
      }
    }
    
    for (i = 0; i < p; ++i){
      for (j = 0; j < p; ++j){
        firstterm[i][j] = sqrt(lambda[k1] * lambda[k2]) * firstterm[i][j];
      }
    }
    
    // Initialising and computing second composite term in the expression of sigma(k1, k2)
    for (i = 0; i < p; ++i){
      for (j = 0; j < p; ++j){
        secondterm[i][j] = 0;
      }
    }
    
    if (k1 == k2){
      for (l1 = 0; l1 < numblocks; ++l1){
        for (l2 = 0; l2 < numblocks; ++l2){
          // Initialising all elements of block4 to 0
          for (i = 0; i < p; ++i){
            for (j = 0; j < p; ++j){
              block4[i][j] = 0;
            }
          }
          
          // Computing fourth block
          for (i = 0; i < p; ++i){
            for (j = 0; j < p; ++j){
              counter = cumsumblocksizes[k1];
              
              while (counter < cumsumblocksizes[k1 + 1]){
                block4[i][j] = block4[i][j] + 
                  (spatialdistarray[counter][l1][i] * spatialdistarray[counter][l2][j]);
                counter = counter + 1;
              }
              block4[i][j] = block4[i][j] / (blocksizes[k1] - 1);                  /////////////////
              
              block4[i][j] = block4[i][j] - (spatialdistav[k1][l1][i] * spatialdistav[k1][l2][j]);
            }
          }
          
          // Adding block4 over l1, l2
          for (i = 0; i < p; ++i){
            for (j = 0; j < p; ++j){
              secondterm[i][j] = secondterm[i][j] + ( lambda[l1] * lambda[l2] * block4[i][j] );
            }
          }
        }
      }
    }
    
    // Filling up tempmatrix
    for (i = 0; i < p; ++i){
      for (j = 0; j < p; ++j){
        tempmatrix[i][j] = firstterm[i][j] + secondterm[i][j];
      }
    }
    
    // Adding trace of tempmatrix to traceCov
    for (i = 0; i < p; ++i){
      traceCov = traceCov + tempmatrix[i][i];
    }
  }
  
  return traceCov;
}


// Computing fast trace of Covariance Operator (adjusted) of SS statistic

// [[Rcpp::export]]

double SSfasttraceCovOperatoradjusted(NumericMatrix Data, NumericVector blocksizes, double resolution) {
  int n = Data.nrow(), p = Data.ncol();
  int numblocks = blocksizes.size();
  double spatialsumcoord;
  double vector[p], vectortemp[p], spatialdistvectorsum[p];
  double spatialdistarray[n][numblocks][p];
  double spatialdistav[numblocks][numblocks][p];
  
  double traceCov = 0;
  
  int i, j, k, l, currentblocksize, counter;
  double distance;
  
  for (i = 0; i < n; ++i){
    for (l = 0; l < p; ++l){
      vector[l] = Data(i, l);
    }
    
    counter = 0;
    for (j = 0; j < numblocks; ++j){
      currentblocksize = blocksizes[j];
      
      for (l = 0; l < p; ++l){
        spatialdistvectorsum[l] = 0;
      }
      for (k = 0; k < currentblocksize; ++k){
        for (l = 0; l < p; ++l){
          vectortemp[l] = Data(counter, l);
        }
        counter = counter + 1;
        
        // L2 distance calculation between vector and vectortemp: sqrt(sum(vector - vectortemp))
        distance = 0;
        for (l = 0; l < p; ++l){
          distance = distance + (pow(vector[l] - vectortemp[l], 2) / resolution);
        }
        distance = sqrt(distance); // Euclidean distance computation complete
        
        if (distance > 1e-15){
          for (l = 0; l < p; ++l){
            spatialdistvectorsum[l] = spatialdistvectorsum[l] + ((vector[l] - vectortemp[l]) / distance);
          }
        }
      }
      
      for (l = 0; l < p; ++l){
        // spatialdistarray[i][j][l] = spatialdistvectorsum[l] / (currentblocksize - 1);   /////////////////
        spatialdistarray[i][j][l] = spatialdistvectorsum[l] / (currentblocksize);   /////////////////
      }
    }
  }
  
  // spatialdistav computation
  for (j = 0; j < numblocks; ++j){
    for (k = 0; k < p; ++k){
      counter = 0;
      for (i = 0; i < numblocks; ++i){
        currentblocksize = blocksizes[i];
        
        spatialsumcoord = 0;
        for (l = 0; l < currentblocksize; ++l){
          spatialsumcoord = spatialsumcoord + spatialdistarray[counter][j][k];
          
          counter = counter + 1;
        }
        
        // spatialdistav[i][j][k] = spatialsumcoord / (currentblocksize - 1);              /////////////////
        spatialdistav[i][j][k] = spatialsumcoord / (currentblocksize);              /////////////////
      }
    }
  }
  
  // Covariance operator trace computation
  int l1,l2;
  int cumsumblocksizes[(numblocks + 1)];
  double lambda[numblocks];
  
  double firsttermdiag[p], secondtermdiag[p];
  double block1diag[p], block2diag[p];
  double block3diag[p], block4diag[p];
  
  for (k = 0; k <= numblocks; ++k){
    cumsumblocksizes[k] = 0;
    for (l = 0; l < k; ++l){
      cumsumblocksizes[k] = cumsumblocksizes[k] + blocksizes[l];
    }
  }
  
  for (k = 0; k < numblocks; ++k){
    lambda[k] = blocksizes[k] / n;
  }
  
  for (k = 0; k < numblocks; ++k){
    // Computating first composite term diagonal in the expression of sigma(k, k)
    for (i = 0; i < p; ++i){
      firsttermdiag[i] = 0;
    }
    
    for (l = 0; l < numblocks; ++l){
      // Initialising all elements of diagonals of block1, block2 and block3 to 0
      for (i = 0; i < p; ++i){
        block1diag[i] = 0;
        block2diag[i] = 0;
        block3diag[i] = 0;
      }
      
      // Computing diagonal of first block
      for (i = 0; i < p; ++i){
        counter = cumsumblocksizes[l];
        
        while (counter < cumsumblocksizes[l + 1]){
          block1diag[i] = block1diag[i] + 
            (spatialdistarray[counter][k][i] * spatialdistarray[counter][k][i]);
          counter = counter + 1;
        }
        block1diag[i] = block1diag[i] / (blocksizes[l] - 1);                        /////////////////
        
        block1diag[i] = block1diag[i] - (spatialdistav[l][k][i] * spatialdistav[l][k][i]);
      }
      
      // Computing diagonal of second block
      for (i = 0; i < p; ++i){
        counter = cumsumblocksizes[k];
        
        while (counter < cumsumblocksizes[k + 1]){
          block2diag[i] = block2diag[i] + 
            (spatialdistarray[counter][l][i] * spatialdistarray[counter][k][i]);
          counter = counter + 1;
        }
        block2diag[i] = block2diag[i] / (blocksizes[k] - 1);                       /////////////////
        
        block2diag[i] = block2diag[i] - (spatialdistav[k][l][i] * spatialdistav[k][k][i]);
        
        block2diag[i] = - block2diag[i];
      }
      
      // Computing diagonal of third block
      for (i = 0; i < p; ++i){
        counter = cumsumblocksizes[k];
        
        while (counter < cumsumblocksizes[k + 1]){
          block3diag[i] = block3diag[i] + 
            (spatialdistarray[counter][k][i] * spatialdistarray[counter][l][i]);
          counter = counter + 1;
        }
        block3diag[i] = block3diag[i] / (blocksizes[k] - 1);                       /////////////////
        
        block3diag[i] = block3diag[i] - (spatialdistav[k][k][i] * spatialdistav[k][l][i]);
        
        block3diag[i] = - block3diag[i];
      }
      
      // Adding diagonals of block1, block2, block3 over l
      for (i = 0; i < p; ++i){
        firsttermdiag[i] = firsttermdiag[i] + 
          ( lambda[l] * (block1diag[i] + block2diag[i] + block3diag[i]) );
      }
    }
    
    for (i = 0; i < p; ++i){
      firsttermdiag[i] = lambda[k] * firsttermdiag[i];
    }
    
    // Initialising and computing diagonal of second composite term in the expression of sigma(k, k)
    for (i = 0; i < p; ++i){
      secondtermdiag[i] = 0;
    }
    
    for (l1 = 0; l1 < numblocks; ++l1){
      for (l2 = 0; l2 < numblocks; ++l2){
        // Initialising all elements of diagonal of block4 to 0
        for (i = 0; i < p; ++i){
          block4diag[i] = 0;
        }
        
        // Computing diagonal of fourth block
        for (i = 0; i < p; ++i){
          counter = cumsumblocksizes[k];
          
          while (counter < cumsumblocksizes[k + 1]){
            block4diag[i] = block4diag[i] + 
              (spatialdistarray[counter][l1][i] * spatialdistarray[counter][l2][i]);
            counter = counter + 1;
          }
          block4diag[i] = block4diag[i] / (blocksizes[k] - 1);                  /////////////////
          
          block4diag[i] = block4diag[i] - (spatialdistav[k][l1][i] * spatialdistav[k][l2][i]);
        }
        
        // Adding diagonal of block4 over l1, l2
        for (i = 0; i < p; ++i){
          secondtermdiag[i] = secondtermdiag[i] + ( lambda[l1] * lambda[l2] * block4diag[i] );
        }
      }
    }
    
    // Addition to traceCov
    for (i = 0; i < p; ++i){
      traceCov = traceCov + (firsttermdiag[i] + secondtermdiag[i]);
    }
  }
  
  return traceCov;
}


// Computing SS statistic

// [[Rcpp::export]]

double SSstatistic(NumericMatrix Data, NumericVector blocksizes, double resolution) {
  int n = Data.nrow(), p = Data.ncol(), numblocks = blocksizes.size();
  int i, j, k, ik;
  int currentsize, cumsumblocksizes[(numblocks + 1)];
  double SSsummand = 0;
  double distance, normsquare;
  double vector[p], vectortemp[p], spatialdistvectorsum[p];
  
  cumsumblocksizes[0] = 0;
  for (k = 1; k <= numblocks; ++k){
    cumsumblocksizes[k] = cumsumblocksizes[k - 1] + blocksizes[k - 1];
  }
  
  for (k = 0; k < numblocks; ++k){
    // Initializing spatialdistvectorsum vector to 0
    for (j = 0; j < p; ++j){
      spatialdistvectorsum[j] = 0;
    }
    
    currentsize = blocksizes[k];
    
    for (ik = 0; ik < currentsize; ++ik){
      // Fixing current vector
      for (j = 0; j < p; ++j){
        vector[j] = Data(cumsumblocksizes[k] + ik, j);
      }
      
      for (i = 0; i < n; ++i){
        for (j = 0; j < p; ++j){
          vectortemp[j] = Data(i, j);
        }
        
        distance = 0;
        for (j = 0; j < p; ++j){
          distance = distance + (pow(vector[j] - vectortemp[j], 2) / resolution);
        }
        distance = sqrt(distance);
        
        if (distance > 1e-15){
          for (j = 0; j < p; ++j){
            spatialdistvectorsum[j] = spatialdistvectorsum[j] + 
              ((vector[j] - vectortemp[j]) / (distance * n));
          }
        }
      }
    }
    
    for (j = 0; j < p; ++j){
      spatialdistvectorsum[j] = spatialdistvectorsum[j] / currentsize;
    }
    
    normsquare = 0;
    for (j = 0; j < p; ++j){
      normsquare = normsquare + (pow(spatialdistvectorsum[j], 2) / resolution);
    }
    
    SSsummand = SSsummand + (normsquare * currentsize);
  }
  
  return SSsummand;
}


// Computing pvalue for SS test based on nonparametric bootstrap

// [[Rcpp::export]]

double SSbootstrappvalue(NumericMatrix Data, NumericVector blocknumbers, double resolution,
                         int bootstrapnum) {
  int n = Data.nrow(), p = Data.ncol(), numblocks = blocknumbers.size();
  int i, j, k, ik;
  int currentsize, cumsumblocknumbers[(numblocks + 1)];
  double SSsummand = 0;
  double distance, normsquare;
  double vector[p], vectortemp[p], spatialdistvectorsum[p];
  bool wrwor = TRUE;
  
  cumsumblocknumbers[0] = 0;
  for (k = 1; k <= numblocks; ++k){
    cumsumblocknumbers[k] = cumsumblocknumbers[k - 1] + blocknumbers[k - 1];
  }
  
  for (k = 0; k < numblocks; ++k){
    // Initializing spatialdistvectorsum vector to 0
    for (j = 0; j < p; ++j){
      spatialdistvectorsum[j] = 0;
    }
    
    currentsize = blocknumbers[k];
    
    for (ik = 0; ik < currentsize; ++ik){
      // Fixing current vector
      for (j = 0; j < p; ++j){
        vector[j] = Data(cumsumblocknumbers[k] + ik, j);
      }
      
      for (i = 0; i < n; ++i){
        for (j = 0; j < p; ++j){
          vectortemp[j] = Data(i, j);
        }
        
        distance = 0;
        for (j = 0; j < p; ++j){
          distance = distance + (pow(vector[j] - vectortemp[j], 2) / resolution);
        }
        distance = sqrt(distance);
        
        if (distance > 1e-15){
          for (j = 0; j < p; ++j){
            spatialdistvectorsum[j] = spatialdistvectorsum[j] + 
              ((vector[j] - vectortemp[j]) / (distance * n));
          }
        }
      }
    }
    
    for (j = 0; j < p; ++j){
      spatialdistvectorsum[j] = spatialdistvectorsum[j] / currentsize;
    }
    
    normsquare = 0;
    for (j = 0; j < p; ++j){
      normsquare = normsquare + (pow(spatialdistvectorsum[j], 2) / resolution);
    }
    
    SSsummand = SSsummand + (normsquare * currentsize);
  }
  
  // Generating bootstrap samples and computing Rkbarstarmatrix
  int bootstrapindex;
  double SSbootstrapsummand;
  NumericMatrix Databootstrap(n, p);
  NumericVector SSbootstrap(bootstrapnum);
  
  for (bootstrapindex = 0; bootstrapindex < bootstrapnum; ++bootstrapindex){
    IntegerVector bootstrapsampleindices = sample(n, n, wrwor) - 1;
    for (i = 0; i < n; ++i){
      for (j = 0; j < p; ++j){
        Databootstrap(i, j) = Data(bootstrapsampleindices[i], j);
      }
    }
    
    SSbootstrapsummand = 0;
    
    // Computing Rkbarstarmatrix for Databootstrap
    for (k = 0; k < numblocks; ++k){
      // Initializing spatialdistvectorsum vector to 0
      for (j = 0; j < p; ++j){
        spatialdistvectorsum[j] = 0;
      }
      
      currentsize = blocknumbers[k];
      
      for (ik = 0; ik < currentsize; ++ik){
        // Fixing current vector
        for (j = 0; j < p; ++j){
          vector[j] = Databootstrap(cumsumblocknumbers[k] + ik, j);
        }
        
        for (i = 0; i < n; ++i){
          for (j = 0; j < p; ++j){
            vectortemp[j] = Databootstrap(i, j);
          }
          
          distance = 0;
          for (j = 0; j < p; ++j){
            distance = distance + (pow(vector[j] - vectortemp[j], 2) / resolution);
          }
          distance = sqrt(distance);
          
          if (distance > 1e-15){
            for (j = 0; j < p; ++j){
              spatialdistvectorsum[j] = spatialdistvectorsum[j] + 
                ((vector[j] - vectortemp[j]) / (distance * n));
            }
          }
        }
      }
      
      for (j = 0; j < p; ++j){
        spatialdistvectorsum[j] = spatialdistvectorsum[j] / currentsize;
      }
      
      normsquare = 0;
      for (j = 0; j < p; ++j){
        normsquare = normsquare + (pow(spatialdistvectorsum[j], 2) / resolution);
      }
      
      SSbootstrapsummand = SSbootstrapsummand + (normsquare * currentsize);
    }
    
    SSbootstrap[bootstrapindex] = SSbootstrapsummand;
  }
  
  int counter = 0;
  for (bootstrapindex = 0; bootstrapindex < bootstrapnum; ++bootstrapindex){
    if (SSbootstrap[bootstrapindex] >= SSsummand){
      counter = counter + 1;
    }
  }
  double pvalue = (double)counter / (double)bootstrapnum;
  
  return pvalue;
}


// Computing pvalue for SS test based on random permutations

// [[Rcpp::export]]

double SSpermpvalue(NumericMatrix Data, NumericVector blocknumbers, double resolution,
                    int permnum) {
  int n = Data.nrow(), p = Data.ncol(), numblocks = blocknumbers.size();
  int i, j, k, ik;
  int currentsize, cumsumblocknumbers[(numblocks + 1)];
  double SSsummand = 0;
  double distance, normsquare;
  double vector[p], vectortemp[p], spatialdistvectorsum[p];
  bool wrwor = FALSE;
  
  cumsumblocknumbers[0] = 0;
  for (k = 1; k <= numblocks; ++k){
    cumsumblocknumbers[k] = cumsumblocknumbers[k - 1] + blocknumbers[k - 1];
  }
  
  for (k = 0; k < numblocks; ++k){
    // Initializing spatialdistvectorsum vector to 0
    for (j = 0; j < p; ++j){
      spatialdistvectorsum[j] = 0;
    }
    
    currentsize = blocknumbers[k];
    
    for (ik = 0; ik < currentsize; ++ik){
      // Fixing current vector
      for (j = 0; j < p; ++j){
        vector[j] = Data(cumsumblocknumbers[k] + ik, j);
      }
      
      for (i = 0; i < n; ++i){
        for (j = 0; j < p; ++j){
          vectortemp[j] = Data(i, j);
        }
        
        distance = 0;
        for (j = 0; j < p; ++j){
          distance = distance + (pow(vector[j] - vectortemp[j], 2) / resolution);
        }
        distance = sqrt(distance);
        
        if (distance > 1e-15){
          for (j = 0; j < p; ++j){
            spatialdistvectorsum[j] = spatialdistvectorsum[j] + 
              ((vector[j] - vectortemp[j]) / (distance * n));
          }
        }
      }
    }
    
    for (j = 0; j < p; ++j){
      spatialdistvectorsum[j] = spatialdistvectorsum[j] / currentsize;
    }
    
    normsquare = 0;
    for (j = 0; j < p; ++j){
      normsquare = normsquare + (pow(spatialdistvectorsum[j], 2) / resolution);
    }
    
    SSsummand = SSsummand + (normsquare * currentsize);
  }
  
  // Generating permutation samples and computing Rkbarstarmatrix
  int permindex;
  double SSpermsummand;
  NumericMatrix Dataperm(n, p);
  NumericVector SSperm(permnum);
  
  for (permindex = 0; permindex < permnum; ++permindex){
    IntegerVector permsampleindices = sample(n, n, wrwor) - 1;
    for (i = 0; i < n; ++i){
      for (j = 0; j < p; ++j){
        Dataperm(i, j) = Data(permsampleindices[i], j);
      }
    }
    
    SSpermsummand = 0;
    
    // Computing Rkbarstarmatrix for Dataperm
    for (k = 0; k < numblocks; ++k){
      // Initializing spatialdistvectorsum vector to 0
      for (j = 0; j < p; ++j){
        spatialdistvectorsum[j] = 0;
      }
      
      currentsize = blocknumbers[k];
      
      for (ik = 0; ik < currentsize; ++ik){
        // Fixing current vector
        for (j = 0; j < p; ++j){
          vector[j] = Dataperm(cumsumblocknumbers[k] + ik, j);
        }
        
        for (i = 0; i < n; ++i){
          for (j = 0; j < p; ++j){
            vectortemp[j] = Dataperm(i, j);
          }
          
          distance = 0;
          for (j = 0; j < p; ++j){
            distance = distance + (pow(vector[j] - vectortemp[j], 2) / resolution);
          }
          distance = sqrt(distance);
          
          if (distance > 1e-15){
            for (j = 0; j < p; ++j){
              spatialdistvectorsum[j] = spatialdistvectorsum[j] + 
                ((vector[j] - vectortemp[j]) / (distance * n));
            }
          }
        }
      }
      
      for (j = 0; j < p; ++j){
        spatialdistvectorsum[j] = spatialdistvectorsum[j] / currentsize;
      }
      
      normsquare = 0;
      for (j = 0; j < p; ++j){
        normsquare = normsquare + (pow(spatialdistvectorsum[j], 2) / resolution);
      }
      
      SSpermsummand = SSpermsummand + (normsquare * currentsize);
    }
    
    SSperm[permindex] = SSpermsummand;
  }
  
  int counter = 0;
  for (permindex = 0; permindex < permnum; ++permindex){
    if (SSperm[permindex] >= SSsummand){
      counter = counter + 1;
    }
  }
  double pvalue = (double)counter / (double)permnum;
  
  return pvalue;
}


// Computing pvalue for SS test based on random permutations (studentized)

// [[Rcpp::export]]

double SSpermstudentizedpvalue(NumericMatrix Data, NumericVector blocknumbers, double resolution,
                    int permnum) {
  int n = Data.nrow(), p = Data.ncol(), numblocks = blocknumbers.size();
  int i, j;
  bool wrwor = FALSE;
  
  double SSstatisticvalue = SSstatistic(Data, blocknumbers, resolution);
  
  NumericMatrix SSCov = SSCovOperatoradjusted(Data, blocknumbers, resolution);
  
  double SSCovtrace = 0;
  for (i = 0; i < SSCov.nrow(); ++i){
    SSCovtrace = SSCovtrace + SSCov(i, i);
  }
  
  double Statistic;
  if (SSCovtrace == 0){
    Statistic = 0;
  }else{
    Statistic = SSstatisticvalue / SSCovtrace;
  }
  
  // Generating permutation samples and computing Rkbarstarmatrix
  int permindex;
  double SSpermstatistic, SSpermCovtrace;
  NumericMatrix SSpermCov((numblocks * p), (numblocks * p));
  NumericMatrix Dataperm(n, p);
  NumericVector Statisticperm(permnum);
  
  for (permindex = 0; permindex < permnum; ++permindex){
    IntegerVector permsampleindices = sample(n, n, wrwor) - 1;
    for (i = 0; i < n; ++i){
      for (j = 0; j < p; ++j){
        Dataperm(i, j) = Data(permsampleindices[i], j);
      }
    }
    
    SSpermstatistic = SSstatistic(Dataperm, blocknumbers, resolution);
    
    SSpermCov = SSCovOperatoradjusted(Dataperm, blocknumbers, resolution);
    
    SSpermCovtrace = 0;
    for (i = 0; i < SSpermCov.nrow(); ++i){
      SSpermCovtrace = SSpermCovtrace + SSpermCov(i, i);
    }
    
    if (SSpermCovtrace == 0){
      Statisticperm[permindex] = 0;
    }else{
      Statisticperm[permindex] = SSpermstatistic / SSpermCovtrace;
    }
  }
  
  int counter = 0;
  for (permindex = 0; permindex < permnum; ++permindex){
    if (Statisticperm[permindex] >= Statistic){
      counter = counter + 1;
    }
  }
  double pvalue = (double)counter / (double)permnum;
  
  return pvalue;
}


// Computing pvalue for SS test based on random permutations (studentized and optimized)

// [[Rcpp::export]]

double SSpermstudentizedfastpvalue(NumericMatrix Data, NumericVector blocknumbers, double resolution,
                               int permnum) {
  int n = Data.nrow(), p = Data.ncol();
  int i, j;
  bool wrwor = FALSE;
  
  double SSstatisticvalue = SSstatistic(Data, blocknumbers, resolution);
  
  double SSCovtrace = SSfasttraceCovOperatoradjusted(Data, blocknumbers, resolution);
  
  double Statistic;
  if (SSCovtrace == 0){
    Statistic = 0;
  }else{
    Statistic = SSstatisticvalue / SSCovtrace;
  }
  
  // Generating permutation samples and computing Rkbarstarmatrix
  int permindex;
  double SSpermstatistic, SSpermCovtrace;
  NumericMatrix Dataperm(n, p);
  NumericVector Statisticperm(permnum);
  
  for (permindex = 0; permindex < permnum; ++permindex){
    IntegerVector permsampleindices = sample(n, n, wrwor) - 1;
    for (i = 0; i < n; ++i){
      for (j = 0; j < p; ++j){
        Dataperm(i, j) = Data(permsampleindices[i], j);
      }
    }
    
    SSpermstatistic = SSstatistic(Dataperm, blocknumbers, resolution);
    
    SSpermCovtrace = SSfasttraceCovOperatoradjusted(Dataperm, blocknumbers, resolution);
    
    if (SSpermCovtrace == 0){
      Statisticperm[permindex] = 0;
    }else{
      Statisticperm[permindex] = SSpermstatistic / SSpermCovtrace;
    }
  }
  
  int counter = 0;
  for (permindex = 0; permindex < permnum; ++permindex){
    if (Statisticperm[permindex] >= Statistic){
      counter = counter + 1;
    }
  }
  double pvalue = (double)counter / (double)permnum;
  
  return pvalue;
}


