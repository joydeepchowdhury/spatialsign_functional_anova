require(Rcpp)
Rcpp.package.skeleton("SSinternal", cpp_files = "SSinternal.cpp")
install.packages("SSinternal", repos=NULL, type="source")