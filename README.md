# FKmL: A package for Fréchet distance-based K-means and extensions for Longitudinal data

The **`FKmL`** package provides shape-based clustering algorithms for multidimensional longitudinal data, using the generalized Fréchet distance to compare trajectory shapes. Traditional clustering approaches often have difficulty accounting for time shifts or speed variations in trajectories, while shape-based clustering methods, like those implemented in this package, offer greater flexibility in handling such variation.

This package implements two main methods:

MFKmL (Multidimensional Fréchet distance-based K-means for Longitudinal data) extends the K-means algorithm to handle multidimensional longitudinal data by using the Fréchet distance. It includes internal functions for computing Fréchet distance-related metrics for multidimensional trajectories.

SFKmL (Sparse Fréchet distance-based K-medoids for Longitudinal data) is a K-medoids-based clustering algorithm that incorporates variable selection. This method uses the permutation-based gap statistic to implement parameter tuning (uch as the L1-bound and time scale) and to automatically select relevant variables that contribute to clustering.


## Installation

You can install the **`FKmL`** package directly from GitHub using the `devtools` package:

```r
# Install devtools first if needed
install.packages("devtools")

devtools::install_github("jhp115/FKmL")
```

## Usage
You can explore example code for each function in the `examples` folder of the package.
The folder contains scripts that demonstrate how to use the MFKmL and SFKmL methods, including how to prepare data, apply the clustering algorithms, and interpret the results. Simply click on the files in the `examples` folder to see the examples in action.

### References
Kang et al. (2023). "Fréchet distance-based cluster analysis for multi-dimensional functional data." *Statistics and Computing*.
