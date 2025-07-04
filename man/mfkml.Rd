% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfkml.R
\name{mfkml}
\alias{mfkml}
\title{Multidimensional Fréchet Distance-Based K-means for Longitudinal Data}
\usage{
mfkml(dt, clt_n, scales, weight, maxIter = 50)
}
\arguments{
\item{dt}{A long-format data.frame containing the following columns in the specified order:
\itemize{
\item \code{ID}: An identifier for each trajectory.
\item \code{Time}: The time points at which measurements were recorded (numeric or integer vector).
\item \code{Variable1}, \code{Variable2}, ... : The measured variables over time (numeric values).
The data.frame should not include any missing values.
See 'Details' for structure requirements.
}}

\item{clt_n}{An integer specifying the number of clusters.
The number of unique trajectories must be greater than or equal to \code{clt_n}.}

\item{scales}{A numeric vector used for scaling the time and variable columns. The length of \code{scales} must be equal to \code{ncol(dt) - 1},
where each value in \code{scales} corresponds to the scaling factor for the respective column (excluding the ID column).
See 'Details' for structure requirements.}

\item{weight}{Specifies the weights used for calculating the weighted Fréchet mean. It can take one of the following forms:
\itemize{
\item A data.frame with two columns: \code{ID} and \code{Weight}, where each \code{Weight} value indicates the importance of the corresponding trajectory.
\item A numeric value of 1, indicating equal weights for all trajectories.
See 'Details' for structure requirements.
}}

\item{maxIter}{The maximum number of iterations allowed before stopping if convergence is not reached. The default value is 50.}
}
\value{
A list with the following components:
\describe{
\item{\code{Cluster}}{A data.frame containing the \code{ID} and \code{Cluster} columns, which indicate the final cluster assignment for each trajectory.}
\item{\code{Center}}{A data.frame representing the final cluster centers, with columns for the cluster IDs, time points, and variable values.}
\item{\code{Iteration}}{The number of iterations the algorithm performed before reaching convergence.}
}
}
\description{
Extends \code{kmlShape} to multidimensional (p ≥ 2) longitudinal data.
It performs scale adjustment and trajectory alignment across all variables prior to clustering
to reduce distortions caused by differences in time grids and amplitude scales.
When variables exhibit substantially different ranges, standardization is required to prevent any single variable
from disproportionately influencing the clustering outcome.

The clustering process follows an iterative K-means framework, where cluster assignments are updated based on Fréchet distances.
Cluster centers are computed using the weighted Fréchet mean, which accounts for variable weights assigned to individual trajectories.
This allows the mean to be adjusted according to the relative importance of each trajectory in the clustering process.
}
\details{
The input dataset (\code{dt}) must contain only numeric values (except for the ID column)
and must not include any missing values.
Each variable should be measured at least three times per trajectory,
since the method relies on trajectory shapes.
Two observations per trajectory are insufficient to capture shape trends (e.g., increasing, decreasing, or stable).

Because the Fréchet distance is sensitive to measurement units, proper scaling is essential when applying the \code{mfkml} function.
The \code{scales} vector contains scaling factors for time and each variable,
which are used to rescale the corresponding columns.
This scaling prevents distortion due to differences in the units of time and variables,
allowing for more accurate shape-based comparisons.

This function involves random sampling internally.
For reproducible results, set the random seed before calling the function using \code{set.seed()}.
}
