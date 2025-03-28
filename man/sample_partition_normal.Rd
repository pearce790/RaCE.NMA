% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_partition_normal.R
\name{sample_partition_normal}
\alias{sample_partition_normal}
\title{Partition sampling in Bayesian RaCE-NMA models (internal use only)}
\usage{
sample_partition_normal(
  mu_hat,
  J,
  nu,
  g,
  K,
  mu0,
  sigma0,
  sigma_hat,
  tau = tau,
  b_g = 0.5,
  d_g = 0.5
)
}
\arguments{
\item{mu_hat}{A vector of estimated average relative intervention effects based on a previous NMA. The jth entry is the effect of intervention j.}

\item{J}{A numeric indicating the total number of interventions being compared.}

\item{nu}{A vector indicating current values for nu in the Gibbs sampler.}

\item{g}{A vector indicating current values for g in the Gibbs sampler.}

\item{K}{A vector indicating current values for K in the Gibbs sampler.}

\item{mu0}{The hyperparameter mu0, usually specified as the grand mean of the average intervention effects.}

\item{sigma0}{The hyperparameter sigma_0, usually a large number as to be minimally informative.}

\item{sigma_hat}{A vector of the estimated standard deviations of each intervention. The jth entry is the standard deviation of intervention j.}

\item{tau}{The standard deviation of the Metropolis Hastings proposal distribution.}

\item{b_g}{The probability of "birth"ing a new partition cluster, if possible. Default is 0.5.}

\item{d_g}{The probability of "death"ing an existing partition cluster, if possible. Default is 0.5.}
}
\value{
A list containing updated values for g, nu, and K.
}
\description{
This function implements a reversible jump MCMC procedure for updating the parameter partition in Bayesian Rank-Clustered Estimation for Network Meta-Analysis models. For internal use only.
}
