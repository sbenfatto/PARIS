% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RF_feature_selection.R
\name{RF_exp}
\alias{RF_exp}
\title{Random Forest Feature Selection - Genes expression -> Dependent gene}
\usage{
RF_exp(data, exp, imp.score = c("raw", "Gini", "Ginicorrected"), outdir)
}
\arguments{
\item{data}{Dataframe with dependecy scores}

\item{exp}{Dataframe with expression data}

\item{imp.score}{importance score methods; NB raw/Gini for Boruta; Ginicorrected RF}

\item{outdir}{Directory where to write Models}
}
\value{
Write Boruta/Random Forest Feature selection models and selected gene paris in the outdir
}
\description{
Random Forest Feature Selection - Genes expression -> Dependent gene
}
