% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sound_DMR.R
\name{run_binomial}
\alias{run_binomial}
\title{Run Binomial Model}
\usage{
run_binomial(LM, i = int, formula, optimizer_func = "optimizer")
}
\arguments{
\item{LM}{(df) - the information to put in the model, usually this is
in a "long" format}

\item{i}{(int) - the row number}

\item{formula}{(formula) - the formula to use in the model}

\item{optimizer_func}{(str) - the optimizer function to use in the model}

\item{Output_Frame}{(df) - all the summary information}
}
\value{
Output_Frame (df) - the data frame containing updated summary information
}
\description{
This function runs the glmer function using a binomial model for a given
optimizer function
}
