% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sound_DMR.R
\name{create_fixed_effects}
\alias{create_fixed_effects}
\title{Create Fixed Effects}
\usage{
create_fixed_effects(fixed = c("effect1", "effect2"))
}
\arguments{
\item{fixed}{(list of strings) - the fixed effects elements}
}
\value{
fixed_effects (str) - the fixed effects properly formatted
to be passed into \code{create_function()}
}
\description{
A function to create the string combining fixed effects that will be passed
into \code{create_function()}
}
\examples{
# create multiple independent fixed effects
create_fixed_effects(c('Group', 'Individual'))
# [1] "Group + Individual"

# Create fixed effects with interactions
create_fixed_effects(c('Group * Individual'))
# [1] "Group * Individual"
}
