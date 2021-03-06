% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random-forest.R
\name{rf}
\alias{rf}
\title{Random Survival Forest Wrapper Function}
\usage{
rf(data, seed)
}
\arguments{
\item{data}{Tibble of data containing the variables to be regressed.}

\item{seed}{The random seed to standardize random forest creation}
}
\value{
\describe{
\item{rand.tree}{The random forest grown in this algorithm}
\item{tune.tree}{The tuning forest and parameters derived}
\item{vars.tree}{{The variables used to grow this forest}}}
}
\description{
This function wraps the \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc()}} function. It first runs an
initial tuning algorithm, identifying the most appropriate parameters for
random forest model fit. Then, the model runs the more detailed random forest
tree generation step. The tree is then returned. The time variable must be
stored in a variable named time, and outcome in variable named outc. These
variable translations are done automatically in the functions accessed by
users.
}
\examples{

\dontrun{

tree <- rf(visit.data, seed=101010)

TRUE USAGE:

tree.list <- filter.rf(visit.data, time="fuptime", outc="hfdiag",
                       iter = 30, adjust = adjustment.variables)
}
}
