% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random-forest.R
\name{valid.rf}
\alias{valid.rf}
\title{Validation Tree Creation}
\usage{
valid.rf(data, time, outc, valid.terms)
}
\arguments{
\item{data}{The tibble containing the Validation dataset}

\item{time}{The time to event variable for this analysis}

\item{outc}{The event indicator variable for this analysis}

\item{valid.terms}{The pre-calculated important terms needed to generate the
tree.}
}
\value{
\describe{ \item{tree}{This is the Random Forest grown using the
supplied terms} \item{vsel}{This is the calcuation of relative variable
importance, showing how the provided terms replicate in the validation
dataset.}}
}
\description{
This function generates a random forest from a set of derived important
variables. For example, if a set of important random forest variables were
calculated from the ARIC Visit 5 Data, this function then uses the HUNT data
to create a tree from those variables. This tree is then returned, evaluating
how these variables perform in informativeness in the new dataset.
}
\examples{

\dontrun{

hunt.data <- haven::read_dta('hunt_data.dta')
time.var <- 'fuptime'
outc.var <- 'hfdiag'

valid.terms <- readRDS('ARIC_random_forest_derived_terms.RDS')

validation.tree <- valid.rf(data = hunt.data, time = time.var,
                            outc = outc.var, valid.terms = valid.terms)

}

}
