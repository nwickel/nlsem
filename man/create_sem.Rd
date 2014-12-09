% File src/library/base/man/create_sem.Rd
\name{create_sem}
\alias{create_sem}
\title{Create a structural equation model from a data frame}
\description{
    Create model matrices from a data frame with columns \code{label} (for
    parameter labels) and \code{group1} to \code{groupX}.
}
\usage{
    create_sem(dat)
}
\arguments{
    \item{dat}{Dataframe with first column \code{label} and one column for
    each group labeled \code{group1}, \code{group2} and so on. See
    Details.}
}
\details{
    Labels on column \code{label} need to be labeled in a certain way.
    Labels can be looked at when an object created with
    \code{\link{specify_sem}} is transformed to a data frame with
    \code{as.data.frame}. See examples below.
}
\value{
    Gives back an object of class 'lms', 'semm', or 'nsemm' which can be
    used to estimate parameters using \code{\link{em}}.
}
\seealso{\code{\link{specify_sem}}}
\examples{
# specify model
model <- specify_sem(4, 1, 2, 1, xi="x1-x2,x3-x4", eta="y1")
# create data frame
dat <- as.data.frame(model)
# recreate model
create_sem(dat)
}
