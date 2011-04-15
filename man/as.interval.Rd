\name{as.interval}
\alias{as.interval}
\title{
Specify a continous interval
}
\description{
Specifies a continous interval by supplying a lower and upper bound. Such intervals can be supplied as parameter value ranges in \code{\link{tunePareto}}.
}
\usage{
as.interval(lower, upper)
}
\arguments{
  \item{lower}{
The lower bound of the interval
}
  \item{upper}{
The upper bound of the interval
}
}
\value{
A list of class \code{Interval} specifying the lower and upper bound.
}

\seealso{
\code{\link{tunePareto}}
}

\keyword{interval}
\keyword{continuous parameter}
