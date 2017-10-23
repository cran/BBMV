\name{trans_from_fixed}
\alias{trans_from_fixed}
\alias{trans_to_fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear transformations
}
\description{
Internal function used to convert back and forth between the actual trait interval and [-1.5;1.5]}
\usage{
trans_from_fixed(x, bounds)
trans_to_fixed(x, bounds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A single value or vector of trait values.
}
  \item{bounds}{
The actual bounds on the trait interval
}
}

\value{
A single value or vector of trait values transformed to the other interval.
}
\author{
F. C. Boucher
}