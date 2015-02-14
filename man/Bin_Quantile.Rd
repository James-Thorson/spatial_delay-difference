\name{Bin_Quantile}
\alias{Bin_Quantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bin regions into quantiles
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Bin_Quantile(Obj, Nregions = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Obj}{
%%     ~~Describe \code{Obj} here~~
}
  \item{Nregions}{
%%     ~~Describe \code{Nregions} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Obj, Nregions = 4) 
{
    Nbreaks = Nregions - 1
    Breaks = quantile(Obj, prob = seq(0, 1, length = Nbreaks + 
        2))
    Region = sapply(Obj, FUN = function(Num) {
        sum(Num >= Breaks)
    })
    Region = ifelse(Region == (Nregions + 1), Nregions, Region)
    if (is.array(Obj)) 
        Region = array(Region, dim = dim(Obj))
    Return = list(Nregions = Nregions, Region = Region, Lwr = Breaks[-c(Nbreaks + 
        2)], Upr = Breaks[-1])
    return(Return)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
