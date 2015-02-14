\name{Legend}
\alias{Legend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add a legend to plots
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Legend(Bin, Col, RowSet = NULL, Pcex = 3, Tcex = 3, Digits = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Bin}{
%%     ~~Describe \code{Bin} here~~
}
  \item{Col}{
%%     ~~Describe \code{Col} here~~
}
  \item{RowSet}{
%%     ~~Describe \code{RowSet} here~~
}
  \item{Pcex}{
%%     ~~Describe \code{Pcex} here~~
}
  \item{Tcex}{
%%     ~~Describe \code{Tcex} here~~
}
  \item{Digits}{
%%     ~~Describe \code{Digits} here~~
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
function (Bin, Col, RowSet = NULL, Pcex = 3, Tcex = 3, Digits = 3) 
{
    if (is.null(RowSet)) 
        RowSet = 1:Bin$Nregions
    plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "", 
        ylab = "", frame.plot = FALSE, xaxt = "n", yaxt = "n", 
        xaxs = "i", yaxs = "i", mar = c(0, 0, 0, 0))
    Y = seq(0, 1, length = length(RowSet) + 2)[-c(length(RowSet) + 
        1:2)] + 1/length(RowSet)/2
    points(x = rep(0.1, length(RowSet)), y = Y, col = Col(Bin$Nregions)[RowSet], 
        pch = 20, cex = Pcex)
    text(x = rep(0.2, length(RowSet)), y = Y, labels = paste(formatC(Bin$Lwr[RowSet], 
        format = "f", digits = Digits), "to", formatC(Bin$Upr[RowSet], 
        format = "f", digits = Digits)), pos = 4, cex = Tcex)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
