\name{ggsave.latex}
\alias{ggsave.latex}
\title{Saves a ggplot graphic to a file and creates the code to include it in a
LaTeX document.}
\usage{
  ggsave.latex(..., caption = NULL, label = NULL,
    figure.placement = "hbt", floating = TRUE,
    caption.placement = "bottom",
    latex.environments = "center")
}
\arguments{
  \item{...}{arguments passed to the ggsave function
  (\code{\link[ggplot2]{ggsave}})}

  \item{caption}{The caption. Default to NULL, indicating
  no caption.}

  \item{label}{The label. Default to NULL, indicating no
  label.}

  \item{figure.placement}{The placement of the figure.
  Default to "hbt".}

  \item{floating}{Logical. Indicates if the figure should
  be placed in a floating environment. Default to TRUE}

  \item{caption.placement}{Should the caption be on top or
  bottom of the figure. Default to "bottom"}

  \item{latex.environments}{Alignment of the figure.
  Default to "center".}
}
\value{
  The graphic will be saved to a plot and the relevant
  LaTeX code is printed.
}
\description{
  Saves a ggplot graphic to a file and creates the code to
  include it in a LaTeX document.
}
\examples{
require(ggplot2)
 data(cars)
	p <- ggplot(cars, aes(x = speed, y = dist)) + geom_point()
	ggsave.latex(p, filename = "test.pdf", label = "fig:Cars",
   caption = "A demo plot", height = 5, width = 4)
}
\author{
  Thierry Onkelinx \email{Thierry.Onkelinx@inbo.be}, Paul
  Quataert
}
\seealso{
  \code{\link[ggplot2]{ggsave}}
}
\keyword{graphs}
\keyword{hplot}

