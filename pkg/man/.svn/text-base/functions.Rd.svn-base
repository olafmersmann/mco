\name{functions}

\alias{belegundu}
\alias{belegundu.constr}
\alias{binh1}
\alias{binh2}
\alias{binh2.constr}
\alias{binh3}
\alias{deb3}
\alias{fonseca1}
\alias{fonseca2}
\alias{gianna}
\alias{hanne1}
\alias{hanne1.constr}
\alias{hanne2}
\alias{hanne2.constr}
\alias{hanne3}
\alias{hanne3.constr}
\alias{hanne4}
\alias{hanne4.constr}
\alias{hanne5}
\alias{hanne5.constr}
\alias{jimenez}
\alias{jimenez.constr}
\alias{vnt}
\alias{zdt1}
\alias{zdt2}
\alias{zdt3}

\title{MCO test problems}
\description{
  Collection of functions implementing various MCO test problems.
}
\usage{
belegundu(x)
belegundu.constr(x)
binh1(x)
binh2(x)
binh2.constr(x)
binh3(x)
deb3(x)
fonseca1(x)
fonseca2(x)
gianna(x)
hanne1(x)
hanne1.constr(x)
hanne2(x)
hanne2.constr(x)
hanne3(x)
hanne3.constr(x)
hanne4(x)
hanne4.constr(x)
hanne5(x)
hanne5.constr(x)
jimenez(x)
jimenez.constr(x)
vnt(x)
zdt1(x)
zdt2(x)
zdt3(x)
}
\arguments{
  \item{x}{Input vector}
}
\value{
}
\references{ ~put references to the literature/web site here ~ }
\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}
\examples{

nsga2(belegundu, 2, 2,
      constraints=belegundu.constr, cdim=2,
      lower.bounds=c(0, 0), upper.bounds=c(5, 3))

nsga2(binh1, 2, 2,
      lower.bounds=c(-5, -5), upper.bounds=c(10, 10))
nsga2(binh2, 2, 2,
      lower.bounds=c(0, 0), upper.bounds=c(5, 3),
      constraints=binh2.constr, cdim=2)

nsga2(binh3, 2, 3,
      lower.bounds=c(10e-6, 10e-6), upper.bounds=c(10e6, 10e6))

nsga2(deb3, 2, 2,
      lower.bounds=c(0, 0), upper.bounds=c(1, 1),
      generations=500)

nsga2(fonseca1, 2, 2,
      lower.bounds=c(-100, -100), upper.bounds=c(100, 100))

nsga2(fonseca2, 2, 2,
      lower.bounds=c(-4, -4), upper.bounds=c(4, 4))

nsga2(gianna, 1, 2,
      lower.bounds=5, upper.bounds=10)

nsga2(hanne1, 2, 2,
      lower.bounds=c(0, 0), upper.bounds=c(10, 10),
      constraints=hanne1.constr, cdim=1)

nsga2(hanne2, 2, 2,
      lower.bounds=c(0, 0), upper.bounds=c(10, 10),
      constraints=hanne2.constr, cdim=1)

nsga2(hanne3, 2, 2,
      lower.bounds=c(0, 0), upper.bounds=c(10, 10),
      constraints=hanne3.constr, cdim=1)

nsga2(hanne4, 2, 2,
      lower.bounds=c(0, 0), upper.bounds=c(10, 10),
      constraints=hanne4.constr, cdim=1)

nsga2(hanne5, 2, 2,
      lower.bounds=c(0, 0), upper.bounds=c(10, 10),
      constraints=hanne5.constr, cdim=1)

nsga2(jimenez, 2, 2,
      lower.bounds=c(0, 0), upper.bounds=c(100, 100),
      constraints=jimenez.constr, cdim=4)

nsga2(vnt, 2, 3,
      lower.bounds=rep(-3, 2), upper.bounds=rep(3, 2))

nsga2(zdt1, 30, 2,
      lower.bounds=rep(0, 30), upper.bounds=rep(1, 30))

nsga2(zdt2, 30, 2,
      lower.bounds=rep(0, 30), upper.bounds=rep(1, 30))

nsga2(zdt3, 30, 2,
      lower.bounds=rep(0, 30), upper.bounds=rep(1, 30))
}
