% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/001functions.R
\name{grad_hess_eval}
\alias{grad_hess_eval}
\title{\code{grad_hess_eval} forms function of gradient and Hessian of log-likelihood produced 
by \code{f_create}.}
\usage{
grad_hess_eval(mn, parnl, hessian = FALSE, fixed = 0, data = "", cheqs0 = NULL)
}
\arguments{
\item{mn}{Expression, can be a list of equations.}

\item{parnl}{Names of parameters.}

\item{hessian}{if \code{TRUE}, returns hessian function, otherwise gradient.}

\item{fixed}{Integer, which parameter is fixed to be 0.}

\item{data}{Name of the data frame with which the function will be evaluated.}

\item{cheqs0}{If continuous are supplied, include the expressions of errors.}
}
\value{
A function for evaluation of gradient or Hessian.
}
\description{
\code{grad_hess_eval} forms function of gradient and Hessian of log-likelihood produced 
by \code{f_create}.
}
\examples{
eq_d <- c("ASC1 * 1 + B11_dur * dur_1" , "ASC2 * 1 + B12_dur * dur_2",
"ASC3 * 1 + B13_dur * dur_3 + B20_cost * cost_3 + B53_parkman * PbAvl_3",
"ASC4 * 1 + B14_dur * dur_4 + B20_cost * cost_4 + B34_serv * servIdx_4 + B44_stop * stopUs1R1_4")
parl <- c(paste0("ASC", 1:4), paste0("B1", 1:4, "_dur"), "B20_cost", "B53_parkman", "B34_serv",
 "B44_stop")
disc_par <- get_par(parl, eq_d)
ffor <- disc_par$cheqs0
parld <- disc_par$parld
res <- MNlogitf(ffor, separatenmm=FALSE, transform=FALSE)
parnl <- paste0("par", 1:length(parld))
gf <- grad_hess_eval (res, parnl, data="data", fixed=1)
hf <- grad_hess_eval (res, parnl, data="data", fixed=1, hessian=TRUE)
}
