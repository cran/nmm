#' \code{maxle} returns expression of log-likelihood (LL) of joint normal distribution.
#' @importFrom stats as.formula deriv logLik na.omit nls pnorm qnorm rnorm var AIC BIC
#' @param cheqs0 Strings defining equations of errors. Systems of Nonlinear Regressions (SNR) variant.
#' @param fixed_term if \code{TRUE} fixed term -(k/2)*log(2*pi) (k number of equations) is included
#' @return List with LL expressions of joint normal distribution, first element is string with 
#' expression for derivative calculations, the second - string for evaluation.
#' @examples
#' # normal distribution
#' eq_c <- c("Tw ~ ((((PH) + (tw)) * (ta - Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) -(1 + (PH))) +
#' sqrt((((PH) + (tw)) * (ta - Tc + 2) + (1 +(tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 - 4 * (1 + (PH) +
#'  (tw)) *(-(PH) * (ta - Tc + 2) + (1 - (tw) * (ta - Tc + 2)) * (2/w -Ec/w))))/(2 * (1 + (PH) +
#'   (tw)))",
#' "Tf1 ~ (th1) * (ta - (((((PH) + (tw)) * (ta - Tc + 2) + (1 + (tw)) *(Ec/w - 2/w) - (1 + (PH))) +
#'  sqrt((((PH) + (tw)) * (ta -Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 - 4 *(1 + (PH) +
#'  (tw)) * (-(PH) * (ta - Tc + 2) + (1 - (tw) *(ta - Tc + 2)) * (2/w - Ec/w))))/(2 * (1 + (PH) +
#'   (tw)))) -Tc + 2) - 1",
#' "Ef1 ~ (ph1)/(PH) * (w * (((((PH) + (tw)) * (ta - Tc + 2) + (1 +(tw)) * (Ec/w - 2/w) -
#' (1 + (PH))) + sqrt((((PH) + (tw)) *(ta - Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 -4 *
#'  (1 + (PH) + (tw)) * (-(PH) * (ta - Tc + 2) + (1 - (tw) *(ta - Tc + 2)) * (2/w - Ec/w))))/(2 *
#'   (1 + (PH) + (tw)))) -Ec + 2) - 1")
#' parl <- c("tw","PH","th1","ph1")
#' para_cont <- get_par(parl, eq_c)
#' cheqs0 <- para_cont$cheqs0
#' res <- maxle(cheqs0=cheqs0)
#' @export
maxle <- function(cheqs0, fixed_term = TRUE)
{
  . <- NULL
  npar <- get_npar(cheqs0)
  neq <- length(cheqs0)
  n1 <- 1:neq
  if(neq>4){
    rez <- ff_generate4maxle(neq)
  }else{
    rez <- datmaxle[[neq]]
  }
  names(rez$expr) <- c('dets', 'invs', 'expt')
  detS <- rez$expr$dets
  invS <- rez$expr$invs
  expt <- rez$expr$expt # ratsimp(transpose(Y) . invs . Y)

  # replace eps with equation expressions for residuals
  exptr <- expt
  y <- n1 %>% paste0('eps[', ., ']')
  yr <- y %>% gsub('[[]', '[[]', .) %>% gsub('[]]$', '[]]', .)
  repdat <- data.frame(old = yr, new = paste0('(', unlist(cheqs0), ')'), stringsAsFactors = FALSE)
  yyr <- replace_par_wrap(repdat , exptr)

  # should constant "(-neqt/2)*log(2*pi)" added, this is fixed to optimization solution is not affected
  fixedp <- if (fixed_term == TRUE) {
    paste0('(-', neq, '/2)*log(2*pi)-(1/2)*log(', detS, ')')
  } else{
    paste0('-(1/2)*log(', detS, ')')
  }

  #for derivatives
  rez <- yyr %>% paste0('-1/2*(', ., ')') %>% paste0(fixedp, ' ', .)
  reze <- formula2string(rez)
  list(expr = reze, formula = rez)
}

#' \code{maxle_p} returns expression of partitioned log-likelihood.
#'  f(y1,y2,..,yn)=f(y1)f(y2|y1)f(y3|y2y1)...f(yn|y1..y(n-1))
#' @param cheqs0 Strings defining equations of errors.
#' @param fixed_term if \code{TRUE} fixed term -(k/2)*log(2*pi) (k number of equations) is included
#' @param version2 another formulation of log-likelihood
#' @return List. First element is expression of joint distribution for derivatives,
#' second for evaluation, third latex, fourth marginal distributions for each variable.
#' @examples
#' # joint normal distribution
#' eq_c <- c("Tw ~ ((((PH) + (tw)) * (ta - Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) -(1 + (PH))) +
#'  sqrt((((PH) + (tw)) * (ta - Tc + 2) + (1 +(tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 - 4 * (1 + (PH) +
#'   (tw)) *(-(PH) * (ta - Tc + 2) + (1 - (tw) * (ta - Tc + 2)) * (2/w -Ec/w))))/(2 * (1 + (PH) +
#'    (tw)))",
#' "Tf1 ~ (th1) * (ta - (((((PH) + (tw)) * (ta - Tc + 2) + (1 + (tw)) *(Ec/w - 2/w) - (1 + (PH))) +
#'  sqrt((((PH) + (tw)) * (ta -Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 - 4 *(1 + (PH) +
#'   (tw)) * (-(PH) * (ta - Tc + 2) + (1 - (tw) *(ta - Tc + 2)) * (2/w - Ec/w))))/(2 * (1 + (PH) +
#'    (tw)))) -Tc + 2) - 1",
#' "Ef1 ~ (ph1)/(PH) * (w * (((((PH) + (tw)) * (ta - Tc + 2) + (1 +(tw)) * (Ec/w - 2/w) - (1 +
#'  (PH))) + sqrt((((PH) + (tw)) *(ta - Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 -4 *
#'   (1 + (PH) + (tw)) * (-(PH) * (ta - Tc + 2) + (1 - (tw) *(ta - Tc + 2)) * (2/w - Ec/w))))/(2 *
#'   (1 + (PH) + (tw)))) -Ec + 2) - 1")
#' parl <- c("tw","PH","th1","ph1")
#' para_cont <- get_par(parl, eq_c)
#' cheqs0 <- para_cont$cheqs0
#' res <- maxle_p(cheqs0=cheqs0)
#' @export
maxle_p <- function(cheqs0, fixed_term = TRUE, version2=TRUE)
{
  neqt <- length(cheqs0)
  n1 <- 1:neqt
  # get conditional means and variances
  # normal partitioned distribution
  
  if(version2){
    ff_generate4maxle_p_all <- ff_generate4maxle_p_v2
    expr_ll_norm_all <- expr_ll_norm_v2
  }else{
    ff_generate4maxle_p_all <- ff_generate4maxle_p
    expr_ll_norm_all <- expr_ll_norm
  }
  
  if(neqt>4){
    outl <- c()
    for(j in n1){
      sdv <- rep(NA, j)
      mv <- rep(0, j)
      outl <- c(outl, list(cond_expr(j, sdv, mv)))
    }
    aa <- ff_generate4maxle_p_all(outl, neqt, fixed_term )
  }else{
    aa <- expr_ll_norm_all[[which(names(expr_ll_norm_all) == paste0("ll_", neqt, "_", fixed_term))]]
  }
  
  llphi <-  gsub('%(pi)', '\\1', aa$expr[-(neqt + 1)])# replace pi expression
  llphi <- paste0('wc_', n1, '*(', llphi, ')')# add weights for cont. eq
  
  #insert expressions of equations
  repdat <- data.frame(old = paste0('eps', n1, ''), new = paste0('(', unlist(cheqs0), ')'),
                       stringsAsFactors = FALSE)
  #formula for each  eq. with weights
  llphi <- replace_par_wrap(repdat, llphi)
  #joint LL of all  eq.
  rez <- paste0(llphi, collapse = '+')
  reze <- rez %>% formula2string
  list(exprs = reze, formula = rez, latex = aa$latex, formulap = llphi)
}


#' \code{mlsem} returns expression of log-likelihood for joint normal distribution, 
#' for maximum likelihood (ML), Simultaneous Equations Models (SEM) variant. 
#' @param cheqs0 Strings defining equations of errors.
#' @param fixed_term if \code{TRUE} fixed term -(k/2)*log(2*pi) (k number of equations) is included
#' @return List with LL expressions of joint normal distribution, first element is string with 
#' expression for derivative calculations, the second - string for evaluation.
#' @examples
#' # normal distribution
#' eq_c <- c("Tw ~ ((((PH) + (tw)) * (ta - Tc) + Ec/w * (1 + (tw)) + sqrt((Ec/w *(1 + (tw)) +
#' (ta - Tc) * ((PH) + (tw)))^2 - 4 * Ec/w * (ta -Tc) * (tw) * (1 + (PH) + (tw))))/(2 * (1 + (PH)
#'  + (tw))))",
#' "Tf1 ~ (th1) *(Tw - Tc)",
#' "Ef1 ~ (ph1)/(PH) * (w * Tw - Ec)")
#' parl <- c("tw","PH","th1","ph1")
#' para_cont <- get_par(parl, eq_c)
#' cheqs0 <- para_cont$cheqs0
#' res <- mlsem(cheqs0, fixed_term=FALSE)
#' @export
mlsem <- function(cheqs0, fixed_term = TRUE)
{
  . <- NULL
  
  npar <- get_npar(cheqs0)
  neq <- length(cheqs0)
  ht <- cheqs0 %>% formula2string
  #nnpar <- para_cont$parn%>%formula2string
  endog <- cheqs0 %>% gsub("-.*", "", .) %>% trimws
  Jt <- lapply(ht, function(x)
    x %>% parse(text = .) %>%
      deriv(., namevec = endog) %>%
      extract_attr_deriv(., 'grad'))

  if(neq > 4){
    rez <- ff_generate4mlsem(neq)
  }else{
    rez <- datmlsem[[neq]]
  }
  names(rez$expr) <- c('dets', 'invs', 'expt', 'detJ')
  detS <- rez$expr$dets
  detJ <- rez$expr$detJ
  invS <- rez$expr$invs
  expt <- rez$expr$expt

  n1 <- 1:neq
  exptr <- expt
  y <- n1 %>% paste0('eps[', ., ']')
  yr <- y %>% gsub('[[]', '[[]', .) %>% gsub('[]]$', '[]]', .)
  repdat <- data.frame(old = yr, new = paste0('(', unlist(cheqs0), ')'), stringsAsFactors = FALSE)
  yyr <- replace_par_wrap(repdat, exptr)

  #form Jacobian
  dd <- n1 %>% expand.grid(., .) %>% apply(., 1, function(x)
      paste0(x, collapse = ',')) %>%
    paste0('deriv', '[', . , ']') %>% sort
  ddr <- dd %>% gsub('[[]', '[[]', .) %>% gsub('[]]$', '[]]', .)
  repdat <- data.frame(old = ddr, new = paste0('(', unlist(Jt), ')'), stringsAsFactors = FALSE)
  detJr <- replace_par_wrap(repdat, detJ)


  fixedp <- if(fixed_term == TRUE){
    paste0('(-', length(cheqs0), '/2)*log(2*pi)-(1/2)*log(', detS, ')', '+log(', detJr, ')')
  }else{
    paste0('-(1/2)*log(', detS, ')', '+log(', detJr, ')')
  }

  #for derivatives
  rez <- yyr %>% paste0('-1/2*(', ., ')') %>% paste0(fixedp, ' ', .)
  reze <- formula2string(rez)
  list(expr = reze, formula = rez)
}




#' Internal helper function for MNlogitf
#' @keywords internal
help_MNlogitf <- function(exps, transform = FALSE, separatenmm = FALSE, 
                          neq = NULL, weight_disc = FALSE)
{
  . <- NULL
  #return LL expression and string
  exp1 <- exps %>% paste0('log(', ., ')')
  n1 <- 1:neq
  if(transform == FALSE){
    exps <- n1 %>% paste0('avl_', .) %>% paste0(., ' * ', exps)
    sume <- exps %>% paste0(., collapse = '+')
    sume <- paste0( paste0('chc_', n1), ' * ',paste0('avl_', n1), ' * ( ', paste0(exp1,
                                                                          '- log(', sume, ')'),')')
  }else{
    sume <- exp1%>%paste0('avl_', n1, ' * ',.)%>%paste0('chc_', n1, ' * ',.)
  }
  if(weight_disc==TRUE){
    sume <- n1 %>% paste0('wd_', ., '*(', sume, ')')
  }

  loglik <- if(separatenmm == FALSE){
    paste0(sume, collapse = '+')
  }else{
    sume
  }
  loglike <- loglik %>% gsub('\\[(\\d{1,10})\\]', '\\1', .)
  list(loglik = loglik, loglike = loglike)
}

#' Internal helper function for MNdogitf
#' @keywords internal
help_MNdogitf <- function(exps, transform = FALSE, separatenmm = FALSE, neq = NULL,
                          ppardogit=NULL, weight_disc = FALSE)
{
  . <- NULL
  
  #return LL expression and string
  #exp1 <- exps #%>% paste0('log(', ., ')')
  #exps <- n1 %>% paste0('avl_', .) %>% paste0(., ' * ', exps)
  sume <- n1 %>% paste0('avl_', .) %>% paste0(., ' * ', exps) %>% paste0(., collapse = '+')
  exp1 <- paste0(exps, " + ", ppardogit, " * (", sume,")")%>% paste0('log(', ., ')')
  n1 <- 1:neq
  if(transform == FALSE){
    
    pjis <- ppardogit %>%
      paste0(., collapse = ' + ') %>% paste0('(1 + ', ., ')')
    sume <- paste0( paste0('chc_', n1), ' * ',paste0('avl_', n1), ' * ( ', paste0(exp1,
                                                                                  '-log(', pjis,')',
                                                                                  '- log(', sume, ')'),')')
  }else{
    sume <- exp1%>%paste0('avl_', n1, ' * ',.)%>%paste0('chc_', n1, ' * ',.)
  }
  if(weight_disc==TRUE){
    sume <- n1 %>% paste0('wd_', ., '*(', sume, ')')
  }
  
  loglik <- if(separatenmm == FALSE){
    paste0(sume, collapse = '+')
  }else{
    sume
  }
  loglike <- loglik %>% gsub('\\[(\\d{1,10})\\]', '\\1', .)
  list(loglik = loglik, loglike = loglike)
}

#' \code{MNlogitf} or \code{MNdogitf} returns log-likelihood(LL) expression for discrete equations of "logit" or "dogit" model.
#' @param ffor Discrete choice equations.
#' @param transform if \code{TRUE}, quantile transformation (normal) is applied.
#' @param separatenmm if \code{TRUE}, equation specific LL is calculated
#' @param weight_disc if \code{TRUE}, equations will include equation specific weights.
#' @describeIn MNlogitf returns log-likelihood(LL) expression for discrete equations of "logit" model.
#' @return
#' \item{formula}{simplified LL for each equation or joint}
#' \item{probs}{probability expression unsimplified}
#' \item{expr}{LL string simplified }
#' \item{formulat}{unsimplified log(P), only if transform is TRUE}
#' \item{probt}{unsimplified P with quantile transformation qnorm(ifelse(P)), only if transform is TRUE}
#' \item{expreteso}{unsimplified pnorm((qnorm(P)-mean)/sd), only if transform is TRUE}
#' \item{probte}{qnorm(P)}
#' \item{tval}{If quantile transformation is applied}
#' \item{sume}{Denominator of logit probability}
#' @examples
#' eq_d <- c("ASC1 * 1 + B11_dur * dur_1" , "ASC2 * 1 + B12_dur * dur_2",
#' "ASC3 * 1 + B13_dur * dur_3 + B20_cost * cost_3 + B53_parkman * PbAvl_3",
#' "ASC4 * 1 + B14_dur * dur_4 + B20_cost * cost_4 + B34_serv * servIdx_4 + B44_stop * stopUs1R1_4")
#' parl <- c(paste0("ASC", 1:4), paste0("B1", 1:4, "_dur"), "B20_cost", "B53_parkman", "B34_serv",
#' "B44_stop")
#' disc_par <- get_par(parl, eq_d)
#' ffor <- disc_par$cheqs0
#' res_l <- MNlogitf(ffor, separatenmm=FALSE, transform=FALSE)
#' res_d <- MNdogitf(ffor, separatenmm=FALSE, transform=FALSE)
#' @export
MNlogitf <- function(ffor, transform = FALSE, separatenmm = FALSE, weight_disc = FALSE)
{
  . <- NULL
  exps <- ffor %>% paste0('exp(', ., ')')
  neq <- ffor %>% length
  n1 <- 1:neq
  probs <- n1 %>% paste0('avl_', .) %>% paste0(., ' * ', exps)
  sume <- probs %>% paste0(., collapse = '+')

  #expression for prob
  probs <- probs %>% paste0(., '/(', sume, ')')
  environment(help_MNlogitf) <- environment()
  #if quantile transformation has to be done for probabilities
  #if(transform==TRUE){
  #apply quantile transformation
  probt <- probs %>%
    paste0('qnorm(ifelse((', .,')==0, 1e-100, ifelse((', .,')==1,1-1e-16,(', .,'))))')
  #get cumulative distribution values for the transformed probabilities
  mean_sd_part <- if(transform==TRUE) n1%>%
    paste0('- meanqt[',.,'])/sqrt(covqt[',.,'])') else ")"
  expret <- paste0('pnorm((', probt , mean_sd_part, ')')
  #########w/o ifelse, for derivatives
  #probabilities
  probte <- probs %>% paste0('qnorm(', ., ')') #quantile expression of probability
  #transformed probabilities
  #difference from probt in the construction of expression, probt is a simpliefied version of probte
  expreteso <- paste0('pnorm((', probte , mean_sd_part, ')')
  #another LL expression with ifelse
  expret1 <- expret%>%paste0('ifelse(', ., ' == 0, 1e-100, ifelse(', .,
                             ' == 1, 1-1e-16, ', ., '))')
  expst <- help_MNlogitf(expret1, transform = TRUE, separatenmm = separatenmm, neq = neq, 
                         weight_disc = weight_disc)

  #LL expression without quantile transformation
  exps <- help_MNlogitf(exps, transform = FALSE, separatenmm = separatenmm, neq = neq,
                        weight_disc = weight_disc)
  exprs <- if(length(exps$loglike) > 1){
    paste0(exps$loglike, collapse = ' + ')
  }else{
    exps$loglike
  }
  #probability string
  #probse <- gsub('\\[(\\d{1,10})\\]', '\\1', probs)

  #exps simplified log(P)
  #exps$loglike - for each equation separatenmm LL string expression
  #exps$loglik - exps$loglike but as an object which can be evaluated
  #probs - probability expression unsimplified P
  #exprs - LL string simplified

  #expst unsimplified log(P)
  #expst$loglike - for each equation separatenmm LL string expression
  #expst$loglik - expst$loglike but as an object which can be evaluated
  #probt unsimplified P with quantile transformation qnorm(ifelse(P))
  #exprst - LL string unsimplified with ifelse

  #expret pnorm((qnorm(ifelse(P))-mean)/sd)
  #exprete LL unsimplified
  #expretes for each eq LL separatenmm unsimplified
  #expreteso P unsimplified pnorm((qnorm(P)-mean)/sd)
  #probte - qnorm(P)
  list(formula = exps$loglik, probs = probs, expr = exprs,
       formulat = expst$loglik, probt = probt, 
       expreteso = expreteso, probte = probte, tval = transform, sume = sume)
}

#' 
#' @describeIn MNlogitf returns log-likelihood(LL) expression for discrete equations of "dogit" model.
#' @param ppardogit "dogit" parameters only used in MNdogitf
#' @export
MNdogitf <- function(ffor, transform = FALSE, separatenmm = FALSE, 
                     weight_disc = FALSE, ppardogit = NULL)
{
  . <- NULL
  exps <- ffor %>% paste0('exp(', ., ')')
  neq <- ffor %>% length
  n1 <- 1:neq
  probs <- n1 %>% paste0('avl_', .) %>% paste0(., ' * ', exps)
  sume <- probs %>% paste0(., collapse = '+')
  pjis <- c(ppardogit) %>%
    paste0(., collapse = ' + ') %>% paste0('(1 + ', ., ')')

  #expression for prob dogit
  probst1 <- probs %>% paste0(., ' + ', ppardogit, ' * (', sume, ')')
  probst2 <-  paste0('', pjis,' * (', sume, ')')
  
  probs <- paste0('(', probst1, ')/(', probst2, ')')
  #probs <- probs %>% paste0(., '/(', sume, ')')
  environment(help_MNdogitf) <- environment()
  #if quantile transformation has to be done for probabilities
  #if(transform==TRUE){
  #apply quantile transformation
  probt <- probs %>%
    paste0('qnorm(ifelse((', .,')==0, 1e-100, ifelse((', .,')==1,1-1e-16,(', .,'))))')
  #get cumulative distribution values for the transformed probabilities
  mean_sd_part <- if(transform==TRUE) n1%>%
    paste0('- meanqt[',.,'])/sqrt(covqt[',.,'])') else ")"
  expret <- paste0('pnorm((', probt , mean_sd_part, ')')
  #########w/o ifelse, for derivatives
  #probabilities
  probte <- probs %>% paste0('qnorm(', ., ')') #quantile expression of probability
  #transformed probabilities
  #difference from probt in the construction of expression, probt is a simpliefied version of probte
  expreteso <- paste0('pnorm((', probte , mean_sd_part, ')')
  #another LL expression with ifelse
  expret1 <- expret%>%paste0('ifelse(', ., ' == 0, 1e-100, ifelse(', .,
                             ' == 1, 1-1e-16, ', ., '))')
  expst <- help_MNdogitf(expret1, transform = TRUE, 
                         separatenmm = separatenmm, neq = neq, ppardogit=ppardogit, 
                         weight_disc= weight_disc)
  
  #LL expression without quantile transformation
  exps <- help_MNdogitf(exps, transform = FALSE, 
                        separatenmm = separatenmm, neq = neq, ppardogit=ppardogit, 
                        weight_disc= weight_disc)
  exprs <- if(length(exps$loglike) > 1){
    paste0(exps$loglike, collapse = ' + ')
  }else{
    exps$loglike
  }
  #probability string
  #probse <- gsub('\\[(\\d{1,10})\\]', '\\1', probs)
  
  #exps simplified log(P)
  #exps$loglike - for each equation separatenmm LL string expression
  #exps$loglik - exps$loglike but as an object which can be evaluated
  #probs - probability expression unsimplified P
  #exprs - LL string simplified
  
  #expst unsimplified log(P)
  #expst$loglike - for each equation separatenmm LL string expression
  #expst$loglik - expst$loglike but as an object which can be evaluated
  #probt unsimplified P with quantile transformation qnorm(ifelse(P))
  #exprst - LL string unsimplified with ifelse
  
  #expret pnorm((qnorm(ifelse(P))-mean)/sd)
  #exprete LL unsimplified
  #expretes for each eq LL separatenmm unsimplified
  #expreteso P unsimplified pnorm((qnorm(P)-mean)/sd)
  #probte - qnorm(P)
  list(formula = exps$loglik, probs = probs, expr = exprs,
       formulat = expst$loglik, probt = probt, 
       expreteso = expreteso, probte = probte, tval = transform, sume = sume)
}

#' \code{f_create} creates functions for log-likelihood of different models.
#' @param mn Expression, can be a list of equations.
#' @param data Name of the data frame with which the function will be evaluated.
#' @param fixed Integer, which parameter is fixed to be 0.
#' @param cheqs0 If continuous are supplied, include the expressions of errors.
#' @param probt Expressions of un-simplified probabilities with quantile transformation(qnorm(ifelse(P))).
#' @param tformula unsimplified log(P)
#' @param hessian Adds lines to check the Hessian, hessian should be the name of hessian function.
#' @param transform if \code{TRUE}, adds lines to check conditional means
#' @param separatenmm if \code{TRUE}, separate log-likelihood for each equations is produced.
#' @param sume Expression of summed likelihoods.
#' @return Function.
#' @examples
#' eq_d <- c("ASC1 * 1 + B11_dur * dur_1" , "ASC2 * 1 + B12_dur * dur_2",
#' "ASC3 * 1 + B13_dur * dur_3 + B20_cost * cost_3 + B53_parkman * PbAvl_3",
#' "ASC4 * 1 + B14_dur * dur_4 + B20_cost * cost_4 + B34_serv * servIdx_4 + B44_stop * stopUs1R1_4")
#' parl <- c(paste0("ASC", 1:4), paste0("B1", 1:4, "_dur"), "B20_cost", "B53_parkman", "B34_serv",
#' "B44_stop")
#' obj <- get_par(parl, eq_d)
#' ffor <- obj$cheqs0
#' res <- MNlogitf(ffor, separatenmm=FALSE, transform=FALSE)
#' ff <- f_create(res$formula, data="data", fixed=1)
#' @export
f_create <- function(mn, data, fixed = 0, cheqs0 = NULL, separatenmm = FALSE, probt = NULL,
                     tformula = NULL, hessian = NULL, transform = TRUE, sume = NULL)
{
  . <- NULL
  ffunc <- function(par){}
  neq <- length(mn)
  n1 <- 1:neq
  llf1 <- f_help(mn, data)
  #if for quantile transformation
  if(!is.null(tformula)){
    llf1t <- f_help(tformula, data)
    llf1t <- c('  if(transform==TRUE){',
               gsub('^  ', '    ', llf1t),
               '  }else{', gsub('^  ', '    ', llf1), '  }')
    llf1 <- llf1t
  }
  #add errors from cont. equations
  t22 <- cheqs0 %>% paste0('"', ., '"') %>% paste0(., collapse = ", ") %>%
    paste0('  err_eq <- c(', ., ')')
  t2 <- paste0('  err <- sapply(err_eq, function(x) eval(parse(text=x),', data, '))')
  t2 <- if(is.null(cheqs0)) t2 else c(t22, t2)
  t3 <- '  sigma <-  cov(err, use = "pairwise.complete.obs")'
  #if expressions include correlations instead of covariances
  if(any(grepl('sigma\\[', mn))){
    if(any(grepl('rho\\[', mn))){
      t3 <- c(t3, '  rho <- cor(err, use = "pairwise.complete.obs")',
              '  sigma <- sqrt(diag(sigma))')
    }else{
      t3 <-  c(t3, '  sigma <- sqrt(diag(sigma))')
      # does this happen that only variances and not correlations are used??
    }
    
  }
  t4 <- '  return(res)'
  #if some parameters is set to be equal to 0, for discrete model
  if(fixed > 0){
    t0 <- if(fixed == 1)'  par <- c(0, par)' else paste0('  par <- c(par[1:', (fixed - 1),
                                                         '], 0,  par[', fixed, ':length(par)])')
  }else{
    t0 <- ''
  }

  if(any(grepl('sigma\\[', mn))){
    tt <- c(t0, t2, t3)
  }else{
    tt <- c(t0)
  }
  #if cond. moments need to be included
  if(!is.null(probt)){
    neq <- length(probt)
    t5 <- n1 %>% paste0('    meanqt', ., " <- mean(eval(parse(text='", probt[.], "'), ", data,
                         '), na.rm=FALSE)')
    t5 <- c(t5, n1 %>% paste0('meanqt', ., collapse = ', ') %>% paste0('    meanqt <- c(', ., ')'))
    t6 <- n1%>%paste0('    covqt', ., " <- var(eval(parse(text='", probt[.], "'), ", data,
                         '), na.rm=FALSE)')
    t6 <- c(t6, n1 %>%paste0('covqt', ., collapse = ', ') %>% paste0('    covqt <- c(', ., ')')  )
    tt0 <- '  if(transform == TRUE){'
    tt <-  c(tt[1], tt0, t5, t6, if (length(tt) > 1) tt[2:length(tt)], '  }')
  }
  #if transform is TRUE check if cond. means have defined values
  if(transform==TRUE){
    t5 <- paste0(
      'if(transform==TRUE) if(any(meanqt%in%c(NaN, Inf, -Inf))|ifelse(any(is.na(covqt)), TRUE, any(covqt < 0.001))) res <- if(separatenmm==TRUE) rep(-Inf, ',
      length(mn), ') else  -Inf')
  }else{
    t5 <- ''
  }
  # include check for hessian function
  if(!is.null(hessian)){
    t05 <- c(paste0('if(check_hess==TRUE){', ' hess <- try(qr.solve(-', hessian, '(par[',
                    if(fixed!=0) -fixed else "", '])), silent=TRUE)'),
             'if(any(class(hess)%in%"try-error")) hess <- matrix(-0.01,1,1)', '}else{',
             'hess <- matrix(0.01,1,1)}')
    t15 <- c('nrepp <- 1',
      paste0('nrepp <- if(methodopt%in%c("BHHH")) nrow(' ,data ,') else nrepp'),
             paste0('nrepp <- if(separatenmm==TRUE)',length(mn), ' else nrepp'),
    paste0('if(any(diag(hess)<0)) res <- rep(-Inf, nrepp)'))
    t5<-c(t5, t05, t15)
  }
  tffunc <-  c('{', tt, llf1, t5, t4, '}')
  if(!is.null(sume)){
    tffunc %<>% gsub(sume, "dalyba", ., fixed = TRUE)
    xs <- paste0("   dalyba <- eval(parse(text='",sume, "'), ", data,")")
    xs1 <- paste0("   dalyba <- eval(parse(text='replace(dalyba, dalyba==Inf, 99999999999)'), ", data,")")
    xs2 <- paste0("   dalyba <- eval(parse(text='replace(dalyba, dalyba==-Inf, -99999999999)'), ", data,")")
    tffunc <- c(tffunc[1:2], xs, xs1, xs2, tffunc[3:length(tffunc)])
  }else{
    #ts <- NULL
  }
  body(ffunc) <- parse(text = tffunc)
  ffunc
}

#' Internal helper function
#' @keywords internal
#' @param mn expression to be converted to a function
#' @param data name of data frame
f_help <- function(mn, data)
{
  . <- NULL
  neq <- length(mn)
  n1 <- 1:neq
  if(neq > 1){
    #for partial conditional distributions
    llf0 <- n1 %>% paste0('    res', . , " <- eval(parse(text='", mn, "'), ", data, ')')
    llf1 <- c('  if(methodopt%in%c("BHHH")){',
              n1 %>% paste0('    res', ., ' <- res', .),
              '  }else{',
              n1 %>% paste0('    res', ., ' <- sum(res', ., ', na.rm=FALSE)'),
              '  }')
    llf1 <- c(llf1, n1 %>% paste0('res', ., collapse = ', ') %>% paste0('  res <- cbind(', ., ')'))
    llf1 <- c(llf1, c('  if(separatenmm==FALSE){','    res <- rowSums(res, na.rm=FALSE)', '  }'))
    llf1 <- c(llf0, llf1)
  }else{
    llf0 <-  paste0("    res <- eval(parse(text='", mn , "'), ", data, ')')
    llf00 <- "res <- replace(res, res%in%-Inf, -999999999999)"
    llf1 <- 'res <- if(methodopt%in%c("BHHH")) res else sum(res, na.rm=FALSE)'
    llf1 <- c(llf0, llf00, llf1)#, '  if(separatenmm==FALSE){ }')
  }
  llf1
}

#' \code{grad_hess_eval} forms function of gradient and Hessian of log-likelihood produced 
#' by \code{f_create}. 
#' @param mn Expression, can be a list of equations.
#' @param parnl Names of parameters.
#' @param hessian if \code{TRUE}, returns hessian function, otherwise gradient.
#' @param fixed Integer, which parameter is fixed to be 0.
#' @param data Name of the data frame with which the function will be evaluated.
#' @param cheqs0 If continuous are supplied, include the expressions of errors.
#' @return A function for evaluation of gradient or Hessian.
#' @examples
#' eq_d <- c("ASC1 * 1 + B11_dur * dur_1" , "ASC2 * 1 + B12_dur * dur_2",
#' "ASC3 * 1 + B13_dur * dur_3 + B20_cost * cost_3 + B53_parkman * PbAvl_3",
#' "ASC4 * 1 + B14_dur * dur_4 + B20_cost * cost_4 + B34_serv * servIdx_4 + B44_stop * stopUs1R1_4")
#' parl <- c(paste0("ASC", 1:4), paste0("B1", 1:4, "_dur"), "B20_cost", "B53_parkman", "B34_serv",
#'  "B44_stop")
#' disc_par <- get_par(parl, eq_d)
#' ffor <- disc_par$cheqs0
#' parld <- disc_par$parld
#' res <- MNlogitf(ffor, separatenmm=FALSE, transform=FALSE)
#' parnl <- paste0("par", 1:length(parld))
#' gf <- grad_hess_eval (res, parnl, data="data", fixed=1)
#' hf <- grad_hess_eval (res, parnl, data="data", fixed=1, hessian=TRUE)
#' @export
grad_hess_eval <- function(mn, parnl,  hessian = FALSE, fixed = 0, data = '', cheqs0 = NULL)
{
  . <- NULL
  #get derivatives
  gradMN <- deriv(parse(text = mn$expr) , parnl, hessian = hessian, function.arg = TRUE)
  mnf <- deparse(gradMN)
  #if line ends with , combine with previous line
  mnfo <- mnf
  ind <- grep(', $', mnfo)
  while(length(ind) > 0){
    mnfo[ind[1]] <- paste0(mnfo[ind[1]], mnfo[ind[1] + 1])
    mnfo <- mnfo[-(ind[1] + 1)]
    ind <- grep(', $', mnfo)
  }
  mnf <- mnfo
  llf1 <- string2formula(mnf)
  llf1 %<>% gsub('expression\\(', '', .)
  #remove first line
  llf1 <- llf1[!grepl('^function|(par\\[\\d{1,10}\\], par\\[\\d{1,10}\\])', llf1)]
  llf1 %<>% gsub('(})(\\))$', '\\1', .)
  ffunc <- function(par){}
  textv <- paste0('llf1 <- c(', paste0("'", llf1[2:(length(llf1) - 1)], "'", collapse = ', '), ')')
  #insert parameter vector
  inspar <- if(fixed == 1)'    par <- c(0, par)' else if(fixed > 1)  paste0('    par <- c(par[1:',
                                          (fixed -1), '], 0,  par[', fixed,
                                          ':length(par)])') else ''

  #add errors from cont. eq
  t22 <- cheqs0 %>% paste0('"', ., '"') %>% paste0(., collapse = ", ") %>%
    paste0('  err_eq <- c(', ., ')')
  t2 <- paste0('  err <- sapply(err_eq, function(x) eval(parse(text=x),', data, '))')
  t2 <- if(is.null(cheqs0)) t2 else c(t22, t2)

  t3 <- '  sigma <-  cov(err, use = "pairwise.complete.obs")'

  if(any(grepl('sigma\\[', textv))){
    if(any(grepl('rho\\[', textv))){
      t3 <- c(t3, '  rho <- cor(err, use = "pairwise.complete.obs")',
              '  sigma <- sqrt(diag(sigma))')
    }else{
      t3 <-  c(t3, '  sigma <- sqrt(diag(sigma))')
      # does this happen that only covariances and not correlations are used??
    }
  }  
  #if correlation(rho) is used instead of covariance
  if(any(grepl('sigma\\[', mn$formula))){
    tt <- c(t2, t3)
  }else{
    tt <- ''
  }

  llf1 <- c('{',tt, inspar, textv,  paste0('    res <- eval(parse(text=llf1), ',data, ')'))
  #add summation of gradient and hessian
  if(hessian){
    llf1 <- c(llf1, '    res <- attributes(res)$hessian',
              '    res <- if(methodopt%in%c("BHHH"))res else apply(res, 3, function(x)colSums(x, na.rm=FALSE))',
              if(fixed>0)paste0('    if(methodopt%in%c("BHHH")) return(res[,-', fixed, ',-', fixed,
                                ']) else return(res[-', fixed, ',-', fixed, '])') else '    return(res)')
  }else{
    llf1 <- c(llf1, '    res <- attributes(res)$gradient',
              '    res <- if(methodopt%in%c("BHHH"))res else apply(res, 2, function(x)sum(x, na.rm=FALSE))',
              if(fixed > 0) paste0('  if(methodopt%in%c("BHHH"))  return(res[, -',fixed,
                                 ']) else return(res[-', fixed, '])') else '    return(res)')
  }
  llf1 <- c(llf1, '}')
  llf1 <- parse(text = llf1)
  body(ffunc) <- llf1
  ffunc
}


#' \code{add_variable} adds columns to the data matrix
#' @param data data.frame, if \code{dname=="chc"} columns "chc_i" has to be in the data.
#' @param dname if \code{dname=="chc"} (dummy for chosen alternative) dummy for the choice 
#' alternative added, if "weights" weights added
#' @param weights Matrix with weights to be added to the data
#' @return data.frame
#' @examples
#' chc <- c(1,2,1,4,3,1,4)
#' data <- data.frame(choice=chc, x=rnorm(length(chc)), y=rnorm(length(chc)))
#' add_variable(data, dname="chc")
#' ww <- c(1,1,1,2,2,2,3)
#' add_variable(data, dname="weights", weights=ww)
#' @export
add_variable <- function(data = data, dname = "chc", weights = NULL)
{
  . <- NULL
  #check if all arguments supplied
  mf <- match.call()
  varl <- paste0(match.call()[[1]], '%>%args%>%as.list') %>% parse(text = .) %>% eval
  no_def <- varl[varl == ""] %>% names#no default values
  m <- match(no_def, mf %>% names, 0L)
  if(any(m==0)){
    paste0(no_def[m == 0], collapse = ", ") %>% paste0("No value for: ", .) %>%
      stop
  }

  if(!any(grepl(paste0(dname, "_"), names(data)))){
    if(dname=="chc"){
      misdata <- matrix(0, ncol = length(unique(data[, "choice"])), nrow = nrow(data))
      misdata <- data.frame(choice = data[, "choice"], misdata, stringsAsFactors = FALSE)
      names(misdata)[-1] <- paste0(paste0(dname, "_"), unique(data[, "choice"]))
      for(i in unique(data[,"choice"])){
        misdata[misdata$choice == i, paste0(paste0(dname, "_"), i)] <- 1
      }
      misdata <- misdata[-1]
      data <- cbind(data, misdata)
      data <- data[, c(names(data)[1:3], sort(names(data)[-(1:3)]))]
    }
    if(dname == "weights" & !is.null(weights)){
      data <- cbind(data, weights)
    }
  }
  return(data)
}


#' \code{get_start_cont} get starting values for continuous equations.
#' @import systemfit
#' @import DEoptim
#' @import maxLik
#' @import AER
#' @param of Function to me optimized.
#' @param hf Hessian function.
#' @param datcont Data-set used in the optimization.
#' @param eq_c Equations that will be estimated with \code{nlsystemfit}.
#' @param par_c Names of parameters.
#' @param startvals Supplied starting values.
#' @param DEoptim_run if \code{TRUE}, \code{DEoptim} is run for optimization.
#' @param deconst Absolute value of lower and upper bounds for \code{DEoptim}.
#' @return Starting values for continuous equations.
#' @keywords internal
get_start_cont <- function(of, hf, datcont, eq_c, par_c, startvals = NULL, DEoptim_run = TRUE,
                           deconst = 2)
{
  . <- NULL
  #check if all arguments supplied
  mf <- match.call()
  varl <- paste0(match.call()[[1]],'%>%args%>%as.list')%>%parse(text=.)%>%eval
  no_def <- varl[varl==""]%>%names#no default values
  m <- match(no_def, mf%>%names, 0L)
  if(any(m==0)){
    paste0(no_def[m==0], collapse = ", ")%>%paste0("No value for: ", .)%>%stop
  }

  npar <- length(par_c)
  #function for DEoptim
  of_DEoptim <- function(par){
    hfr <- try(diag(solve(-hf(par))), silent=TRUE)
    if((any(class(hfr)%in%"try-error")|any(hfr%in%c(Inf, -Inf, NA, NaN)))){
      res <- Inf
    }else{
      if(any(hfr<0)){
        res <- Inf
      }else{
        if(is.null(hf)){
          hfr1 <- try(maxLik(of, start = par, method="BFGS", iterlim=10000), silent=TRUE)
        }else{
          hfr1 <- try(maxLik(of, hess = hf, start = par, method="BFGS", iterlim=10000), silent=TRUE)
        }

        if(any(class(hfr1)%in%"try-error") )browser()
        res <- -of(par)
      }
    }
    res
  }
  startv <-  rep(0.45, npar)
  names(startv) <- par_c

  if(DEoptim_run==TRUE){
    a <- DEoptim(of_DEoptim, lower=rep(-deconst, npar),upper=rep(deconst, npar),
                 control=list(trace=FALSE, VTR=10000000))
  }
  startv <- if(!is.null(startvals))  startvals else startv

  start1 <- startv
  a1 <- try(nlsystemfit(method="SUR", eqns=lapply(eq_c, as.formula), startvals = start1,
                        data=datcont), silent=TRUE)
  if(any(class(a1)%in%"try-error")){
    a <- DEoptim(of_DEoptim, lower=rep(-deconst, npar),upper=rep(deconst, npar),
                 control=list(trace=FALSE, VTR=10000000))
    startv <-  a$optim$bestmem
    names(startv ) <- par_c
    
    start1 <- startv
    a1 <- try(nlsystemfit(method="SUR", eqns=lapply(eq_c, as.formula), startvals = start1,
                      data=datcont), silent=TRUE)
    if(any(class(a1)%in%"try-error")){
      a1 <- c()
      a1$b <- startv
    }
  }else{
    
  }
  aL1 <- suppressWarnings(of_DEoptim(a1$b))
  if(DEoptim_run==TRUE){
    aL <- of_DEoptim(a$optim$bestmem)
    aa <- if(aL < aL1) a$optim$bestmem else a1$b
  }else{
    aa <- a1$b
  }
  aa
}

#' \code{prepare_data} prepare data for the estimation.
#' @import data.table
#' @importFrom abind abind
#' @importFrom tidyr spread gather unite separate
#' @param data \code{data.frame} 
#' @param choice Name of variable with modes.
#' @param dummy Name of variable indicating, if the mode was chosen.
#' @param PeID Name of variable with individual identification numbers.
#' @param WeID Name of variable with trip identification.
#' @param type Type of data. If "long", then modifications are done.
#' @param mode_spec_var Used if format "long", mode specific variables.
#' @param avl if \code{TRUE}, includes dummies for mode availability.
#' @param chc if \code{TRUE}, includes dummies for choice of mode.
#' @param wc if \code{TRUE}, creates weights 1 for continuous equations.
#' @param wd if \code{TRUE}, creates weights 1 for discrete equations.
#' @param nc Integer, number of continuous equations.
#' @param weights Data matrix with weights, 
#' column names have to be $wc_i$(continuous), $wd_i$(discrete).
#' @param weight_paths if \code{TRUE}, weight according to number of trips per person, discrete part.
#' @param weight_paths_cont if \code{TRUE}, weight continuous part.
#' @param mode_factors if choice is not factor or numeric, this is important to supply.
#' @return data.frame used for modeling.
#' @examples
#' data("TravelMode", package = "AER")
#' mode_spec_var <- c("wait", "vcost", "travel", "gcost")
#' res <- prepare_data(TravelMode, choice="mode", dummy="choice", PeID="individual", WeID="",
#' type="long", mode_spec_var =mode_spec_var, nc=3)
#' @export
prepare_data <- function(data, choice = "", dummy = "", PeID = "", WeID = "",  type = "",
                         mode_spec_var = "",  avl = TRUE, chc = TRUE, wc = TRUE, wd = TRUE, nc = 0,
                         weights = NULL, weight_paths = FALSE, weight_paths_cont = FALSE,
                         mode_factors = NULL)
{
  . <- NULL
  choicen <- NULL
  npaths <- NULL
  npaths_cont <- NULL
  data_weight <- NULL
  mode_code <- c()
  if(choice!=""){
    if(!is.null(mode_factors)){
      data$choicen <- data[, choice]%>%factor(., levels=mode_factors)%>%as.numeric
    }else{
      data$choicen <- data[, choice]%>%as.factor%>%as.numeric
    }
    mode_code <- data[, c(choice, "choicen")] %>% data.frame%>%unique
    if(!is.null(mode_factors)){
      cat("Used mode coding:", "\n")
      print(mode_code[order(mode_code$choicen),])
    }
    data[, choice] <- NULL 
  }
  if(dummy!=""){
    data$dummy <- if(any(class(data[, dummy])%in%
                         c("character", "factor", "logical"))) data[, dummy]%>%
      as.factor %>% as.numeric - 1 else data[, dummy]
    if(length(data$dummy %>% unique)==1){
      data$dummy <- 1
    }
    data[, dummy] <- NULL
    
  }
  data$PeID <- if(!any(names(data)%in%"PeID")&PeID==""){
    1:nrow(data)
  }  else{
    xx <- data[, c(PeID, "PeID")[c(PeID, "PeID")%in%names(data)]]
    if(is.null(ncol(xx))){
      xx
    }else{
      if(ncol(xx)==2){
        xx[, 1]
      }else{
        browser()
      }
    }

  }
  if(PeID!=""&PeID!="PeID")data[, PeID] <- NULL
  

  data$WeID <- if(!any(names(data)%in%"WeID")&WeID=="") {1
  } else{
    xx <- data[, c(WeID, "WeID")[c(WeID, "WeID")%in%names(data)]]
    if(is.null(ncol(xx))){
      xx
    }else{
      if(ncol(xx)==2){
        xx[, 1]
      } else{
        browser()
      }
    }
  }
  if(WeID!=""&WeID!="WeID")data[, WeID] <- NULL
  
  if(type=="long"){
    const_var <- names(data)[!names(data)%in%c("PeID", "WeID", mode_spec_var, "choicen", "dummy")]
    ctext <- const_var%>%paste0(., collapse=' + ')
    ctext <- if(ctext!="") paste0(" + ", ctext) else ctext
    vtext <- mode_spec_var%>%paste0('"',., '"')%>%paste0(., collapse=', ')
    setDT(data)
    choiced <- data[, .(PeID, WeID, dummy, choicen)]
    choiced <- choiced[dummy==1,]
    choiced[, dummy:=NULL]
  
    data <- paste0('dcast(data, PeID + WeID ', ctext,' ~ choicen, value.var=c(',vtext,'))')%>%
      parse(text=.)%>%eval
    data <- merge(data, choiced, all=TRUE, by=c("PeID","WeID"))
    #setnames(data, "choicen", "choice")
    data <- data.frame(data, stringsAsFactors = FALSE)
  }
  data <- data.table(data, stringsAsFactors = FALSE)

  if(weight_paths_cont==FALSE|weight_paths==TRUE){
    if(!any(names(data)%in%"npaths_cont")) data[, "npaths_cont"] <- 1
  }
  if(weight_paths==FALSE){
    if(!any(names(data)%in%"npaths")) data[, "npaths"] <- 1
  }
  if(!any(names(data)%in%"npaths")) data[, npaths:=length(WeID), by="PeID"]
  if(!any(names(data)%in%"npaths_cont")) data[, npaths_cont:=length(WeID), by="PeID"]

  data <- data.frame(data, stringsAsFactors = FALSE)
  if(chc==TRUE&!any(grepl("^chc_", names(data)))){
    xx <- try(setnames(data, "choicen", "choice"), silent = TRUE)
    data <- add_variable(data, "chc")
    if(!any(class(xx)%in%"try-error"))setnames(data, "choice", "choicen")
  }
  if(is.null(mode_code)&any(names(data)%in%"choice")){
    mode_code <- data.frame(choicen=1:length(unique(data[, "choice"])))
  }
  if(avl==TRUE&!any(grepl("^avl_", names(data)))){
    avld <-  matrix(1, nrow=nrow(data), ncol=nrow(mode_code))
    colnames(avld) <- paste0("avl_", mode_code$choicen)
    data <- cbind(data, avld)
  }


  if(!is.null(weights)){
    data <- add_variable(data = data, dname = "weights", weights = data_weight)
  }
  if(wc==TRUE&!any(grepl("^wc_", names(data)))&nc!=0){
    avld <-  matrix(1, nrow=nrow(data), ncol=nc)
    colnames(avld) <- paste0("wc_", 1:nc)
    data <- add_variable(data, dname = "weights",weights = avld)
  }
  if(wd==TRUE&!any(grepl("^wd_", names(data)))){
    avld <-  matrix(1, nrow=nrow(data), ncol=nrow(mode_code))
    colnames(avld) <- paste0("wd_",mode_code$choicen)
    data <- add_variable(data, dname = "weights",weights = avld)
  }
  data
}


#' Function for finding starting values for nls equations with errors
#' @keywords internal
generate_disc_coefs <- function(coefs, eqest, i = 1, type_select = "smallest") 
{
  . <- NULL
  para <- NULL
  parai <- NULL
  value <- NULL
  inder <- which(lapply(coefs, class)%in%"try-error")
  if(length(inder)==0) {
    inder <- 1
    mci <- coefs[[1]]
  }else{
    mci <- coefs[-inder][[1]]
  }
  
  for(ij in 1:length(inder)){
    #print(paste0("nls for equation: ", inder[ij], " produced an error."))
    get_par_sep <- eqest[[inder[ij]]] %>%deparse() %>% paste0(., collapse="") %>%
      gsub(".*~", "", .) %>%trimws %>%  gsub("I\\(", "", .) %>%
      gsub("\\)$", "", .) %>%  strsplit(  . ,"\\+|\\*|\\^|\\*\\*|\\-") %>%unlist %>% 
      gsub("\\)|\\(", "", .) %>% 
      trimws() %>% grep("^par\\d{1,10}", ., value=TRUE)
    if(length(mci)!=length(get_par_sep)){
      all_coef <- unlist(coefs[-inder])
      acd <- data.table(para = all_coef %>% names %>% gsub(".*\\.", "", .),
                        value = all_coef)
      acd <- acd[order(para)]
      acd[, parai := gsub("par", "", para) %>% as.numeric()]
      mci <- lapply(get_par_sep, function(x){
        xx <-  acd[para==x, value][i]
        indxx <- gsub("par", "", x) %>% as.numeric()
        if(type_select == "smallest"){
          if(is.na(xx)){
            xx <- acd[parai<indxx, value][i]
          }
          if(is.na(xx)){
            xx <- acd[parai>indxx, value][i]
          }
        }else if(type_select == "biggest"){
          if(is.na(xx)){
            xx <- acd[parai>indxx, value][i]
          }
          if(is.na(xx)){
            xx <- acd[parai<indxx, value][i]
          }
        }else{
          #set.seed(1234)
          xx <- sample(acd[, value], 1)
        }
        xx
      })
      mci <- unlist(mci)
    }
    names(mci) <- get_par_sep
    coefs[[inder[ij]]] <- mci
  }
  coefs
}

#' \code{get_start_disc} get starting values for discrete choice model.
#' @import AER
#' @import mlogit
#' @param MNf1 Function to me optimized
#' @param hff Hessian function
#' @param eqest Equations that will be estimated with nls
#' @param sdat Data-set used for in the optimization
#' @param ppardogit "dogit" parameters, NULL for "logit"
#' @param type_select method for filling missing starting values. Default is "smallest"(closest by the index to the left value).
#' "biggest" takes closest by index to the right value and "random" - samples randomly values.
#' @return Starting values for discrete equations.
#' @keywords internal
get_start_disc <- function(MNf1, hff, eqest, sdat, ppardogit = NULL, type_select = "smallest")
{
  . <- NULL
  #set.seed(123124)
  #check if all arguments supplied
  mf <- match.call()
  varl <- paste0(match.call()[[1]],'%>%args%>%as.list')%>%parse(text=.)%>%eval
  no_def <- varl[varl==""]%>%names#no default values
  m <- match(no_def, mf%>%names, 0L)
  if(any(m==0)){
    paste0(no_def[m==0], collapse = ", ")%>%paste0("No value for: ", .)%>%stop
  }

  coefs <- suppressWarnings(lapply(eqest, function(m) try(summary(nls(m, sdat))$coef[,1], silent = TRUE)))
  # if error for come equations write down the values of other equations
  
  coefs_o <- coefs
  if(any(lapply(coefs, class)%in%"try-error")){
    coefs <- generate_disc_coefs(coefs, eqest, i=1, type_select=type_select)
    # only works if the same number of parameters
  }


  help_start <- help_get_disc_start(coefs_o, coefs, hff, MNf1, eqest, ppardogit, i=1)
  help_start_o <- help_start
  
  gind <- help_start[[3]]
  
  if(sum(gind)==0){
    coefs <- generate_disc_coefs(coefs_o, eqest, i=1, type_select="biggest")
    help_start <- help_get_disc_start(coefs_o, coefs, hff, MNf1, eqest, ppardogit,  
                                      type_select = "biggest",  i=1)
    gind <- help_start[[3]]
    if(sum(gind)==0){
      coefs <- generate_disc_coefs(coefs_o, eqest, i=1, type_select="random")
      help_start <- help_get_disc_start(coefs_o, coefs, hff, MNf1, eqest, ppardogit,  
                                        type_select = "random",  i=1)
      gind <- help_start[[3]]
    }
    
    if(sum(gind)==0){
      help_start <- help_start_o 
    }
  }
  
  goodv <- help_start[[1]]
  valf <- help_start[[2]]
  gind <- help_start[[3]]
  startv <- help_start[[4]]
  
  if(length(goodv)>1&sum(gind)>0){
    goodv <- goodv[which.min(unlist(valf[gind]))]
  }else if(length(goodv)==0){
    print("No good starting values automatically found: can not invert the Hessian.")
    print("Optimization will be longer. You can stop the program and input starting values by yourself or wait.")
    goodv <- startv
  }
  goodv <- goodv[[1]][,2]
  names(goodv) <-  startv[[1]][,1]
  goodv
}

#' @keywords internal
help_get_disc_start <- function(coefs_o, coefs, hff, MNf1, eqest, ppardogit, 
                                type_select="smallest", i=1)
{
  . <- NULL
  coefs <- data.frame(parlp=gsub(".*[.]","",names(unlist(coefs))), stv=unlist(coefs),
                      stringsAsFactors = FALSE,row.names = NULL)
  coefs[, "parlp"] <- coefs[, "parlp"] %>% gsub(".*~", "", .) %>%trimws %>%  gsub("I\\(", "", .) %>%
    gsub("\\)$", "", .) %>% gsub("\\s.*", "", .) %>% gsub("\\*.*", "", .)
  coefs <- coefs[order(nchar(coefs$parlp), coefs$parlp),]
  duplc <- coefs[duplicated(coefs$parlp),"parlp"]%>%unique
  
  # if no same coefficients
  
  if(length(duplc)==0){
    startv <- list(coefs)
  }else{
    duf <- table(coefs$parlp)[duplc]
    unc <- coefs[!coefs$parlp%in%duplc,]
    sac <- coefs[coefs$parlp%in%duplc,]
    sac$ind <- paste0(sac$parlp,"_", unlist(lapply(duf, function(i)1:i)))
    
    
    inder <- TRUE
    ninder <- length(duf)
    nmodes <- length(eqest)

    # might be to many combinations, try removing the possibilities
    while(inder&any(duf>=2)){
      poscomb <- choose(length(sac$ind), ninder)
      if(poscomb<=10000){
        #possible combs
        inder <- FALSE
        combn <- try(combn(sac$ind, ninder), silent = TRUE)

      }else{
        # remove 2 par from each block
        for(ij in sac$parlp %>% unique){
          xx <- sac[sac$parlp==ij,]
          if(nrow(xx)>=2){
            xx <- xx[sample(1:nrow(xx), 2), ]
            sac <- rbind(sac[sac$parlp!=ij,], xx)
          }
        }
      }
      duf <- sac$parlp %>% table
    }
    
    combnn <- combn
    
    indc <- apply(combnn, 2, function(x){
      xx <- unlist(strsplit(x,split = "_"))%>%table
      if(any(xx[grep("par", names(xx))]>1)) FALSE else TRUE
    })
    
    combnn <- combnn[,indc]%>%matrix(., ncol=sum(indc))
    if(length(combnn)==0){
      combnn <- matrix(sac[,"ind"], nrow=1)
    }
    
    startv <- apply(combnn, 2, function(x){
      xx <- rbind(unc, sac[sac$ind%in%x, 1:2])
      xx[order(nchar(xx$parlp), xx$parlp),]
    })
    
    if(length(startv)>20){
      startv <- startv[sample(1:length(startv), 20)]
    }
  }
  # add for each alternative if dogit 1/n, ((1 + par[8] + par[9] + par[10] + par[11]) * (dalyba)
  if(!is.null(ppardogit)){
    startv <- lapply(startv, function(x){
      xd <- data.frame(parlp=ppardogit %>% gsub("\\[|\\]", "", .), stv=1/(length(ppardogit)))
      rbind(x, xd)
    })
  }
  
  valf <- lapply(startv,function(par) sum(MNf1(par[-1,2])))
  
  valh <- lapply(startv,function(par) try( tryCatch({
    qr.solve(-hff(par[-1,2]))},
    error = function(err){
      r <- solve(-hff(par[-1,2]))
      if(any(is.na(diag(r)))){
        diag(r) <- -1
      }
      return(r)}),
    silent = TRUE))
  
  gind <- !sapply(valh, function(hes){
    if(any(class(hes)%in%"try-error")) hes <- matrix(-1)
    any(diag(hes)<0)
  } )
  goodv <- startv[gind]
  
  if(length(goodv)>1){
    goodv <- goodv[which.min(unlist(valf[gind]))]
    res_goodv <- list(goodv, valf, gind, startv)
  }else if(length(goodv)==0&i<=length(coefs_o)){
    #print("No good starting values automatically found: can not invert the Hessian.")
    #print("Starting another combination")
    i <- i+1
    coefs <- generate_disc_coefs(coefs_o,eqest, i, type_select = type_select)
    res_goodv <- help_get_disc_start(coefs_o, coefs, hff, MNf1, eqest, ppardogit, type_select, i)
  }else{
    goodv <- startv
    res_goodv <- list(goodv, valf, gind, startv)
  }
  return(res_goodv)
}


#' @keywords internal
parf <- function(fname)
{
  . <- NULL
  paste0("args(",fname,")%>%as.list%>%names%>%paste0('# @param ', ., '\n')%>%cat")%>%
    parse(text=.)%>%eval
}

#' \code{get_start} get starting values for discrete or continuous choice model.
#' @param eq_c Continuous equations errors.
#' @param eq_d Discrete equations.
#' @param data \code{data.frame} is used in the optimization.
#' @param part Type of estimation: "joint", "cont", "disc".
#' @param datan Name of \code{data.frame} used in the optimization.
#' @param fixed_term if \code{TRUE}, includes fixed term in log-likelihood.
#' @param weight_paths if \code{TRUE}, weights paths of the whole system.
#' @param weight_paths_cont if \code{TRUE}, weight paths only in continuous part.
#' @param data_weight  \code{data.frame} with weights for continuous and discrete equations,
#'  same dim as \code{data}.
#' @param par_c Names of parameters in continuous equations.
#' @param par_d Names of parameters in discrete equations.
#' @param startvals Starting values, can be also \code{NULL}.
#' @param DEoptim_run if \code{TRUE}, runs DEoptim for the optimization.
#' @param hessian Name of hessian function.
#' @param best_method if \code{TRUE}, try all possible optimization methods and 
#' choose the one with the smallest likelihood.
#' @param transform if \code{TRUE}, quantile transformation is applied.
#' @param MNtypef "dogit" or "logit".
#' @param pardogit "dogit" parameters.
#' @param opt_method optimization method to use.
#' @param numerical_deriv if \code{TRUE}, numerical derivatives are calculated in nmm function
#' @return Starting values for discrete or continuous blocks.
#' @examples
#' # Example of discrete choice model
#' data("TravelMode", package = "AER")
#' eq_d <- c("ASC1 * 1 + B2_t * travel_1 + B3_v * vcost_1" ,
#'           "ASC2 * 1  + B2_t * travel_2 + B3_v * vcost_2",
#'           "ASC3 * 1  + B2_t * travel_3 + B3_v * vcost_3",
#'           "ASC4 * 1  + B2_t * travel_4 + B3_v * vcost_4")
#' parl <- c(paste0("ASC", 1:4), "B2_t", "B3_v")
#' obj <- get_par(parl, eq_d)
#'
#' mode_spec_var <- c("wait", "vcost", "travel", "gcost")
#' data <- TravelMode
#' data$wait <- as.numeric(data$wait)
#' data[data$wait==0,"wait"] <- 0.000001 # add a small number to 0
#' data$travel <- as.numeric(data$travel)
#' data[data$travel==0,"travel"] <- 0.000001
#' data$vcost <- as.numeric(data$vcost)
#' data[data$vcost==0,"vcost"] <- 0.000001
#' data <- prepare_data(data, choice="mode", dummy="choice", PeID="individual", WeID="",
#' type="long", mode_spec_var =mode_spec_var, wc=FALSE)
#' stv <- get_start(eq_d=eq_d, data=data, datan="data", part="disc", par_d = parl, 
#' transform = FALSE)
#' 
#' #example, system of equations
#' data("CreditCard", package="AER")
#' cdat <- CreditCard
#' cdat$income2 <- cdat$income^2
#' cdat$d_selfemp <- as.numeric(cdat$selfemp)
#' eq_c <- c("expenditure ~  b1*age + b2*income + b3*income2",
#' "income ~ a1*age + a2*d_selfemp + a3*dependents + a4*majorcards")
#' parl <- c(paste0("b", 1:3), paste0("a", 1:4))
#' para_cont <- get_par(parl, eq_c)
#' cheqs0 <- para_cont$cheqs0
#' \donttest{
#' stv <- get_start(eq_c = eq_c, data=cdat, datan="cdat", part="cont", par_c=parl)
#' }
#' @export
get_start <- function(eq_c = NULL, eq_d = NULL, data = NULL,  part = 'joint', datan = 'data',
                      fixed_term = FALSE, weight_paths = TRUE, weight_paths_cont = FALSE,
                      data_weight = NULL, par_c = NULL, par_d = NULL, best_method = FALSE,
                      startvals = NULL, DEoptim_run = FALSE, hessian = NULL, transform = TRUE,
                      MNtypef="logit", pardogit = NULL, opt_method = "BFGS", 
                      numerical_deriv = FALSE)
{
  . <- NULL
  #set.seed(1134252)
  if(!is.null(eq_c)){
    para_cont <- get_par(par_c, eq_c)
    cheqs0 <- para_cont$cheqs0
  }else{
    cheqs0 <- NULL
  }
  if(!is.null(eq_d)){
    para_disc <- get_par(par_d, eq_d)
    ffor <- para_disc$cheqs0
  }else{
    ffor <- NULL
  }
  indp <- !part%in%c("cont")
  data <- prepare_data(data, avl= indp, chc=indp, wd=indp, nc=length(cheqs0), wc=!indp)
  methodopt <- "BFGS"
  if(!is.null(cheqs0)&part!='disc'){
    npar <- length(par_c)
    datac <- 'datcont'
    
    # version 1
    mn <- maxle_p(cheqs0 = cheqs0, fixed_term = fixed_term)
    varcont <- strsplit(mn$formulap%>%gsub('_', 'XXX', .),'[[:punct:]]')%>%unlist%>%unique%>%
      trimws%>%gsub('XXX', '_', .)
    ncont <- c('PeID',paste0('wc_', 1:length(cheqs0)),names(data)[names(data)%in%varcont])%>%unique
    datcont <- data[, ncont]
    datcont <- unique(datcont)
    of <- f_create(mn$formula,data = datac,  cheqs0=cheqs0, transform = FALSE)
    hf <- grad_hess_eval(mn = mn, parnl = paste0('par',1:npar), hessian = TRUE,  data=datac,
                         cheqs0=cheqs0)
    environment(of)  <- environment(hf) <- environment()
    s_cont <- get_start_cont(of, hf,  datcont, eq_c, par_c, startvals=startvals,
                             DEoptim_run=DEoptim_run)
    # version 2
    mn <- maxle_p(cheqs0 = cheqs0, fixed_term = fixed_term, version2 =  TRUE)
    of <- f_create(mn$formula,data = datac,  cheqs0=cheqs0, transform = FALSE)
    hf <- grad_hess_eval(mn = mn, parnl = paste0('par',1:npar), hessian = TRUE,  data=datac,
                         cheqs0=cheqs0)
    environment(of)  <- environment(hf) <- environment()
    s_cont2 <- get_start_cont(of, hf,  datcont, eq_c, par_c, startvals=startvals,
                             DEoptim_run=DEoptim_run)
    
    v1 <- of(s_cont)
    v2 <- of(s_cont2)

    res_cont <- nmm(data = datcont, start_v=s_cont, eq_c=eq_c, par_c=par_c,
                         eq_type = "cont", opt_method = opt_method,
                         best_method =  best_method, hessian=hessian, 
                    numerical_deriv = numerical_deriv)
    if(v1!=v2){
      res_cont_v2 <- nmm(data = datcont, start_v=s_cont2, eq_c=eq_c, par_c=par_c,
                         eq_type = "cont", opt_method = opt_method,
                         best_method =  best_method, hessian=hessian, 
                         numerical_deriv = numerical_deriv)
      if(res_cont_v2$maximum>res_cont$maximum){
        res_cont <- res_cont_v2
      }
    }

    s_cont <- res_cont$estimate
  }else{
    s_cont <- NULL
  }
  if(!is.null(ffor)&part!='cont'){
    methodopt <- 'NA'
    #transform <- TRUE
    fforp <- gsub("\\[|\\]","", ffor)
    eqest <- sapply(paste0(paste0("chc_", 1:length(fforp))," ~ ", "I(",fforp,")" ), as.formula)
    if(MNtypef=="logit"){
      mnf1 <- MNlogitf(ffor, transform = transform, separatenmm = TRUE, weight_disc = TRUE)
      ppardogit <- NULL
    }else{
      ppardogit <- c(paste0('par[',(length(para_disc$parld)+1):(length(para_disc$parld)+length(pardogit)),']'))
      mnf1 <- MNdogitf(ffor, transform = transform, separatenmm = TRUE,
                       weight_disc = TRUE, ppardogit)
    }
    
    datad <- 'datdisc'

    vardisc <- strsplit(mnf1$formula%>%gsub('_', 'XXX', .),'[[:punct:]]')%>%unlist%>%unique%>%
      trimws%>%gsub('XXX', '_', .)
    ndisc <- c('PeID','WeID','choice',names(data)[names(data)%in%vardisc])
    datdisc<- data[, names(data)[names(data)%in%ndisc]]
    datdisc <- unique(datdisc)
    MNf1 <- f_create(mn=mnf1$formula, data=datad, cheqs0=cheqs0, fixed=1, probt=mnf1$probt,
                     tformula=mnf1$formulat, separatenmm = TRUE, sume=mnf1$sume)
    dnpar <- ffor%>%get_npar
    if(!is.null(pardogit)){
      dnpar <- dnpar+length(pardogit)
    }

    hf <- grad_hess_eval(mn=mnf1, parnl = paste0('par', 1:dnpar),
                         fixed=1, hessian = TRUE, data=datad, cheqs0=cheqs0)
    
    environment(MNf1) <- environment(hf) <- environment()
    separatenmm <- FALSE
    s_disc <- get_start_disc(MNf1, hf, eqest, datdisc, ppardogit=ppardogit)#calc hessian gradient
    
    if(!is.null(pardogit)){
      names(s_disc) <- c(para_disc$parld, pardogit)
    }else{
      names(s_disc) <- para_disc$parld
    }
    s_disc <- s_disc[-1]
    res_disc <- nmm(data=datdisc, start_v = s_disc, eq_d=eq_d, eq_type = "disc",
                         weight_paths = FALSE, par_d = par_d, opt_method = opt_method,
                            best_method =  FALSE, hessian=hessian, MNtypef = MNtypef, 
                    numerical_deriv = numerical_deriv)
    s_disc <- res_disc$estimate
  }else{
    s_disc <- NULL
  }
  s_joint <- c(s_cont, s_disc)
  s_joint
}