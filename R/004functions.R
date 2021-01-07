#' Maximum likelihood estimation of nonlinear multivariate models (NMM).
#' 
#' \code{nmm}, \code{nmm_sigma} and \code{summary} are the main functions used for the estimation of NMM.
#' \itemize{
#' \item{\code{nmm}}{ - Maximum likelihood estimation of nonlinear multivariate models (NMM)}
#' \item{\code{nmm_sigma}}{ - Optimizes the covariance matrix}
#' \item{\code{summary}}{ - returns summary of \code{nmm} object with "normal", "robust" or "clustered" standard errors. With option \code{new_coef} one can supply new coefficients and test their significance.}
#' }
#' 
#' @param data \code{data.frame} used in the optimization.
#' @param eq_type Possible options "joint", "cont", "disc".
#' @param eq_d Discrete equations.
#' @param eq_c Continuous equations.
#' @param par_c Parameters from continuous equations.
#' @param par_d Parameters from discrete equations.
#' @param start_v Starting values for optimization. 
#' If \code{NULL}, starting values are found by \code{get_start}  function.
#' @param check_hess If \code{TRUE}, check the Hessian.
#' @param corrl If \code{TRUE},  correlation between blocks (continuous and discrete).
#' @param weight_paths If \code{TRUE}, weight according to the number of choices made by individual i will be added to 
#' the whole system.
#' @param weight_paths_cont If \code{TRUE}, if only to continuous part should be weighted.
#' @param data_weight Data weight matrix.
#' @param estimate If \code{TRUE}, estimation is performed.
#' @param fixed_term If \code{TRUE}, includes fixed term to continuous equation block.
#' @param best_method If \code{TRUE}, all optimizers are checked. 
#' @param best_method4start If \code{TRUE}, all optimizers are checked for starting values.
#' @param DEoptim_run If \code{TRUE}, runs \code{DEoptim} in generation of starting values.
#' @param hessian String, name of the Hessian function.
#' @param print_out If \code{TRUE}, prints out log-likelihood for each equation.
#' @param diff_hessian If \code{TRUE}, for changing hessian and gradient with weights.
#' @param DEoptim_run_main If \code{TRUE}, run \code{DEoptim} in the main optimization.
#' @param deconst absolute value of lower and upper bound in \code{DEoptim} optimization.
#' @param bayesian_random If \code{TRUE}, than par[1] is changed to par[,1] to be used for 
#' optimization of random parameters in Bayesian estimation.
#' @param eqsys "sur" or "sem".
#' @param opt_method optimization method for \code{maxLik}.
#' @param try_last_DEoptim If \code{TRUE}, in case of error in \code{maxLik} should \code{DEoptim} be run.
#' @param numerical_deriv If \code{TRUE}, uses numerical derivative instead of the analytical.
#' @param miterlim Number many iterations passed to \code{maxLik} function
#' @param transform if \code{TRUE}, quantile transformation is applied to discrete equations.
#' @param MNtypef estimate "logit", or "dogit"
#' @param nmm_object \code{nmm} object created by \code{nmm} function.
#' @return 
#' \item{\code{nmm}}{returns \code{nmm} object with estimated parameters, functions, and data.}
#' \item{\code{nmm_sigma}}{returns estimated parameters, functions, data.}
#' \item{\code{summary}}{returns summary of \code{nmm} object.}
#' @examples
#' # estimation of System of Nonlinear Equations based on example from 'systemfit'
#' library(systemfit)
#' data( ppine , package="systemfit")
#' hg.formula <- hg ~ exp( h0 + h1*log(tht) + h2*tht^2 + h3*elev)
#' dg.formula <- dg ~ exp( d0 + d1*log(dbh) + d2*hg + d3*cr)
#' labels <- list( "height.growth", "diameter.growth" )
#' model <- list( hg.formula, dg.formula )
#' start.values <- c(h0=-0.5, h1=0.5, h2=-0.001, h3=0.0001,
#'                   d0=-0.5, d1=0.009, d2=0.25, d3=0.005)
#' model.sur <- nlsystemfit( "SUR", model, start.values, data=ppine, eqnlabels=labels )
#' eq_c <- as.character(c(hg.formula, dg.formula))
#' parl <- c(paste0("h", 0:3),paste0("d", 0:3))
#' res <- nmm(ppine, eq_c=eq_c, par_c=parl, start_v = start.values,
#' eq_type = "cont", best_method = FALSE, numerical_deriv=TRUE)
#' summary(res)
#' res_sigma_cont <- nmm_sigma(res,estimate=TRUE) # Estimation of the Variance-Covariance matrix
#' summary(res_sigma_cont)
#' 
#' #example discrete choice
#' library(mlogit)
#' data("Fishing", package = "mlogit")
#' Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")
#' ## a pure "conditional" model
#' mres <- summary(mlogit(mode ~ price + catch, data = Fish))
#' data <- prepare_data(Fish %>% data.frame %>% dplyr::select(-idx),
#' choice="alt", dummy="mode", PeID="chid", mode_spec_var = c("price", "catch"),
#'  type="long")
#' eq_d <- c("a1 + p1 * price_1 + p2 * catch_2", "a2 + p1 * price_2 + p2 * catch_2",
#'              "a3 + p1 * price_3 + p2 * catch_3", "a4 + p1 * price_4 + p2 * catch_4")
#' par_d <- c(paste0("a", 1:4), paste0("p", 1:2))
#' res <- nmm(data, eq_d=eq_d, par_d = par_d,  eq_type="disc", fixed_term=FALSE,
#' best_method=FALSE)
#' summary(res)
#' 
#' # joint estimation mockup example
#' data(dataM)
#' dataMp <- dataM %>% data.frame %>% prepare_data(. , choice="DR_Course", 
#' PeID = "Student")
#' eq_c <- c("PlcmtScore ~ exp(a0 + a1 * PSATM + a2 * Rank + a3 * Size)",
#' "ACTM ~ exp(c0 + c1 * GPAadj)")
#' par_c <- c(paste0("a", 0:3), paste0("c", 0:1))
#' eq_d <- c("ASC1" ,
#' "ASC2 + b1_2 * SATM + b2_2 * PlcmtScore",
#' "ASC3 + b1_3 * SATM + b2_3 * PlcmtScore")
#' par_d <- c(paste0("ASC", 1:3), paste0("b", rep(1:2, rep(2,2)), "_", 2:3))
#' \donttest{
#' nmm_joint_res <- nmm(dataMp, eq_type = "joint", eq_d = eq_d, 
#' par_d = par_d,  eq_c = eq_c, par_c = par_c, 
#' start_v = c(a0=3.394, a1=0.001, a2=-0.001, a3=0, c0=3.583, c1=-0.008, 
#' ASC2=-1.452, ASC3=3.047, b1_2=0.145, b1_3=0.102, b2_2=-0.133, b2_3=-0.168))
#' summary(nmm_joint_res)
#' }
#' @export
nmm <- function(data, eq_type = c("joint", "cont", "disc"), eq_d = NULL, eq_c = NULL, 
                     par_c = NULL, par_d = NULL, start_v = NULL,  
                     check_hess = TRUE, corrl = TRUE, weight_paths = TRUE,
                     weight_paths_cont = FALSE, data_weight = NULL, estimate = TRUE,
                     fixed_term = FALSE, best_method = FALSE, DEoptim_run = FALSE,
                     hessian = "joint_hess", print_out = FALSE, diff_hessian = FALSE,
                     bayesian_random = FALSE, DEoptim_run_main = FALSE, deconst = 2,
                     numerical_deriv = FALSE, best_method4start = FALSE, eqsys = "sur",
                     miterlim = 10000, opt_method = "BFGS", try_last_DEoptim = TRUE,
                     transform=NULL, MNtypef = "logit", nmm_object = NULL)
{
  . <- NULL
  PeID <- NULL
  mutate <- NULL
  WeID <- NULL
  objm <- NULL
  mnp <- NULL
  if(MNtypef=="dogit"&eq_type!="cont"){
    pardogit <- paste0("capt", 1:length(eq_d))
  }else{
    pardogit <- NULL
  }
  if(eq_type=="cont"){
    pardogit <- NULL
  }
  if(eq_type!="joint"){
    corrl <- FALSE
  }
  data %<>% data.frame(.)
  data %<>% replace(., is.na(.), 0)
  if(any(unlist(data)%in%c(NA, Inf, -Inf, NaN))){
    indx <- apply(data, 1, function(x) any(x%in%c(NA, Inf, -Inf, NaN)) ) %>% which
    nindx <- length(indx)
    indx <- sample(indx, min(6, nindx))
    print(paste0("Rows (sampled, in total ", nindx," rows): ", paste0(indx, collapse = ", "), " have non-normal values! Stopping."))
    stop()
  }
  
  
  if(!is.null(eq_c)){
    par_c %<>% sort 
    para_cont <- get_par(par_c, eq_c)
    cheqs0 <- para_cont$cheqs0
  }else{
    cheqs0 <- NULL
  }
  if(!is.null(eq_d)){
    #par_d_o <- par_d
    par_d %<>% sort
    ffor <-  get_par(par_d, eq_d)$cheqs0
  }else{
    ffor <- NULL
  }
  if(is.null(transform)){
    transform <- if(eq_type == "joint") TRUE else FALSE
  }
  
  indp <- eq_type!="cont"
  indc <- eq_type!="disc"
  data <- prepare_data(data, avl = indp, chc = indp, wd = indp, nc = length(cheqs0), wc = indc,
                       weights = data_weight, weight_paths = weight_paths,
                       weight_paths_cont = weight_paths_cont)
  
  if(!is.null(start_v)){
    start_v <- check_names_start(start_v)
    start_v <- start_v[c(par_c,par_d[-1])]
  }
  if(is.null(start_v)&estimate&MNtypef=="dogit"){
    mc <- match.call()
    mc[["MNtypef"]] <- "logit"
    #mc %>% deparse %>% gsub("dogit", "logit", .)
    start_v <- eval(mc)$estimate
    stp <- rep(1/length(pardogit), length(pardogit))
    names(stp) <- pardogit
    start_v <- c(start_v,stp)
  }else{
    if(is.null(start_v)&estimate){
      start_v <- get_start(eq_c = eq_c, par_c = par_c, par_d = par_d, eq_d = eq_d,
                           data = data, datan="data", part = eq_type, DEoptim_run = DEoptim_run,
                           hessian = NULL, best_method = best_method4start, transform=transform,
                           MNtypef=MNtypef, pardogit=pardogit,
                           opt_method = opt_method, numerical_deriv=numerical_deriv)
    }
  }
  
  

  methodopt <- opt_method
  data_name <- "data"
  separatenmm <- FALSE
  if(diff_hessian==TRUE){
    datan <- copy(data)
    wnames <- c(paste0('wd_', 1:length(ffor)), paste0('wc_', 1:length(cheqs0)))
    data_weightn <- data.frame(matrix(1, ncol=length(wnames), nrow=nrow(datan)),
                               stringsAsFactors = FALSE)
    names(data_weightn) <- wnames
    datan <- add_variable(data = datan, dname = "weights", weights = data_weightn)
    datan <- ddply(datan, .(PeID), mutate, npaths = length(WeID))
    datan$npaths_cont <- 1
    data_name <- "datan"
  }
  if(is.null(nmm_object)){
    if(is.null(start_v)){
      start_v <- rep(0, length(c(par_c, par_d)))
    }
    #correlation between cont and discrete
    cdmr <- get_nmm_functions(corrl, eq_type, cheqs0, ffor, objm, print_out,
                              bayesian_random, fixed_term, mnp, start_v, data_name, check_hess,
                              hessian, eqsys, diff_hessian, MNtypef=MNtypef)
  }else{
    cdmr <- attributes(nmm_object)$functions
    nn_cdmr <- c(paste0(rep(c("joint_", "cont_","text_"), rep(3, 3)), c("func", "grad", "hess")),paste0("text_", c("bd", "gd", "edd")))
    names(cdmr) <- nn_cdmr[1:length(names(cdmr))]
  }

  joint_func <- cdmr$joint_func
  joint_grad <- cdmr$joint_grad
  joint_hess <- cdmr$joint_hess
  separatenmm <- cdmr$separatenmm
  if(corrl==TRUE&any(eq_type=="joint")){
    cont_func <- cdmr$cont_func
    cont_grad <- cdmr$cont_grad
    cont_hess <- cdmr$cont_hess
    text_func <- cdmr$text_func
    text_grad <- cdmr$text_grad
    text_hess <- cdmr$text_hess
    text_bd <- cdmr$text_bd
    text_gd <- cdmr$text_gd
    text_edd <- cdmr$text_edd
    environment(cont_func) <- environment(cont_grad) <- environment(cont_hess) <- environment()
  }
  if(corrl==FALSE){
    separatenmm <- FALSE
  }
  environment(joint_func) <- environment(joint_grad) <- environment(joint_hess) <- environment()
  
  # list of functions, likelihood, gradient, hessian, cont. gradient, cont.hessian
  functionsl <- if(corrl == TRUE&any(eq_type == "joint")){
    list(jfunc = joint_func, jgrad = joint_grad, jhess = joint_hess,
         cfunc = cont_func, cgrad = cont_grad, chess = cont_hess,
         text_func = text_func, text_grad = text_grad, text_hess=text_hess,
         text_bd=text_bd, text_gd=text_gd, text_edd=text_edd)
  } else{
    list(jfunc = joint_func, jgrad = joint_grad, jhess = joint_hess)
  }

  # if estimation is needed, default FALSE
  if(estimate == TRUE){
    # numerical derivative instead of analytical, default FALSE
    if(numerical_deriv){
        joint_grad <- NULL
        joint_hess <- NULL
        check_hess <- FALSE # used in jfunc to check hessian for inversion
    }
    # possible maxLik optimization methods

    possible_m <- c("BFGS", "NR",  "CG", "NM", "SANN", "BHHH")
    # search for best results between all methods
    if(best_method == TRUE){

      res <- nmm_best(joint_func, joint_grad, joint_hess, start_v, "NA", miterlim, possible_m,
                      eq_type, corrl)
      resif <- try_if(res)
    }else{
      # if only run 1 optimization method, default BFGS
      if(DEoptim_run_main){
        # if DEoptim should be used for the main optimization of log-likelihood
        res <- nmm_DEoptim(joint_func, joint_grad, joint_hess, start_v, "NA", miterlim,
                                possible_m, deconst, eq_type, corrl)
        resif <- try_if(res)
        MLerror <- resif
      }else{
        # DEoptim not used for main optimization
        MLerror <- FALSE
      }
      # if DEoptim not used for main optimization (nothing yet optimized) or
      # ML produced error (DEoptimization done)
      if(DEoptim_run_main == FALSE|MLerror){
        res <- nmm_1good(joint_func, joint_grad, joint_hess, start_v, "NA", miterlim,
                              possible_m, methodopt, eq_type, corrl)
      }
      resif <- try_if(res)
    }
    # if still no solution last attempt with DEoptim
    resif <- try_if(res)
    if(resif&try_last_DEoptim){
      stde <- if(is.null(res$estimate)) start_v
      stde <- if (is.na(sum(joint_func(start_v)))) stde else start_v
      if(!is.na(sum(joint_func(start_v)))&!is.na(sum(joint_func(res$estimate)))){
        stde <- if(joint_func(start_v)>joint_func(res$estimate)) start_v else res$estimate
      }

      res <- nmm_lastDEoptim(joint_func, joint_grad, joint_hess, stde, "NA", miterlim,
                                  possible_m, deconst, eq_type, corrl)
      resif <- try_if(res)
    }

    # calc. log-likelihood
    estj <- joint_func(res$estimate)
    # separate log-likelihood
    estja <- attributes(estj)[[1]]
    if(length(estja)>0){
      names(estja) <- paste0('eq_', 1:length(estja))
    }else{
      estja <- estj
    }

    res$gradient <- gradient.nmm(res)
    res$hessian <- hessian.nmm(res)
  }else{
    # in case of no estimation
    res <- "object"
    estja <- NA
  }

  a <- c(attributes(res), list(functions = functionsl, data = data, LL = estja,
                               type = eq_type, cont_e = cheqs0, disc_e = ffor,
                               eq_c = eq_c, eq_d = eq_d, par_c = par_c,
                               par_d = par_d, corr_joint=corrl))
  
  if(any(class(res)%in%"maxLik")){
    res <- res
  }else{
    res <- list(res)
  }

  attributes(res) <- a
  class(res) <- c("nmm", class(res))
  res$call$data <- "data"
  if(data %>% dplyr::select(dplyr::starts_with("chc")) %>% ncol>0){
    res$freq <- try(data %>% dplyr::select(dplyr::starts_with("chc")) %>% colSums, silent = TRUE)
  }else{
    res$freq <- NULL
  }
  res
}

#' Checks if starting values have names
#' @keywords internal
check_names_start <- function(start_v){
  nstv <- names(start_v)
  if(is.null(nstv)){
    print("Provide names for `start_v`!")
    stop()
  }else{
    start_v[order(names(start_v))]
  }
}

#' searches for best optimization method between possible maxLik algorithms
#' @keywords internal
nmm_best <- function(joint_func, joint_grad, joint_hess, start_v, method = NULL, miterlim,
                          possible_m, eq_type, corrl)
{
  . <- NULL
  print('optimization searches for the best result.')
  # estimation with all methods

  if(is.null(joint_hess)){
    reslist <- lapply(possible_m, function(x) try(maxLik(joint_func, start = start_v,
                                                         method = x, iterlim = miterlim),
                                                  silent=TRUE))
  }else{
    reslist <- lapply(possible_m, function(x) try(maxLik(joint_func, grad = joint_grad,
                                                         hess = joint_hess, start = start_v,
                                                         method = x, iterlim = miterlim),
                                                  silent=TRUE))
  }

  res <- reslist[[1]]
  # check all methods for possible errors
  indr <- sapply(reslist, try_if)
  #remove the ones with errors
  reslist <- reslist[!indr]
  if(length(reslist)>0){
    indr <- sapply(reslist, function(x)x$maximum)
    res <- reslist[[which.max(indr)[1]]]
  }
  return(res)
}

#' searches for 1st good optimization method between possible maxLik algorithms
#' @keywords internal
nmm_1good <- function(joint_func, joint_grad, joint_hess, start_v, method = NULL, miterlim,
                           possible_m, methodopt = "BFGS", eq_type, corrl)
{
  . <- NULL
  res <- try(try_alg(joint_func, joint_grad, joint_hess, methodim = methodopt, miterlim, start_v,
                     ki = 1, k = 2, eq_type=eq_type, corrl=corrl), silent = TRUE)
  im <- 1
  resif <- try_if(res)
  # remove the methods that was already used and produced errors
  possible_m <- possible_m[!possible_m%in%methodopt]
  # if error try another methods
  while(resif&im <= length(possible_m)){
    print(paste0("Algorithm ", ifelse(im == 1, methodopt, possible_m[im-1]),
                 " produced an error. Trying out next algorithm."))
    # repeat algorithm 2 times
    res <- try(try_alg(joint_func, joint_grad, joint_hess, methodim=possible_m[im], miterlim,
                       start_val = start_v, ki = 1, k = 2, eq_type=eq_type, corrl=corrl),#start_val=res$estimate
               silent=TRUE)
    resif <- try_if(res)
    if(resif){
      # model with start values from initial generation
      res <- try(try_alg(joint_func, joint_grad, joint_hess, methodim=possible_m[im], miterlim,
                         start_v, ki = 1, k = 1, eq_type=eq_type, corrl=corrl), silent=TRUE)
      resif <- try_if(res)
    }
    im <- im+1
  }
  return(res)
}

#' searches for 1st good optimization method between possible maxLik algorithms and then applies DEoptim
#' @keywords internal
nmm_DEoptim <- function(joint_func, joint_grad, joint_hess, start_v = NULL, method = NULL,
                             miterlim, possible_m, deconst = 2, eq_type, corrl)
{
  . <- NULL
  print('DEoptim is used for optimization. This may take a while.')
  # generate initialpop with maxLik -> Deoptim with const +-deconst ->
  # maxLik with values from DEoptim -> Deoptim with maxLik values, used if maxLik produced errors
  # starting population


  res <- nmm_1good(joint_func = joint_func, joint_grad=joint_grad, joint_hess=joint_hess, start_v=start_v, method="NA", miterlim=miterlim, possible_m=possible_m,
                     eq_type=eq_type, corrl=corrl, methodopt=method)
  res1 <- nmm_1good(joint_func = joint_func, joint_grad=joint_grad, joint_hess=joint_hess, start_v=res$estimate, method="NA", miterlim=1000000, possible_m=possible_m,
                   eq_type=eq_type, corrl=corrl, methodopt="NM")
  
  resif <- try_if(res)

  # add random error to population
  initialpop <- matrix(c(res1$estimate ) %>% rep(.,length(start_v)), ncol = length(start_v), byrow = TRUE)
  initialpop <- rbind(initialpop - rep(-0.0001, length(start_v)) %>% diag, initialpop + rep(-0.0001, length(start_v)) %>% diag)
  
  initialpop <- rbind(res$estimate,res1$estimate,initialpop)

  
  
  #function used in case of DEoptim optimization in nmm without Hessian check
  jj_f <- function(par){
    res <- joint_func(par)
    -res
  }  
  environment(jj_f) <- environment()

  deconst <- max(ceiling(max(abs(res$estimate))), deconst)
  reso <- DEoptim(jj_f, lower = rep(-deconst, length(start_v)),
                  upper = rep(deconst, length(start_v)),
                  control = list(trace = TRUE, initialpop = initialpop, itermax=1000))
  # rerun maxLik
  st_new <- reso$optim$bestmem
  names(st_new) <- names(start_v)
  if(is.null(joint_hess)){
    res <- try(maxLik(joint_func, start = st_new, iterlim = miterlim, method = "BFGS"), silent = TRUE)
    resnm <- try(maxLik(joint_func, start = res$estimate, iterlim = 1000000, method = "NM"), silent = TRUE)
  }else{
    res <- try(maxLik(joint_func, grad = joint_grad, hess = joint_hess,
                      start = st_new, iterlim = miterlim, method = "BFGS"), silent = TRUE)
    resnm <- try(maxLik(joint_func, grad = joint_grad, hess = joint_hess, 
                        start = res$estimate, iterlim = 1000000, method = "NM"), silent = TRUE)
  }
  
  if(res$maximum<resnm$maximum){
    res <- resnm
  }

  resif <- try_if(res)
  if(resif){
    # rerun DEoptim, values used in case res is an error
    reso <- DEoptim(jj_f, lower = res$estimate-0.5, upper = res$estimate + 0.5,
                    control = list(trace = FALSE))
    st_new <- reso$optim$bestmem
    names(st_new) <- names(start_v)
    # modify res object with reso values
    res$estimate <- st_new
    res$gradient <- joint_grad(st_new)
    res$hessian <- joint_hess(st_new)
    res$type <- "DEoptim"
    res$maximum <- joint_func(st_new)
  }
  return(res)
}

#' last time try out DEoptim if no functions before produced "good" results
#' @keywords internal
nmm_lastDEoptim <- function(joint_func, joint_grad, joint_hess, start_v = NULL, method = NULL,
                                 miterlim, possible_m, deconst = 2, eq_type, corrl)
{
  . <- NULL
  paste0("All maxLik algorithms produced errors. Trying out DEoptim. This may take a while.")
  jj_f <- function(par){
    res <- joint_func(par)
    -res
  }  
  environment(jj_f) <- environment()
  if(!is.null(start_v)){
    if(any(is.na(start_v))){
      start_v[is.na(start_v)] <- 0
    }
    # set.seed(1987)
    # add random error to population
    initialpop <- matrix(c(start_v + rnorm(length(start_v) * 1000, sd = 0.01)),
                         ncol = length(start_v), byrow = TRUE)
    deconst <- start_v %>% abs %>% max + 1
    reso <- try(DEoptim(jj_f, lower = rep(-deconst, length(start_v)),
                    upper = rep(deconst, length(start_v)),
                    control = list(trace = FALSE, itermax = 1000, VTR = 10000000,
                                   initialpop = initialpop)), silent = TRUE)
    if(any(class(reso)%in%"try-error")){
      vals_pop <- apply(initialpop, 1, jj_f) 
      initialpop <- initialpop[!is.nan(vals_pop),]
      if(nrow(initialpop)==0){
        stop("No initial population values for DEoptim found!")
      }else{
        reso <- DEoptim(jj_f, lower = rep(-deconst, length(start_v)),
                        upper = rep(deconst, length(start_v)),
                        control = list(trace = FALSE, itermax = 1000, VTR = 10000000,
                                       initialpop = initialpop))  
      }
    }
    st_new <- reso$optim$bestmem
    names(st_new) <- names(start_v)
    if(is.null(joint_hess)){
      res <- try(maxLik(joint_func, start = st_new, iterlim = miterlim, method = "BFGS"), silent = TRUE)
      resnm <- try(maxLik(joint_func, start = res$estimate, iterlim = 1000000, method = "NM"), silent = TRUE)
    }else{
      res <- try(maxLik(joint_func, grad = joint_grad, hess = joint_hess,
                        start = st_new, iterlim = miterlim, method = "BFGS"), silent = TRUE)
      resnm <- try(maxLik(joint_func, grad = joint_grad, hess = joint_hess, 
                          start = res$estimate, iterlim = 1000000, method = "NM"), silent = TRUE)
    }
    
    if(res$maximum<resnm$maximum){
      res <- resnm
    }

    resif <- try_if(res)
    if(resif){
      reso <- try(DEoptim(jj_f, lower = rep(-deconst, length(start_v)),
                      upper = rep(deconst, length(start_v)),
                      control = list(trace = FALSE, VTR = 10000000)), silent = TRUE)
    }
  }else{
    reso <- try(DEoptim(jj_f, lower = rep(-deconst, length(start_v)),
                    upper = rep(deconst, length(start_v)),
                    control = list(trace = FALSE, VTR = 10000000)), silent = TRUE)
  }

  if(!any(class(reso)%in%"try-error")){
    st_new <- reso$optim$bestmem
    names(st_new) <- names(start_v)
  }
  

  if(is.null(joint_hess)){
    res <- try(maxLik(joint_func, start = st_new, iterlim = miterlim, method = "BFGS"), silent = TRUE)
    resnm <- try(maxLik(joint_func, start = res$estimate, iterlim = 1000000, method = "NM"), silent = TRUE)
  }else{
    res <- try(maxLik(joint_func, grad = joint_grad, hess = joint_hess,
                      start = st_new, iterlim = miterlim, method = "BFGS"), silent = TRUE)
    resnm <- try(maxLik(joint_func, grad = joint_grad, hess = joint_hess, 
                        start = res$estimate, iterlim = 1000000, method = "NM"), silent = TRUE)
  }
  
  if(res$maximum<resnm$maximum){
    res <- resnm
  }

  resif <- try_if(res)
  if(resif){
    res$estimate <- st_new
    if(!is.null(joint_grad)) res$gradient <- joint_grad(st_new)
    if(!is.null(joint_hess)) res$hessian <- joint_hess(st_new)
    res$type <- "DEoptim"
    res$maximum <- joint_func(st_new)
  }
  return(res)
}


#' Checks if nmm object has an error
#' @keywords internal
#' @param res maxLik object for which summary and other statistics have to be checked
try_if <- function(res)
{
  . <- NULL
  sumstat <- try(res%>%summary_check, silent = TRUE)
  tester <- ifelse(any(class(res)%in%"try-error"), TRUE, ifelse(res$code == 100, TRUE, FALSE))
  ststats <- if(any(sumstat%>%class%in%"try-error")|tester) TRUE else
    all(summary_check(res)$estimate[, 4]%>%na.omit %>% round(., 4)>0.999)|any(is.nan(summary_check(res)$estimate[, 4]))
  return(tester|any(sumstat%>%class%in%"try-error")|ststats)
}

#' Test if maxLik produces an error
#' @param methodim - method to try out
#' @param start_v - starting values
#' @param ki - iteration number
#' @param k - how many times to try
#' @keywords internal
try_alg <- function(joint_func, joint_grad, joint_hess, methodim, miterlim, start_val,
                    ki = 1, k = 1, eq_type, corrl = FALSE)
{
  . <- NULL
  env_org <- environment(joint_func)
  method_org <- env_org$methodopt
  env_org$methodopt <- methodim
  
  
  if(is.null(joint_hess)){
    res <- try(maxLik(joint_func, start = start_val,
                      method = methodim, iterlim = miterlim), silent = TRUE)
  }else{
    res <- try(maxLik(joint_func, grad = joint_grad, hess = joint_hess, start = start_val,
                      method = methodim, iterlim = miterlim), silent = TRUE)
  }
  env_org$methodopt <- method_org
  
  if(eq_type=="joint"&corrl==FALSE&!is.null(joint_hess)){
    h <- joint_hess(res$estimate) 
    if(methodim=="BHHH"){
      h %<>% apply(., 3, function(x) colSums(x, na.rm = FALSE))
    }
    res$hessian <- h 
  }
  resif <- try_if(res)
  if(resif&ki<k){
    ki <- ki+1
    res <- try_alg(joint_func, joint_grad, joint_hess, methodim, miterlim,
                   start_val = res$estimate, ki, k, eq_type=eq_type, corrl=corrl)
    return(res)
  }
  return(res)
}

#' @keywords internal
get_nmm_functions <- function(corrl, eq_type, cheqs0, ffor, objm, print_out,
                                   bayesian_random, fixed_term, mnp, start_v, data_name, check_hess,
                                   hessian, eqsys, diff_hessian, MNtypef = NULL)
{
  . <- NULL
  separatenmm <- FALSE
  if(corrl == TRUE&any(eq_type == "joint")){
   # conditional mean and variance expressions

   objm <- cond_mean_cov_expr(neqt = length(cheqs0), kk = length(ffor))
   joint_obj <- LL_joint(ffor, cheqs0, objm,  check_hess = TRUE , cfunc='cont_func',
                            chess = "cont_hess", cgrad = "cont_grad",
                            gfunc = 'joint_grad', hfunc = 'joint_hess', print_out =  print_out,
                            bayesian_random = bayesian_random, MNtypef=MNtypef)
    
    text_func <- joint_obj$text_func
    text_grad <- joint_obj$text_grad
    text_hess <- joint_obj$text_hess
    text_bd <- joint_obj$text_bd
    text_gd <- joint_obj$text_gd
    text_edd <- joint_obj$text_edd
    
    joint_func <- joint_obj$func
    joint_grad <- joint_obj$grad #gradient
    joint_hess <- joint_obj$hessian #hessian

    npar <- get_npar(cheqs0)+1 # number of prameters in cont. equations
    add_pars <- paste0("par[", npar:(length( start_v)+1), "]")
    
    
    
    # get continuous functions
    cont_obj <- LL_joint_no_corr(cheqs0=cheqs0, datan="data", eq_type = "cont",
                                  fixed_term = fixed_term,
                                  bayesian_random = bayesian_random,
                                  eqsys = eqsys, add_pars = add_pars)
    mod_act <- .%>% deparse %>% gsub("npaths_cont", "1", ., fixed = TRUE) %>% 
      gsub("npaths", "1", ., fixed = TRUE) 
    cf <- cont_obj$func %>% mod_act
    intr <- cf %>% grep("separatenmm", .)-1
    if(length(intr)>0){
      cf <- cf[-c(intr:(length(cf)-2))]
    }
    
    
    cont_func <- function(par){}
    body(cont_func)<- parse(text=cf[-c(1)] )

    cf <- cont_obj$grad %>% mod_act
    cont_grad <- function(par){}
    body(cont_grad)<- parse(text=cf[-c(1)] )
    
    cf <- cont_obj$hess %>% mod_act
    cont_hess<- function(par){}
    body(cont_hess)<- parse(text=cf[-c(1)] )

    if(bayesian_random == TRUE){
      actbr <- .%>%gsub("par\\[\\[(\\d{1,10})\\]\\]", "par[,\\1]", .)%>%gsub("par\\[(\\d{1,10})", "par[,\\1",.)%>%
        gsub("par\\[\\((\\d{1,10})", "par[,(\\1",.)%>%
        gsub("par <- c(", "par <- cbind(", ., fixed = TRUE)%>%
        gsub("length(par)", "ncol(par)", ., fixed = TRUE)

      dpactbr <- . %>% deparse %>% actbr
      ffj <- cont_func %>% dpactbr
      ffj <- ffj[-c(1)] %>% parse(text = .)
      cont_func <- function(par) {}
      body(cont_func) <- ffj

      gf1 <- cont_grad %>% dpactbr
      gf1 <- gf1[-c(1)] %>% parse(text = .)
      cont_grad <- function(par){}
      body(cont_grad) <- gf1

      hf1 <- cont_hess %>% dpactbr
      hf1 <- hf1[-c(1)] %>% parse(text = .)
      cont_hess <- function(par){}
      body(cont_hess) <- hf1
    }

    if(diff_hessian == TRUE){
      joint_obj <- LL_joint(ffor, cheqs0, objm,  check_hess = TRUE, cfunc='cont_func',
                            chess = "cont_hess", cgrad = "cont_grad",
                            gfunc = 'joint_grad', hfunc = 'joint_hess', datan = data_name,
                            print_out =  print_out,
                            bayesian_random = bayesian_random,
                            MNtypef=MNtypef)
      joint_grad <- joint_obj$grad #gradient
      joint_hess <- joint_obj$hessian #hessian
    }
  }else{
    #no correlation
    joint_obj <- LL_joint_no_corr(ffor, cheqs0, datan="data", eq_type = eq_type,
                                  fixed_term = fixed_term, hessian = hessian,
                                  bayesian_random = bayesian_random,
                                  eqsys = eqsys, MNtypef=MNtypef)
    joint_func <- joint_obj$func
    joint_grad <- joint_obj$grad #gradient
    joint_hess <- joint_obj$hessian #hessian
  }
  
  if(bayesian_random == TRUE){
    tjf <- joint_func %>% deparse
    istart <- grep("methodopt", tjf)[1]
    if(corrl == TRUE){
      iend <- grep("joint_hess", tjf)[1]+1
      tjf <- tjf[-c(istart:iend)]
      istart <- grep(" res <- if ", tjf)
      iend <- grep("attributes\\(res\\)", tjf)
      tjf <- tjf[-c(istart:iend)]
      tjf <- gsub("contpart <- cbind\\(contpart", "res <- cbind(contpart", tjf)
      body(joint_func) <-  parse(text = tjf[-1])
      
      # modify cont_func
      cjf <- cont_func %>% deparse
      istart <- grep("methodopt", cjf)[1]
      iend <- grep("res <- cbind", cjf)-1
      cjf <- cjf[-c(istart:iend)]
      body(cont_func) <-  parse(text = cjf[-1])
    }else{
      iend <- grep("attributes\\(res\\)", tjf)
      ires <- grep("res <- cbind", tjf, value = TRUE)
      tjf <- tjf[-c(istart:iend)]
      body(joint_func) <-  parse(text = c(tjf[2:(length(tjf)-2)], ires, tjf[(length(tjf)-1):length(tjf)]))
    }


  }
  
  
  if(corrl==TRUE&any(eq_type=="joint")){
    res <-
      list(
        joint_func = joint_func,
        joint_grad = joint_grad,
        joint_hess = joint_hess,
        separatenmm = separatenmm,
        cont_func = cont_func,
        cont_grad = cont_grad,
        cont_hess = cont_hess,
        text_func = text_func,
        text_grad = text_grad,
        text_hess = text_hess,
        text_bd = text_bd,
        text_gd = text_gd,
        text_edd = text_edd
      )
  }else{
    res <-
      list(
        joint_func = joint_func,
        joint_grad = joint_grad,
        joint_hess = joint_hess,
        separatenmm = separatenmm
      )
  }
  res
}

#' @rdname nmm 
#' @param object \code{nmm} object, for \code{summary} and \code{nmm_sigma}
#' @param type Type of standard errors c("robust", "clustered", "normal"), for \code{summary}
#' @param new_coef New coefficients that will be tested, for \code{summary}
#' @param ... additional arguments affecting the summary produced, for \code{summary}
#' @export
summary.nmm <- function(object, type = "normal", new_coef = NULL, ...)
{
  . <- NULL
  new_coef <- new_coef[names(new_coef) %>% sort]
  estimate <- if(!is.null(new_coef)) new_coef else object$estimate

  breadd <- bread.nmm(object, type = type, new_coef = new_coef)
  sde <-  breadd %>% diag %>% sqrt
  t <- estimate / sde
  p <- 2 * pnorm(-abs(t))
  results <- cbind(Estimate = estimate, `Std. error` = sde, `t value` = t, `Pr(> t)` = p)
  results <- apply(results, 2, function(x)round(x, 7))
  ll <- logLik.nmm(object, new_coef = new_coef)

  summary <- list(maximType = object$type, iterations = object$iterations,
                  returnCode = object$code, returnMessage = object$message,
                  loglik = ll, estimate = results,
                  constraints = object$constraints, vcov = breadd)
  if(!is.null(new_coef)){
    summary[["maximType"]] <- "New"
    summary[["iterations"]] <- "NA"
    summary[["returnCode"]] <- "NA"
    summary[["returnMessage"]] <- "New coefficients"
  }
  class(summary) <- "summary.maxLik"
  return(summary)
}

#' Helper functions
#' 
#' Produce function that calculate estimates of endogenous variables from
#' the continuous block and probabilities from discrete part.
#' 
#' @import gsubfn
#' @param eq_c continuous equations
#' @param eq_d discrete equations
#' @param par_c parameters from cont. eq
#' @param par_d parameters from disc. eq
#' @param fixed index of fixed parameter
#' @return Returns functions.
#' @examples
#' eq_d <- c("ASC1 * 1 + B11_dur * dur_1" , "ASC2 * 1 + B12_dur * dur_2",
#' "ASC3 * 1 + B13_dur * dur_3")
#' eq_c <- c("Tw ~ tw*w + ph1*Tc", "Tf1 ~ (1+w)^tw + ph1^3*Tc")
#' parl <- c("tw", "ph1")
#' par_d <- c(paste0("ASC", 1:3), paste0("B1", 1:3, "_dur"))
#' stfunc <- stats_function(eq_c, eq_d, parl,par_d, fixed=3)
#' data <- matrix(runif(1000, min=0.001, max=50), ncol=8)
#' data <- data.frame(data)
#' names(data) <- c("dur_1", "dur_2", "dur_3", "w", "Tc", "avl_1", "avl_2", "avl_3")
#' parv <- c(0.5, 1, 1.5, 2, 1, -0.3, 0.2, -0.8)
#' methodopt <- "BHHH"
#' separatenmm <- TRUE
#' env <- environment()
#' fnames <- c("prob_func", "cont_func")
#' sapply(1:(length(stfunc)-1), function(x)assign(fnames[x], stfunc[[x]], envir=env))
#' eval(parse(text=paste0("environment(", fnames[1:(length(stfunc)-1)], ") <- env")))
#' probs <- prob_func(parv)
#' apply(probs, 2, mean)
#' cont <- cont_func(parv)
#' apply(cont, 2, mean)
#' @export
stats_function <- function(eq_c = NULL, eq_d = NULL, par_c = NULL, par_d = NULL, fixed = 0)
{
  . <- NULL
  if(!is.null(eq_c)){
    para_cont <- get_par(par_c, eq_c)
    npar <- length(par_c)
    cheqs0 <- para_cont$cheqs0
    cheqs0n <- cheqs0%>%sub('^[^-]*-', '', .)
    endog <- cheqs0 %>% sub('-.*', '', .) %>% trimws
    fixed <- if(!is.null(eq_d)) npar + 1 else 0
    cont_func <- f_create(cheqs0n, data = 'data', fixed = fixed, transform = FALSE,
                          cheqs0 = cheqs0n)
  }else{
    cont_func <- NULL
    endog <- NULL
    npar <- 0
    fixed <- if(!is.null(eq_d)) npar+1 else 1
  }
  if(!is.null(eq_d)){
    disc_par <- get_par(par_d, eq_d)
    ffor <- disc_par$cheqs0
    fforn <- gsubfn('(par\\[)(\\d{1,10})(\\])', function(x, y, z) paste0(x, as.numeric(y) + npar, z),
                    ffor)
    mnf1 <- MNlogitf(fforn, transform = TRUE, separatenmm = TRUE, weight_disc = TRUE)
    probexpr <- mnf1$probte %>% gsub('qnorm\\(|\\)$', '', .)
    prob_func <- f_create(probexpr, data = 'data', fixed = fixed, transform = FALSE)
  }else{
    prob_func <- NULL
  }
  list(prob = prob_func, cont = cont_func, endog = endog)
}
