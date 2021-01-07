#' @rdname nmm
#' @param methodopt optimizer from \code{maxLik} package, for \code{nmm_sigma}
#' @param try_1good If \code{TRUE}, stops then first good values are found, for \code{nmm_sigma}
#' @param try_DEoptim If \code{TRUE}, uses \code{DEoptim} for optimization, 
#' then \code{maxLik} optimizers produce errors, for \code{nmm_sigma}
#' @param try_diff_method If \code{TRUE}, stops then first good values with Hessian check are found, for \code{nmm_sigma}
#' @param trace If \code{TRUE}, trace of \code{DEoptim} is printed, for \code{nmm_sigma}
#' @export
nmm_sigma <- function(object, methodopt="BFGS", try_1good = TRUE, try_DEoptim = FALSE,
                      try_diff_method = FALSE, trace = FALSE, estimate = FALSE){
  . <- NULL
  check_hess <- FALSE
  separatenmm <- FALSE
  
  par <- object$estimate
  names(attributes(object))
  eq_type <- attributes(object)$type
  par_c <-  attributes(object)$par_c
  par_d <- attributes(object)$par_d
  data <- attributes(object)$data
  corr_joint <- attributes(object)$corr_joint
  funcl <- attributes(object)$functions
  joint_func <- funcl$jfunc
  cont_func <- funcl$cfunc
  text_func <- funcl$text_func
  eq_c <-  attributes(object)$eq_c
  eq_d <- attributes(object)$eq_d
  joint_hess <- function(par){}
  
  neq <- length(eq_c)
  if(neq==0){
    print("Variance-covariance matrix is not available for logit!")
    stop()
  }
  ndisc <- length(eq_d)
  sigma_func <- calc_sigma_rho(object)
  
  mobj <-  mod_function(object)
  mod_func <- mobj[[1]]
  cont_func1 <- mobj[[2]]
  text_func1 <- mobj[[3]]
  environment(joint_func) <- environment()
  #joint_func(par)
  environment(sigma_func) <- environment()
  environment(mod_func) <- environment()
  try(environment(cont_func1) <- environment(), silent=TRUE)
  
  rr <- sigma_func(par) # in corr_joint==TRUE  sigma is variance not standard deviation!!!
  start_v <- sigma2par(rr, corr_joint)
  if(corr_joint){
    start_v[1:neq] <- sqrt(start_v[1:neq])
  }
  
  if(estimate) {
    # cont_func(par)-cont_func1(start_v[1:(neq+1)])
    # mod_func(start_v) - joint_func(par)
    # mod_func(start_v) == joint_func(par) # should be the same
    res <-  maxLik(mod_func, start = start_v, method = methodopt)
    poss_method <- c("BFGS", "NR", "NM", "CG", "SANN")
    elevels <- c(1e-10, 1e-8, 1e-6, 1e-5, 1e-4)
    
    # function with Hessian checking used in DEoptim
    mod_func_hess <- function(start_v){
      LL <- mod_func(start_v)
      # numerical hessian 
      nmh <- numericNHessian(mod_func, t0 = start_v, eps = 1e-5)
      hh <- try(qr.solve(-nmh), silent = TRUE)
      if(any(class(hh)%in%"try-error")) {
        -Inf
      }else{
        if(any((hh %>% diag)%in%c(NA, NaN))){
          -Inf
        }else{
          if(any((hh %>% diag)<0)){
            -Inf
          }else{
            LL
          }
        }
      }
    }
    
    if(try_if(res)&try_1good){
      print("Trying different optimization methods and tolerance. ")
      res1 <- nmm_sigma_1good(res, mod_func = mod_func, poss_method, elevels, start_v)
      if(!try_if(res1)){
        print("Good values found with nmm_sigma_1good.") 
      }else{
        print("No values found with nmm_sigma_1good.")
      }
      res <- res1
    }
    
    
    
    ############ if still error try DEoptim
    if(try_if(res)&try_DEoptim){
      print("Trying DEoptim.")
      resDE <- nmm_sigma_DEoptim(res, mod_func, mod_func_hess, start_v = res$estimate, 
                                 itermax = 100, trace=trace)
      if(resDE$maxim!=-Inf){
        start_v_DE <- resDE$estimate
        nmh <- numericNHessian(mod_func, t0 = resDE$estimate, eps = 1e-5)
        nmg <- numericGradient(mod_func, t0 = resDE$estimate, eps = 1e-5)
        res1 <- res
        res1$estimate <- resDE$estimate
        res1$hessian <- nmh
        res1$gradient <- nmg
        summary(res1)
        if(!try_if(res1)){
          res <- res1
          print("Good values found with nmm_sigma_DEoptim.") 
        }else{
          print("No values found with DEoptim.") 
        }
      }
    }
    #################################
    # try find values by checking the Hessian in optimization
    
    if(try_if(res)&try_diff_method){
      print("Trying with hessian check.")
      res1 <-  try(maxLik(mod_func_hess, start = start_v, method = "BFGS"), silent=TRUE)
      k <- 1
      while(any(class(res)%in%"try-error")&k<=length(poss_method)){
        k <- k+1
        res1 <- try(maxLik(mod_func_hess, start = start_v, method = poss_method[k]), silent=TRUE)
      }
      if(!try_if(res1)){
        print("Good values found hessian check.") 
        res <- res1
      }else{
        print("No values found hessian check.") 
      }
    }
  }else{
    res <- list()
    res$estimate <- start_v
  }
  res <- get_est_sigma_matrix(res, neq, corr_joint, ndisc, estimate = estimate)
  res
}

#' @keywords internal
nmm_sigma_DEoptim <- function(res, mod_func, mod_func_hess, start_v,
                              itermax=100, trace = FALSE){
  jj_h <- function(start_v){
    rr <- -mod_func(start_v)
    if(any(rr)%in%c(NA, NaN)){
      Inf
    }else{
      rr
    }
  }
  ll <- rep(-1, length(start_v))
  uu <- rep(1, length(start_v))
  ll[names(start_v)!=""] <- 0.0001
  uu[names(start_v)!=""] <- max(abs(start_v))+10
  reso <- DEoptim(jj_h, lower =ll,
                  upper = uu,
                  control = list(trace = trace, itermax=itermax))
  LL <- mod_func_hess(reso$optim$bestmem)
  estimate <- reso$optim$bestmem
  names(estimate) <- names(start_v)
  list(maximum=LL, estimate = estimate)
}

#' @keywords internal
nmm_sigma_1good <- function(res, mod_func, poss_method, elevels, start_v){
  nn <- length(poss_method)*length(elevels)
  k <- 1
  j <- 1
  while(k<=nn&try_if(res)){
    nmh <- numericNHessian(mod_func, t0 = res$estimate, eps = elevels[j])
    nmg <- numericGradient(mod_func, t0 = res$estimate, eps = elevels[j])
    res$hessian <- nmh
    res$gradient <- nmg
    if(!try_if(res)) k <- nn + 1
    if(k%%length(elevels)==0&k<nn){
      j <- 1
      res <-  try(maxLik(mod_func, start = start_v, 
                         method = poss_method[k/length(elevels)+1]), silent = TRUE)
      if(any(class(res)%in%"try-error")) {
        k <- ((k)%/%length(elevels)+1)*(length(elevels))
        while(any(class(res)%in%"try-error")&k<=nn){
          res <-  try(maxLik(mod_func, start = start_v, 
                             method = poss_method[k/length(elevels)+1]), silent = TRUE)
          k <- ((k)%/%length(elevels)+1)*(length(elevels))

        }
      }else{
        warning("Check k!")
      }
    }else{
      j <- j+1
    }
    k <- k+1
  }
  res
}


#' Converts vector of variances and correlations values into a matrix
#' @param par Vector of variances and correlations
#' @param neq Number of continuous equations
#' @param corr_joint Correlation between blocks
#' @param ndisc Number of discrete choice alternatives
#' @return Vector 
#' @keywords internal
par2sigma <- function(par, neq,  corr_joint, ndisc){
  . <- NULL
  if(corr_joint==FALSE&corr_joint==FALSE){
    sigma <- par[1:neq]
    rho <- diag(1, nrow=neq, ncol=neq)
    rho[lower.tri(rho)] <- par[(neq+1):length(par)]
    rho[upper.tri(rho)] <- t(rho)[upper.tri(rho)]
  }
  if(corr_joint==TRUE){
    #check for lower matrix
    ncont <- neq
    sigma <- diag(c(par[1:ncont],rep(1, ndisc)),  nrow=ncont+ndisc, ncol=ncont+ndisc)
    rho <- diag(1, nrow=ncont+ndisc, ncol=ncont+ndisc)
    rrl <- lower.tri(rho)
    rrl[,(ncont+1):nrow(rrl)] <- FALSE
    rho[rrl] <-  par[(ncont+1):length(par)]
    rho[1:ncont, 1:ncont][upper.tri(rho[1:ncont, 1:ncont])] <- t(rho[1:ncont, 1:ncont])[upper.tri(rho[1:ncont, 1:ncont])]
    rho[1:ncont, (ncont+1):nrow(rrl)] <- t(rho[(ncont+1):nrow(rrl), 1:ncont] )
  }
  list(sigma, rho)
}


#' Converts sigma matrix into parameter vector
#' @param obj nmm object
#' @keywords internal
sigma2par <- function(x, corr_joint){
  . <- NULL
  if(corr_joint==FALSE){
    res <- c(x$sigma, x$rho[lower.tri(x$rho)])
  }
  if(corr_joint==TRUE){
    ncont <- x$sigma %>% colnames %>% grep("eta", .) %>% length
    ndisc <- x$sigma %>% colnames %>% grep("transf", .) %>% length
    ss <- x$sigma[1:ncont, 1:ncont]
    if(ncont==1){
      ss <- as.matrix(ss, ncol=1, nrow=1)
      colnames(ss) <- rownames(ss) <- "eta1"
    }
    rr <- x$rho
    rr[(ncont+1):ncol(rr), (ncont+1):ncol(rr)] <- diag(1, ncol = ndisc, nrow=ndisc) # no correlation between discrete equations!
    rrl <- lower.tri(rr)
    rrl[,(ncont+1):nrow(rrl)] <- FALSE
    res <- c(ss%>% diag, rr[rrl])
  }
  res
}



#' Modifies function to optimize sigma and rho
#' @param obj nmm object
#' @return Function 
#' @keywords internal
mod_function <- function(obj){
  . <- NULL
  eq_type <- attributes(obj)$type
  corr_joint <- attributes(obj)$corr_joint
  funcl <- attributes(obj)$functions
  joint_func <- funcl$jfunc
  cont_func <- funcl$cfunc
  text_func <- funcl$text_func
  
  
  mod_func <- function(ppp){}
  if((corr_joint==FALSE&corr_joint==FALSE)|eq_type!="joint"){
    dfc <- deparse(joint_func) 
    ind <- dfc %>% grep("sigma <- sqrt", ., fixed = TRUE)
    body(mod_func) <- c(dfc[2:ind],  
                        "rte <- par2sigma(ppp, neq, corr_joint, ndisc)",
                        "sigma <- rte[[1]]",
                        "rho <- rte[[2]]",
                        dfc[(ind+1):length(dfc)]) %>% parse(text=.)
    cont_func1 <- NULL
    text_func1 <- NULL
  }
  if(corr_joint==TRUE&eq_type=="joint"){
    # modify cont_func
    cont_func1 <- function(ppp){}
    dfc <- deparse(cont_func) 
    ind <- dfc %>% grep("sigma <- sqrt", ., fixed = TRUE)
    body(cont_func1) <- c(dfc[2:ind],  
                          "rte <- par2sigma(ppp, neq, corr_joint=FALSE, ndisc)",
                          "sigma <- rte[[1]]",
                          "rho <- rte[[2]]",
                          dfc[(ind+1):length(dfc)]) %>% parse(text=.)
    # modify text_func
    text_func1 <- text_func
    ind <- text_func1 %>% grep("rho <- cor", ., fixed = TRUE)
    text_func1 <- c(text_func1[1:ind],  
                    " ppp1 <- ppp",
                    " ppp1[1:neq] <- ppp1[1:neq]^2",
                    "rte <- par2sigma(ppp1, neq, corr_joint, ndisc)",
                    "sigma <- rte[[1]]",
                    "rho <- rte[[2]]",
                    text_func1[(ind+1):length(text_func1)])
    
    # modify joint_func
    dfc <- deparse(joint_func) 
    dfc %<>% gsub("text_func", "text_func1",.)
    dfc %<>% gsub("text_func11", "text_func1",.)
    ind <- dfc %>% grep("contpart <- cont_func", ., fixed = TRUE)
    dfc[ind] <-"        contpart <- cont_func(ppp1) * (1/npaths) * (1/npaths_cont)"
    dfc <- c(dfc[1:(ind-1)],
             "        ppp1 <- c(((sigma %>% diag %>% sqrt)[1:neq]), rho[1:neq, 1:neq][lower.tri( rho[1:neq, 1:neq])] )",
             dfc[(ind):length(dfc)])
    dfc %<>% gsub("cont_func", "cont_func1",., fixed = TRUE)
    dfc %<>% gsub("cont_func11", "cont_func1",., fixed = TRUE)
    body(mod_func) <- dfc[2:length(dfc)]%>% parse(text=.)
    
  }
  list(mod_func, cont_func1, text_func1)
}

#' Returns function for sigma and rho calculation
#' @param obj nmm object
#' @return Function that can be evaluated
#' @keywords internal
calc_sigma_rho <- function(obj){
  . <- NULL
  eq_type <- attributes(obj)$type
  corr_joint <- attributes(obj)$corr_joint
  funcl <- attributes(obj)$functions
  joint_func <- funcl$jfunc
  cont_func <- funcl$cfunc
  text_func <- funcl$text_func
  
  sigma_func <- function(par){}
  dfc <- deparse(joint_func) 
  if(corr_joint==FALSE|eq_type!="joint"){
    ind <- dfc %>% grep("sigma <- sqrt", ., fixed = TRUE)
    text_sigma_body <- c(dfc[2:ind],  
                          "colnames(rho) <- NULL",
                          "rownames(rho) <- NULL",
                          "names(sigma) <- NULL",
                          "list(rho=rho, sigma=sigma)","}") 
  }
  if(corr_joint==TRUE&eq_type=="joint"){
    # sigma in text_joint
    indf <- dfc %>% grep("text_func", ., fixed = TRUE)# where text_func starts
    ind <- text_func %>% grep("rho <- cor", ., fixed = TRUE)
    tend <- "})}"
    text_sigma_body <- c(dfc[2:(indf-1)], 
                          text_func[1:ind],
                          "colnames(rho) <- NULL",
                          "rownames(rho) <- NULL",
                          "names(sigma) <- NULL",
                          "list(rho=rho, sigma=sigma)",tend) 
    
  }
  if(attributes(obj)$cont_e %>% unlist %>% length == 1){
    #if only 1 continuous equation this is needed
    ind <- grep("rho", text_sigma_body)[1]
    text_sigma_body <- c(text_sigma_body[1:(ind-1)], "rho <- matrix(1)", text_sigma_body[ind:length(text_sigma_body)])
  }
  body(sigma_func) <- text_sigma_body%>% parse(text=.)
  sigma_func
}

#' @keywords internal
get_est_sigma_matrix <- function(sobj, neq, corr_joint, ndisc, estimate = FALSE){
  pl <- par2sigma(sobj$estimate, neq, corr_joint , ndisc)
  sigma <- pl[[1]]
  if(class(sigma)=="matrix"){
    sigma <- diag(sigma)[1:neq]
  }
  names(sigma) <- paste0("s_", 1:neq)
  rho <- pl[[2]]
  if(corr_joint==FALSE){
    nnmaes <- paste0("rho_", 1:neq)
  }
  if(corr_joint==TRUE){
    nnmaes <- c(paste0("rho_", 1:neq),paste0("d_", 1:ndisc))
  }
  colnames(rho) <- row.names(rho) <- nnmaes
  if(estimate){
    # sd functions
    pls <- par2sigma(summary(sobj)$estimate[, 2], neq, corr_joint , ndisc)
    sigmas <- pls[[1]]
    if(class(sigmas)=="matrix"){
      sigmas <- diag(sigmas)[1:neq]
    }
    names(sigmas) <- paste0("sd_", 1:neq)
    rhos <- pls[[2]]
    colnames(rhos) <- row.names(rhos) <- nnmaes
    #plls <- list(rho=rhos, sigma=sigmas)
    # t values
    plst <- par2sigma(summary(sobj)$estimate[,3], neq, corr_joint , ndisc)
    sigmast <- plst[[1]]
    if(class(sigmast)=="matrix"){
      sigmast <- diag(sigmast)[1:neq]
    }
    names(sigmast) <- paste0("sd_", 1:neq)
    rhost <- plst[[2]]
    colnames(rhost) <- row.names(rhost) <- nnmaes
    # p vals
    plp <- par2sigma(summary(sobj)$estimate[,4], neq, corr_joint , ndisc)
    sigmap <- plp[[1]]
    if(class(sigmap)=="matrix"){
      sigmap <- diag(sigmap)[1:neq]
    }
    names(sigmap) <- paste0("sd_", 1:neq)
    rhop <- plp[[2]]
    colnames(rhop) <- row.names(rhop) <- nnmaes
    
    #make par
    xx <- rbind(data.frame(par=paste0("sd_", 1:neq), est=sigma, sd=sigmas, tval=sigmast, pval=sigmap),
                make_par(rho, rhos, rhost, rhop))
    
    npar <- xx$par[match(sobj$estimate, xx$est)] %>% as.character()
    names(sobj$estimate) <- npar
  }else{
    sobj <-  list(sd=sigma, correlation=rho)
  }
  sobj

}

#' @keywords internal
make_par <- function(rho, rhos, rhost, rhop){
  . <- NULL
  xx <- c()
  for(i in 1:nrow(rhos)){
    for(j in 1:nrow(rhos)){
      if(i <= j) next
      ij <- colnames(rho)[i] %>% gsub("rho_", "c_", .)
      ji <- colnames(rho)[j]%>% gsub("rho_", "c_", .)
      xx <- rbind(xx, data.frame(par=paste0("corr(",ji, "&", ij,")"), est=rho[i, j], sd=rhos[i, j],
                                 tval=rhost[i, j], pval=rhop[i, j]))
    }
  }
  xx
}
