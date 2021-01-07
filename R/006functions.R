#' Add interactions
#'
#' \code{addInter} add interactions into continuous equations.
#' @param eqcont Vector of strings containing equations.  
#' @param par_c Names of coefficients.
#' @param intv Vector of integers corresponding to coefficients to which interactions
#'  should be added.
#' @param inter_parl Names of new coefficients (interactions).
#' @return list: 1 - expressions of errors, equations, parameters to estimate
#' @examples
#' eq_c <- c("Tw ~ tw*w + ph1*Tc", "Tf1 ~ (1+w)^tw + ph1^3*Tc")
#' parl <- c("tw", "ph1")
#' intv <- c(1,0)
#' inter_parl <- c('yytw','yyph1')
#' res <- addInter(eq_c, parl, intv, inter_parl)
#' @export
addInter <- function(eqcont, par_c, intv, inter_parl)
{
  . <- NULL
  para_cont <- get_par(par_c, eqcont)
  cheqs0 <- para_cont$cheqs0

  parnames <- paste0("par[",1:length(par_c), "]" )
  tt <- inter_parl%>%gsub("^y", "", .)
  newN <- paste0("Z", tt)#c("tw", "PH", "t1", "p1")
  newNN <- paste0("y", tt)
  names(parnames) <- par_c
  index <- which(intv==1)
  if(length(index)>0){
    for(ij in index){
      ii <- which(index==ij)
      pn <- parnames[index[ii]]
      pn1 <- names(parnames)[index[ii]]
      cheqs0 <- gsub(pn, paste0("", pn, " * ", newN[index[ii]], "**par[",length(par_c)+1,"]"),
                     cheqs0, fixed = TRUE)
      par_c <- c(par_c,  newNN[index[ii]])
      eqcont <- gsub(pn1, paste0("",pn1, " * ", newN[index[ii]], "**", newNN[index[ii]],""),
                     eqcont, fixed = TRUE)
    }
  }
  list(cheqs0, eqcont, par_c)
}



#' Log-likelihood(LL) with supplied coefficients.
#' @param object Object of class \code{nmm}.
#' @param new_coef "New" coefficients for which LL should be calculated.
#' @param separatenmm if \code{TRUE}, returns separate LL for each equation.
#' @param transform if \code{TRUE}, do quantile transformation (normal quantiles).
#' @param methodopt "NA" means that automatic algorithm was used in \code{maxLik}, if equal to
#'  "BHHH" will return LL for each individual.
#' @param ... some methods for this generic function require additional arguments.
#' @return Returns log-likelihood.
#' @examples
#' #example  continuous nonlinear
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
#' res <- nmm(ppine, eq_c=eq_c, start_v=start.values, par_c=parl, eq_type = "cont",
#'  best_method = FALSE)
#' logLik(res)
#' logLik(res, new_coef=res$estimate)
#' logLik(res, new_coef=model.sur$b)
#' #example discrete
#' library(mlogit)
#' data("Fishing", package = "mlogit")
#' Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")
#' ## a pure "conditional" model
#' mres <- summary(mlogit(mode ~ price + catch, data = Fish))
#' data <- prepare_data(Fish %>% data.frame %>% dplyr::select(-idx),
#' choice="alt", dummy="mode", PeID="chid", mode_spec_var = c("price", "catch"),
#' type="long")
#' eq_d <- c("a1 + p1 * price_1 + p2 * catch_2", "a2 + p1 * price_2 + p2 * catch_2",
#'              "a3 + p1 * price_3 + p2 * catch_3", "a4 + p1 * price_4 + p2 * catch_4")
#' par_d <- c(paste0("a", 1:4), paste0("p", 1:2))
#' res <- nmm(data, eq_d=eq_d,  eq_type="disc", fixed_term=FALSE, par_d=par_d,
#' best_method=FALSE)
#' logLik(res)
#' logLik(res, new_coef=res$estimate)
#' logLik(res, new_coef=mres$coefficients)
#' @export
logLik.nmm <- function(object, new_coef=NULL, separatenmm = FALSE, transform = FALSE,
                        methodopt = "NA", ...)
{
  . <- NULL
  joint_func <- NULL
  if(is.null(new_coef)){
    object$maximum
  }else{
    data <- attributes(object)$data
    ffunc <- attributes(object)$functions
    check_hess <- FALSE
    envr <- environment()
    namesf <- c(paste0("joint_",c("func","grad", "hess")), paste0("cont_",c("func","grad", "hess")),
                paste0("text_",c("func","grad", "hess", "bd", "gd", "edd")))
    aa <- sapply(1:length(ffunc), function(x)assign(namesf[x], ffunc[[x]], envir=envr))
    for(i in 1:length(ffunc)){
      try(eval(parse(text=paste0("environment(", namesf[i], ") <- envr"))), silent=TRUE)
    }
    joint_func(new_coef)
  }
}


#' Hessian with supplied coefficients
#' @keywords internal
#' @param res result from nmm
#' @param new_coef coefficients for which logLik should be calculated
#' @param separatenmm return separatenmm LL for each equation
#' @param transform to do quantile transformation
#' @param methodopt "NA" means that automatic algorithm was used in maxLik, if equal to "BHHH" will
#' return LL for each individual
hessian.nmm <- function(res, new_coef = NULL, separatenmm = FALSE, transform = FALSE,
                         methodopt = "NA")
{
  . <- NULL
  joint_hess <- NULL
  if(is.null(new_coef)){
    res$hessian
  }else{
    data <- attributes(res)$data
    ffunc <- attributes(res)$functions
    env <- environment()
    namesf <- c(paste0("joint_",c("func","grad", "hess")),
                paste0("cont_",c("func","grad", "hess")),
                paste0("text_",c("func","grad", "hess", "bd", "gd", "edd")))
    aa <- sapply(1:length(ffunc), function(x)assign(namesf[x], ffunc[[x]], envir=env))
    for(i in 1:length(ffunc)){
      try(eval(parse(text=paste0("environment(", namesf[i], ") <- envr"))), silent=TRUE)
    }
    joint_hess(new_coef)
  }
}

#' Gradient with supplied coefficients
#' @keywords internal
#' @param res result from nmm
#' @param new_coef coefficients for which logLik should be calculated
#' @param separatenmm return separatenmm LL for each equation
#' @param transform to do quantile transformation
#' @param methodopt "NA" means that automatic algorithm was used in maxLik, if equal to "BHHH" will
#' return LL for each individual
#' @param ind_grad get individual gradient values
gradient.nmm <- function(res, new_coef = NULL, separatenmm = FALSE, transform = FALSE,
                          methodopt = "NA", ind_grad = FALSE)
{
  . <- NULL
  joint_grad <- NULL
  if(is.null(new_coef)&ind_grad==FALSE){
    res$gradient
  }else{
    if(is.null(new_coef)){
      new_coef <- res$estimate
    }
    data <- attributes(res)$data
    ffunc <- attributes(res)$functions
    env <- environment()
    namesf <- c(paste0("joint_",c("func","grad", "hess")),
                paste0("cont_",c("func","grad", "hess")),
                paste0("text_",c("func","grad", "hess", "bd", "gd", "edd")))
    aa <- sapply(1:length(ffunc), function(x)assign(namesf[x], ffunc[[x]], envir=env))
    for(i in 1:length(ffunc)){
      try(eval(parse(text=paste0("environment(", namesf[i], ") <- envr"))), silent=TRUE)
    }
    joint_grad1 <- if(ind_grad == TRUE){
       function(par){
        methodopt <- "BHHH"
        environment(joint_grad) <- environment()
        try(environment(cont_grad) <- environment(), silent=TRUE)
        joint_grad(par)
      }
    }else{
      joint_grad
    }
    environment(joint_grad1) <- environment()
    joint_grad1(new_coef)
  }
}

#' \code{in2nmm} convert some estimation results into \code{nmm} object.
#' @param to \code{nmm} object.
#' @param new_coef New coefficients.
#' @return \code{nmm} object.
#' @examples
#' #example  continuous nonlinear
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
#' res <- nmm(ppine, eq_c=eq_c, start_v=start.values, par_c=parl, eq_type = "cont",
#' best_method = FALSE)
#' aa <- in2nmm(res, model.sur$b)
#' summary(res, new_coef=model.sur$b, type="robust")
#' summary(aa, type="robust")
#' summary(res, type="robust")
#' @export
in2nmm <- function(to, new_coef)
{
  . <- NULL
  new_coef <- new_coef[names(new_coef) %>% sort]
  object <- c()
  attrto <- attributes(to)

  eq_type <- attrto$type
  object$estimate <- new_coef
  object$fixed <- rep(FALSE, length(new_coef))

  object$activePar <- rep(TRUE, length(new_coef))

  estja <- attr(logLik(to, new_coef=new_coef), "LLv") 
  if(length(estja)>0){
    names(estja) <- paste0('eq_', 1:length(estja))
  }
  object$maximum <- estja %>% sum
  object$maximType <- "NA"
  object$iterations <- "NA"
  object$returnCode <- "NA"
  object$returnMessage  <- object$message <- "NA"
  object$code <- "NA"
  object$type <- "NA"
  object$gradient <- gradient.nmm(to, new_coef=new_coef)
  object$hessian <- hessian.nmm(to, new_coef=new_coef)
  object$constraints <- NULL
  class(object) <- c("nmm", "maxLik", "maxim","list")
  
  attributes(object) <- c(attributes(object),list(functions=attrto$functions, data=attrto$data,
                                                  LL=estja, type=attrto$type,
                                                  cont_e=attrto$cont_e,
                                                  disc_e=attrto$disc_e,
                                                  eq_c=attrto$eq_c, eq_d=attrto$eq_d,
                                                  par_c = attrto$par_c, par_d=attrto$par_d,
                                                  corr_joint= attrto$corr_joint))
  
  object$freq <- to$freq 
  object
}



#' Adjusted Akaike's Information Criterion.
#' 
#' Calculates adjusted and Bayesian Information Criterion for \code{nmm} object
#' 
#' @import maxLik
#' @param object Fitted \code{nmm} model. 
#' @param ... Not used.
#' @param k Multiplication factor.
#' @rdname AICc
#' @return a numeric value with the corresponding AIC, AICc, BIC.
#' @examples
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
#' res <- nmm(ppine, eq_c=eq_c, start_v=start.values, par_c=parl,
#' eq_type = "cont", best_method = FALSE)
#' aa <- in2nmm(res, model.sur$b)
#' AICc(res)
#' AICc(aa)
#' AIC(res)
#' AIC(aa)
#' BIC(res)
#' BIC(aa)
#' @export
AICc <- function(object, ..., k = 2)
{
  UseMethod("AICc", object)
}


#' @rdname AICc
#' @method AICc nmm
#' @S3method AICc nmm
#' @export
AICc.nmm <- function(object, ..., k = 2)
{
  np <- nrParam(object, free=TRUE)
  nn <- nrow(attributes(object)$data)
  AIC(object) - (k * np^2+k*np)/(nn-np-1)
}


#' @rdname AICc
#' @method AICc default
#' @S3method AICc default
#' @export
AICc.default <-  function(object, ..., k = 2)
{
  print(paste0("Method AICc not implemented for class ", class(object)[1]))
}

#' @rdname AICc
#' @method BIC nmm
#' @S3method BIC nmm
#' @export
BIC.nmm <- function(object, ..., k = 2)
{
  np <- nrParam(object, free = TRUE)
  nn <- nrow(attributes(object)$data)
  log(nn)*np - k*logLik(object)
}


#' Bread for Sandwiches.
#' @keywords internal
#' @param object Fitted \code{nmm} model. 
#' @param type Type of bread: "normal", "robust", "clustered".
#' @param new_coef New coefficients.
bread.nmm <- function(object, type = "normal", new_coef = NULL)
{
  if(type!="normal"){
    sde <- meat.nmm(object, type = type, new_coef = new_coef)
  }else{
    sde <- NA
  }
  hv <- if(is.null(new_coef)) object$hessian else  hessian.nmm(object, new_coef=new_coef)
  breadd <- tryCatch({
    qr.solve(-hv)},
    error = function(err){
      r <- solve(-hv)
      return(r)})
  breadd <- if(type!="normal") t(breadd) %*% sde %*% breadd else breadd
  breadd
}

#' Meat for Sandwiches
#' @keywords internal
#' @param object Fitted \code{nmm} model. 
#' @param type Type of bread: "normal", "robust", "clustered".
#' @param new_coef New coefficients.
meat.nmm <- function(object, type="normal", new_coef=NULL)
{
  id <- NULL
  ids <- NULL
  data <- attributes(object)$data
  gv <- gradient.nmm(object, new_coef = new_coef, ind_grad = TRUE)
  k <- NCOL(gv)
  n <- NROW(gv)
  sde <- NULL
  if(type=="robust"){
    #robust meat
    rmeat <- gv%>%as.matrix%>%crossprod
    sde <- rmeat
  }
  if(type=="clustered"){
    region <- data.table(region=data[, "PeID"], ids=1:nrow(data))
    rd <- data.table(region=unique(region[, region]), id=1:length(unique(region[, region])),
                     check.names = TRUE)

    region <- merge(region, rd, by="region", all=TRUE)
    region <- region[order(ids),]
    region[, ids:=NULL]
    reg <- region[, id]

    #clustered meat
    ind <- matrix(0, nrow=nrow(region), ncol=length(unique(region[, id])))
    ind[cbind(seq_along(region[,id]), region[,id])] <- 1
    cgv <- apply(gv, 2, function(x)matrix(x, nrow=1)%*%ind)
    cmeat <- cgv%>%as.matrix%>%crossprod
    sde <- cmeat
  }
  sde
}


#' Calculate Heat rate
#' @keywords internal
#' @param x Fitted \code{nmm} model. 
#' @param prob_func function to calculate probabilities.
hitRate <- function(x, prob_func)
{
  . <- NULL
  value <- NULL
  chc <- NULL
  est <- NULL
  err <- NULL
  PeID <- NULL
  eq <- NULL
  data <- attributes(x)$data
  par <- x$estimate

  probs <- prob_func(par)
  f2 <- function(x){
    indmax <- which.max(x)
    if(length(indmax)>1) browser()
    x[x < indmax] <- 0
    x[x = indmax] <- 1
    x
  }
  probst <- apply(probs, 1, f2)  %>%t
  colnames(probst) <- paste0('est_', 1:ncol(probst))
  rez <- cbind(data[, c('PeID', 'WeID', grep('chc', names(data), value=TRUE))], probst)
  rez %<>% data.table(., stringsAsFactors = FALSE, key = c('PeID', 'WeID'))
  rezt <- rez%>%gather(.,var, value,  3:ncol(.))%>%separate(., var, into=c('var', 'eq'),sep='_')%>%
    spread(., var, value)%>%dplyr::mutate(., err=chc==est)%>%data.table(., key=c('PeID', 'WeID'))
  #total
  discr <- rezt[, sum(!err)/length(PeID)*100]
  discr <- data.table(stat='MisClas', eq="Total",V1=discr, stringsAsFactors = FALSE)
  discr_mode <- rezt[, sum(!err)/length(PeID)*100, by=eq]
  discr_mode <- data.table(stat='MisClas', discr_mode, stringsAsFactors = FALSE)
  discr <- rbind(discr, discr_mode)
  setnames(discr, "V1", "value")
  list(stat=discr, probs=probs)
}

#' pseudo R^2
#' 
#' Calculates pseudo R^2 for discrete choice part.
#' 
#' @param x Fitted \code{nmm} model. 
#' @param which which pseudo R^2 to calculate, options are:
#'  "all", "McFadden", "adjMcFadden", "Cox&Snell","Nagelkerke".
#' @param only_total If \code{TRUE}, compute R^2 only for the whole sample.
#' @return matrix with goodness of fit measures
#' @examples
#' library(mlogit)
#' data("Fishing", package = "mlogit")
#' Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")
#' ## a pure "conditional" model
#' mres <- summary(mlogit(mode ~ price + catch, data = Fish))
#' data <- prepare_data(Fish %>% data.frame %>% dplyr::select(-idx),
#' choice="alt", dummy="mode", PeID="chid", mode_spec_var = c("price", "catch"),
#'type="long")
#' eq_d <- c("a1 + p1 * price_1 + p2 * catch_2", "a2 + p1 * price_2 + p2 * catch_2",
#'              "a3 + p1 * price_3 + p2 * catch_3", "a4 + p1 * price_4 + p2 * catch_4")
#' par_d <- c(paste0("a", 1:4), paste0("p", 1:2))
#' res <- nmm(data, eq_d=eq_d, par_d = par_d,  eq_type="disc", fixed_term=FALSE,
#' best_method=FALSE)
#' pseudoR(res, which = c("McFadden"))
#' ncf <- c(mres$coefficients)
#' names(ncf) <- par_d[-1]
#' mress <- in2nmm(res, new_coef =  ncf)
#' pseudoR(mress, which = c("McFadden"))
#' pseudoR(mress)
#' @export
pseudoR <- function(x, which = c("all", "McFadden", "adjMcFadden", "Cox&Snell", "Nagelkerke"), 
                     only_total = FALSE)
{
  . <- NULL
  eq <- NULL
  #https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/
  try((ls()[ls()%in%c("mcfad", "amcfad", "cox", "nagel")]) %>%
        paste0("'", ., "'")%>% paste0(., collapse=", ") %>%
        paste0("rm(list=c(", ., "))") %>% parse(text=.) %>% eval, silent=TRUE)
  # statistics
  h1r <-  attributes(x)$LL[(length( attributes(x)$eq_c)+1):length( attributes(x)$LL)]
  neqd <- length(h1r)
  
  freq <- x$freq
  h2r <- freq * log(prop.table(freq))
  
  Lh1r <- exp(h1r)
  Lh2r <- exp(h2r)
  nn <- attributes(x)$data %>% nrow
  para_cont <- attributes(x)$par_d
  neq_par <- length(para_cont)
  
  rrd <- data.table(stat="XXX", eq=c("Total", 1:neqd), value=NA)
  
  if(any(which%in%c("McFadden", "all"))){
    # McFaddedn R2sq
    tmcr2 <- 1 - sum(h1r)/sum(h2r)
    # separate, can be negative
    stmcr2 <- 1 - h1r/h2r
    mcfad <- copy(rrd)
    mcfad$stat <- "McFadden"
    mcfad$value <- c(tmcr2, stmcr2)
  }
  if(any(which%in%c("adjMcFadden", "all"))){
    # McFadden adjusted
    tmcr2 <- 1 - (sum(h1r)-sum(neq_par))/sum(h2r)
    # separate, can be negative
    stmcr2 <- 1 - (h1r-neq_par)/h2r
    amcfad <- copy(rrd)
    amcfad$stat <- "adjMcFadden"
    amcfad$value <- c(tmcr2, stmcr2)
  }
  if(any(which%in%c("Cox&Snell", "all"))){
    #Cox & Snell
    tmcr2 <- 1 - (sum(Lh2r)/sum(Lh1r))^(2/nn)
    # separate, can be negative
    stmcr2 <- 1 - (Lh2r/Lh1r)^(2/nn)
    cox <- copy(rrd)
    cox$stat <- "Cox&Snell"
    cox$value <- c(tmcr2, stmcr2)
  }
  if(any(which%in%c("Nagelkerke", "all"))){
    #Nagelkerke / Cragg & Uhler’s
    tmcr2 <- (1 - (sum(Lh2r)/sum(Lh1r))^(2/nn))/(1-sum(Lh2r)^(2/nn))
    # separate, can be negative
    stmcr2 <- (1 - (Lh2r/Lh1r)^(2/nn))/(1-Lh2r^(2/nn))
    nagel <- copy(rrd)
    nagel$stat <- "Nagelkerke"
    nagel$value <- c(tmcr2, stmcr2)
  }
  stats <- (ls()[ls()%in%c("mcfad", "amcfad", "cox", "nagel")])  %>%paste0(., collapse=", ") %>%
    paste0("rbind(", ., ")") %>% parse(text=.) %>% eval
  if(only_total) stats <- stats[eq=="Total", ]
  stats
}


#' Goodness of fit measures
#' 
#' Calculate RMSE, MAPE, R^2 and adjusted R^2
#' 
#' @param x Fitted \code{nmm} model.
#' @param which What to calculate. Options: "all", "RMSE", "MAPE", "Rx2", "Rx2adj".
#' @param only_total If \code{TRUE}, calculate statistics only for totals.
#' @return matrix with Goodness of fit measures
#' @examples
#' library(systemfit)
#' data(ppine , package="systemfit")
#' hg.formula <- hg ~ exp( h0 + h1*log(tht) + h2*tht^2 + h3*elev)
#' dg.formula <- dg ~ exp( d0 + d1*log(dbh) + d2*hg + d3*cr)
#' labels <- list( "height.growth", "diameter.growth" )
#' model <- list( hg.formula, dg.formula )
#' start.values <- c(h0=-0.5, h1=0.5, h2=-0.001, h3=0.0001,
#'                   d0=-0.5, d1=0.009, d2=0.25, d3=0.005)
#' model.sur <- nlsystemfit( "SUR", model, start.values, data=ppine, eqnlabels=labels )
#' eq_c <- as.character(c(hg.formula, dg.formula))
#' parl <- c(paste0("h", 0:3),paste0("d", 0:3))
#' res <- nmm(ppine, eq_c=eq_c, start_v=start.values, par_c=parl,
#' eq_type = "cont", best_method = FALSE)
#' cont_stats(res, which = "all")
#' @export
cont_stats <- function(x, which=c("all", "RMSE", "MAPE", "Rx2", "Rx2adj"), only_total = FALSE)
{
  . <- NULL
  rezct <- NULL
  est <- NULL
  eqi <- NULL
  neq_par <- NULL
  kpar <- NULL
  cont_eq <- NULL
  V1 <- NULL
  stat <- NULL
  eq <- NULL
  try((ls()[ls()%in%c("contr", "contrm", "rsq", "rsqa")]) %>% paste0("'", ., "'") %>%
        paste0(., collapse=", ") %>%
        paste0("rm(list=c(", ., "))") %>% parse(text=.) %>% eval, silent=TRUE)
  help_diagnostics(x, env = environment())

  ##ssr
  ssrt <- data.table(stat="ssr", eq="Total", V1=rezct[, (est-eqi)^2] %>% sum, stringsAsFactors = FALSE)
  #by equation:
  ssr_mode <- rezct[, sum((est-eqi)^2), by=eq]
  ssr_mode <- data.table(stat='ssr', ssr_mode, stringsAsFactors = FALSE)
  ssr <- rbind(ssrt, ssr_mode)
  
  if(any(which%in%c("RMSE", "all"))){
    #total
    contr <- data.table(stat="RMSE", eq="Total", V1=rezct[, (est-eqi)^2] %>% sum %>% "/"(nrow(rezct)-sum(neq_par)) %>%
                          sqrt, stringsAsFactors = FALSE)
    #by equation:
    contr_mode <- rezct[, sqrt((sum((est-eqi)^2))/(length(eqi)-unique(kpar))), by=eq]
    contr_mode <- data.table(stat='RMSE', contr_mode, stringsAsFactors = FALSE)
    contr <- rbind(contr, contr_mode)
  }

  if(any(which%in%c("MAPE", "all"))){
    #total
    contrm <- data.table(stat="MAPE", eq="Total", V1=rezct[eqi!=0, abs((eqi-est)/eqi)]%>%
                           mean  %>% "*"(100) ,
                         stringsAsFactors = FALSE)
    #by equation:
    contrm_mode <- rezct[eqi!=0, mean(abs((eqi-est)/eqi))*100, by=eq]
    contrm_mode <- data.table(stat='MAPE', contrm_mode, stringsAsFactors = FALSE)
    contrm <- rbind(contrm, contrm_mode)
  }


  if(any(which%in%c("Rx2", "Rx2adj", "all"))){
    rsq <- data.table(stat="Rx2", eq="Total", 1 - rezct[, sum((est-eqi)^2)/(crossprod(eqi) - mean(eqi)^2 * length(eqi))],
                      stringsAsFactors = FALSE)
    rsq_mode <- rezct[, 1 - sum((est-eqi)^2)/(crossprod(eqi) - mean(eqi)^2 * length(eqi)), by = eq]  
    rsq_mode <- data.table(stat='Rx2', rsq_mode, stringsAsFactors = FALSE)
    rsq <- rbind(rsq, rsq_mode)
  }

  if(any(which%in%c("Rx2adj", "all"))){
    nn <- c(nrow(rezct), rep(nrow(rezct)/length(cont_eq), length(cont_eq)))
    rsqa <- copy(rsq)
    rsqa[, V1:=1-(1-rsq$V1)*((nn-1)/(nn-c(sum(neq_par), neq_par)))]
    rsqa[, stat:="Rx2adj"]
  }

  stats <- (ls()[ls()%in%c("contr", "contrm", "rsq", "rsqa")]) %>% paste0(., collapse=", ") %>%
    paste0("rbind(", ., ")") %>% parse(text=.) %>% eval
  if(only_total) stats <- stats[eq=="Total", ]
  stats
}

#' Generates objects needed for diagnostics functions
#' @param x Fitted \code{nmm} model.
#' @param env To which environment to write the results.
#' @keywords internal
help_diagnostics <- function(x, env=parent.env())
{
  . <- NULL
  value <- NULL
  data <- attributes(x)$data
  cont_eq <- attributes(x)$eq_c
  disc_eq <- attributes(x)$eq_d
  cont_e <- attributes(x)$cont_e #with par[[]]
  disc_e <- attributes(x)$disc_e
  eq_type <- attributes(x)$type
  par_c <- attributes(x)$par_c
  par_d <- attributes(x)$par_d
  
  if(!is.null(disc_eq)&(!is.null(cont_eq))){
    npar <- length(par_c)
    disc_e <- gsubfn('(par\\[)(\\d{1,10})(\\])', function(x, y, z)paste0(x,
                                                                        as.numeric(y) + npar, z), disc_e)
  }
  cd_eqs <- stats_function(cont_eq, disc_eq, par_c = par_c, par_d = par_d)

  methodopt <- "BHHH"
  check_hess <- TRUE
  separatenmm <- TRUE
  namesf <- switch(eq_type, disc="prob_func", cont="cont_func",
                   joint=c("prob_func", "cont_func"))
  cd_eqs <- switch(eq_type, disc=cd_eqs["prob"], cont=cd_eqs["cont"], joint=cd_eqs[-3])

  sapply(1:length(cd_eqs), function(x){
    xx <-  cd_eqs[[x]]
    environment(xx) <- parent.env(environment())
    assign(namesf[x], xx, envir = parent.env(environment()))
    environment(xx) <- env
    assign(namesf[x], xx, envir = env)
  })
  #env <- parent.env(environment())
  eval(parse(text=paste0("environment(", namesf[1:length(cd_eqs)], ") <- env")))

  if(eq_type!="disc"){
    para_cont <- get_par(par_c, cont_eq)
    neq_par <- para_cont$neq_par %>% unlist
    endog <- para_cont[["endog"]]
    exog <- para_cont[["exog"]]
    data <- data[, intersect(names(data), c('PeID', endog, exog))]%>%unique
    environment(cont_func) <- environment()
    conts <- cont_func(x$estimate)
    colnames(conts) <- paste0('est_', 1:ncol(conts))
    rezc <- cbind(data[, c('PeID',  endog)], conts)
    names(rezc)[names(rezc)%in%endog] <- paste0('eqi_', 1:length(endog))
    rezc %<>% data.table(., stringsAsFactors = FALSE, key = c('PeID'))%>%unique(., by=c('PeID'))
    rezct <- rezc %>% gather(.,var, value,  2:ncol(.)) %>%
      separate(., var, into=c('var', 'eq'),sep='_') %>%
      spread(., var, value)  %>% data.table(., key=c('PeID'))# %>% dplyr::filter(., eqi!=0) why this needed?
    # add number of parameters in each equation
    npardat <- data.table(eq=(1:length(unique(rezct$eq))) %>% as.character, kpar=neq_par)
    rezct <- merge(rezct, npardat, by="eq")
  }

  obj_exp <- ls(envir = environment())
  obj_exp <- obj_exp[!obj_exp%in%c("x", "cont_func", "cd_eqs", "env", "obj_exp", "y")]
  for(y in 1:length(obj_exp)){
    xx <- obj_exp[[y]]
    assign(xx, eval(parse(text=xx)), envir = env)
  }
}

#' Goodness of fit measures for both parts
#' 
#' Calculation RMSE, misclassification and other goodness of fit measures.
#' 
#' @param x Fitted \code{nmm} model.
#' @param xdigit rounding number
#' @param which What to calculate. Options: "all", "RMSE", "MAPE", "Rx2", "Rx2adj".
#' @param only_total If \code{TRUE}, calculate statistics only for the whole system
#' @param cPseudoR If \code{TRUE}, calculate pseudo R^2s.
#' @param cRs Include "AIC", "AICc", "BIC"
#' @return matrix with goodness of fit measures. 
#' \code{attribute} \code{corr} holds empirical variance-covariance matrix.
#' @examples
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
#' start.values <- c(h0=-0.5, h1=0.5, h2=-0.001, h3=0.0001,
#'                   d0=-0.5, d1=0.009, d2=0.25, d3=0.005)
#' res <- nmm(ppine, eq_c=eq_c, start_v=start.values, par_c=parl, eq_type = "cont",
#' best_method = FALSE)
#' ressur <- in2nmm(res, new_coef=model.sur$b)
#' diagnostics(res)
#' diagnostics(ressur)
#'
#' #example discrete
#' library(mlogit)
#' data("Fishing", package = "mlogit")
#' Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")
#' ## a pure "conditional" model
#' mres <- summary(mlogit(mode ~ price + catch, data = Fish))
#' data <- prepare_data(Fish %>% data.frame %>% dplyr::select(-idx), 
#' choice="alt", dummy="mode", PeID="chid", mode_spec_var = c("price", "catch"),
#' type="long")
#' eq_d <- c("a1 + p1 * price_1 + p2 * catch_2", "a2 + p1 * price_2 + p2 * catch_2",
#'              "a3 + p1 * price_3 + p2 * catch_3", "a4 + p1 * price_4 + p2 * catch_4")
#' par_d <- c(paste0("a", 1:4), paste0("p", 1:2))
#' res <- nmm(data, eq_d=eq_d, par_d=par_d,  eq_type="disc")
#' ncoef <- mres$coefficients
#' names(ncoef) <- par_d[-1]
#' resdisc <- in2nmm(res, new_coef = ncoef)
#' a <- diagnostics(res, xdigit=2)
#' a2 <- diagnostics(resdisc)
#' attributes(a2)$corr
#' @export
diagnostics <- function(x, xdigit = 4, which = "all", only_total = FALSE,
                        cPseudoR = TRUE, cRs = TRUE)
{
  . <- NULL
  eq_type <- NULL
  disc_eq <- NULL
  prob_func <- NULL
  rezct <- NULL
  eqi <- NULL
  est <- NULL
  V1 <- NULL
  endog <- NULL
  eq <- NULL
  value <- NULL
  par_c <- attributes(x)$par_c
  par_d <- attributes(x)$par_d
  data <- attributes(x)$data
  help_diagnostics(x, env = environment())
  
  data <- attributes(x)$data
  if(eq_type!="cont"){
    # discrete part
    para_cont <- get_par(par_d, disc_eq)
    neq_par <- para_cont$neq_par %>% unlist
    discrl <- hitRate(x, prob_func)
    discr <- discrl$stat
    probs <- discrl$probs

    rrd <- discr
    probsd <- data.table(data[, c("WeID", "PeID")],probs)
    probsd1 <- probsd[, lapply(.SD, function(x)qnorm(ifelse(x==0, 1e-100,
                                                            ifelse(x==1, 1-1e-16,x )))),
                      by=c('PeID', 'WeID'), .SDcols=grep('res',names(probsd), value=TRUE)]
    nn <- nrow(probsd1)
    nchoice <- 1:length(disc_eq)
    cnamesd <- if(eq_type=="disc")paste0('P(',
                                    nchoice[-length(nchoice)],')') else paste0('P(',nchoice,')')
    rnamesd <- paste0('P(', nchoice[-1],')')
    if(cPseudoR){
      psR2 <- pseudoR(x, which=c("McFadden", "adjMcFadden"), only_total = only_total)
    }else{
      psR2 <- NULL
    }
    rrd <- rbind(rrd, psR2)
    rrc <- NULL
  }else{
    rrd <- NULL
    probsd <- probsd1 <- data.table(PeID=data[, "PeID"])
    nchoice <- NA
    cnamesd <- NULL
    rnamesd <- NULL
  }

  if(eq_type!="disc"){
    #cont par
    rrc <- cont_stats(x, which=which, only_total = only_total)

    rezcc <- rezct[, eqi-est, by=c("PeID",  "eq")]%>%spread(., 'eq', V1)
    setnames(rezcc, c("PeID", paste0('est',1:length(endog) )))
    cnamesc <- endog[-length(endog)]
    rnamesc <-  paste0('e(',if(eq_type=="cont") endog[-1] else endog,')')
    setnames(rrc, "V1", "value")
  }else{
    rrc <- NULL
    rezcc <- data.table(PeID=data[, "PeID"])
    cnamesc <- NULL
    rnamesc <- NULL
  }

  if(cRs){
    rrt <- data.table(stat=c("AIC", "AICc", "BIC" ), eq="Total", value=c(AIC(x),
                                                                         AICc(x), BIC(x)))
  }else{
    rrt <- NULL
  }

  rr <- rbind(rrc, rrd, rrt)
  ###correlations:
  
  if(eq_type!="disc"){
    #cont par
    corrd <- merge(probsd, rezcc,
                   all=TRUE, by=c("PeID"))
    corrd1 <- merge(probsd1, rezcc, all=TRUE, by=c("PeID"))
  }else{
    corrd <- probsd 
    corrd1 <- probsd1
  }
  

  #correlation with stat. significance
  a1 <- corrd[,  c(grep('res|est', names(corrd), value=TRUE)), with=FALSE]%>%
    as.matrix%>%corstarsl(., xdigit)
  a1 <- a1[-1, ]

  a11 <- corrd1[,  c(grep('res|est', names(corrd1), value=TRUE)), with=FALSE]%>%
    as.matrix%>%corstarsl(., xdigit)
  a11 <- a11[-1,]
  if(is.null(dim(a1))){
    a1 <- data.frame(a1)
    a11 <- data.frame(a11)
  }

  names(a11) <- names(a1) <- c(cnamesd, cnamesc)
  rownames(a11)<- rownames(a1) <- c(rnamesd, rnamesc)
  a1 %<>% data.table(eq_type = eq_type, probt = "normal", var = row.names(.), .)
  a11 %<>% data.table(eq_type = eq_type, probt = "quantile", var = row.names(.), .)
  a1 <- rbind(a1, a11)
  
  if(only_total) rr <- rr[eq=="Total", ]
  rr[, value:=round(value, xdigit)]
  attributes(rr) <- c(attributes(rr), list(corr=a1, errors=corrd1))
  return(rr)
}


#' Significance of correlation matrix
#' @importFrom Hmisc rcorr
#' @keywords internal
#' @param x correlation matrix
corstarsl <- function(x, xdigit=2){
  . <- NULL
  x <- as.matrix(x)
  R <- rcorr(x)$r
  p <- rcorr(x)$P

  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))

  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), xdigit))[,-1]

  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")

  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)

  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew)
}


#' Takes out only data used in continuous estimation
#' @keywords internal
#' @param data full data matrix
#' @param cont_eq cont. equations
#' @param par_c parameter vector of cont. equations
extract_cont_data <- function(data, cont_eq, par_c)
{
  . <- NULL
  data <- data[, intersect(names(data),
                           c('PeID', get_par(par_c, cont_eq)[c("exog", "endog")] %>% unlist %>%
                               unique))] %>% unique
  data
}



#' @keywords internal
summary_check <- function (object)
{
  . <- NULL
  result <- object$maxim
  coef2obj <- object$estimate
  nParam <- length(coef2obj)
  activePar <- activePar(object)
  # variance covariance
  if (!is.null(object$varcovar)){
    vcov2 <- object$varcovar
  }else{
    vcov2 <- if (!is.null(hess <- hessian(object))) {
      hess <- hessian(object)[activePar, activePar, drop = FALSE]
      hessev <- abs(eigen(hess, symmetric = TRUE, only.values = TRUE)$values)
      varcovar <- matrix(0, nrParam(object), nrParam(object))
      rownames(varcovar) <- colnames(varcovar) <- names(coef2obj)
      varcovar[activePar, activePar] <-  tryCatch({
        qr.solve(-hessian(object)[activePar, activePar])},
        error = function(err){
          r <- solve(-hessian(object)[activePar, activePar])
          return(r)})
      varcovar <- (varcovar + t(varcovar))/2
      varcovar
    } else NULL
  }
  # std error
  stdEr2 <- if (!is.null(vc <- vcov2)) {
    s <- sqrt(diag(vc))
    names(s) <- names(coef2obj)
    s
  } else NULL
  
  if ((object$code < 100) & !is.null(coef2obj)) {
    t <- coef2obj/stdEr2
    p <- 2 * pnorm(-abs(t))
    t[!activePar(object)] <- NA
    p[!activePar(object)] <- NA
    results <- cbind(Estimate = coef2obj, `Std. error` = stdEr2,
                     `t value` = t, `Pr(> t)` = p)
  } else {
    results <- NULL
  }
  summary <- list(maximType = object$type, iterations = object$iterations,
                  returnCode = object$code, returnMessage = object$message,
                  loglik = object$maximum, estimate = results, fixed = !activePar,
                  NActivePar = sum(activePar), constraints = object$constraints)
  class(summary) <- "summary.maxLik"
  summary
}

#' \code{nrParam} return number of parameters used in the estimation.
#' @param x Fitted \code{nmm} model.
#' @param free Sum or length of parameters.
#' @keywords internal
nrParam <- function (x, free = FALSE) 
{
  if (!inherits(x, c("nmm", "maxLik"))) {
    stop("'nrParam' called on non-'nmm' or 'maxLik' object")
  }
  if (free) sum(activePar(x)) else length(x$estimate)
}

