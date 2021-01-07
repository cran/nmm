#' @keywords internal
insert_obj <- function(inde, obj, insert)
{
  obj <- c(obj[1:inde], insert, obj[(inde+1):length(obj)])
  obj
}

#' \code{LL_joint} Function for joint log-likelihood with correlation between
#'  continuous and discrete equations
#' @importFrom gsubfn gsubfn
#' @param ffor Discrete equations.
#' @param cheqs0 Continuous equations.
#' @param objm Expressions from \code{cond_mean_cov_expr}.
#' @param datan Character string. Name of data-set.
#' @param cfunc Character string. Name for continuous block function.
#' @param cgrad Character string. Name for continuous block gradient.
#' @param chess Character string. Name for continuous block hessian.
#' @param gfunc Character string. Name for joint gradient.
#' @param hfunc Character string.  Name for joint Hessian.
#' @param check_hess if \code{TRUE}, Hessian is checked.
#' @param print_out if \code{TRUE}, prints out LL for each equation.
#' @param bayesian_random If \code{TRUE}, than par[1] is changed to par[,1] to be used for 
#' optimization of random parameters in Bayesian estimation.
#' @param MNtypef "logit" or "dogit".
#' @return Function for joint estimation
#' @keywords internal
LL_joint <- function(ffor, cheqs0, objm = NULL, datan = "data", cfunc = "cfunc", cgrad = "cgrad",
                     chess = "chess", gfunc = "gfunc", hfunc = "hfunc", check_hess = FALSE, 
                     print_out = FALSE, bayesian_random = FALSE, MNtypef = "logit")
{
  . <- NULL
  if(is.null(objm )){
    objm <- cond_mean_cov_expr(neqt=length(cheqs0), kk=length(ffor))
  }
  hessian <- if(check_hess==TRUE) hfunc else NULL
  npar <- get_npar(cheqs0) # number of prameters in cont. equations
  
  #in cheqs0 par numeration begins by 1, we need to shift it to begin from #params from cont. eq.
  fforn <- gsubfn('(par\\[)(\\d{1,10})(\\])', function(x, y, z)paste0(x,
                                                                     as.numeric(y) + npar, z), ffor)
  parlc <- get_npar(c(cheqs0, fforn), values=TRUE)
  
  #parameters
  indexp <- gsub("par\\[|\\]", '', parlc)
  parlc <- parlc[indexp%>%as.numeric%>%order]
  parln <- parlc%>%gsub('\\[|\\]', '', .)
  partable <- data.frame(parlp=parln, parnl=parlc, stringsAsFactors = FALSE, row.names = NULL)

  neqt <- length(cheqs0) #number of cont equations
  kk <- length(ffor) # number of disc equations
  fixed <- npar+1 #which parameter is fixed to be 0

  #expressions for conditional mean and cov
  meannewel <- objm$meannewel #mean
  denoml <- objm$denoml # denominator in mean
  nominl <- objm$nominl # nominator in mean
  cel <- objm$cel # variance

  #discrete part
  if(MNtypef=="logit"){
    mnf1 <-  MNlogitf(fforn, transform = TRUE, separatenmm = TRUE, weight_disc = TRUE)
    ppardogit <- NULL
  }else{
    ndd <- fforn %>% gsub("(\\[\\d{1,10}\\])", "XX\\1XX", .) %>% strsplit("XX") %>% unlist %>% unique %>% 
      grep("\\[", ., value = TRUE) %>% gsub("\\[|\\]", "", .) %>% as.numeric() %>% max
    ppardogit <- c(paste0('par[',(ndd+1):(ndd+length(fforn)),']'))
    print(ppardogit)
    mnf1 <- MNdogitf(fforn, transform = TRUE, separatenmm = TRUE,
                     weight_disc = TRUE, ppardogit)
    px111 <- data.frame(parlp=ppardogit %>% gsub("\\[|\\]", "", .), 
                        parnl=ppardogit, stringsAsFactors = FALSE, row.names = NULL)
    partable <- rbind(partable, px111)
  }
  
  #substitute equations into the expression of mean
  mele <- meannewel%>%gsub("meanqt[[:digit:]] <- ", "", .)
  #substitute eta expressions
  cheqs0p <- gsub('\\[|\\]', '', cheqs0)
  repdat <- data.frame(old=paste0('eta\\[\\[', 1:neqt,'\\]\\]'), new=paste0('(',cheqs0,')'),
                       stringsAsFactors = FALSE)
  meleta <- replace_par_wrap(repdat, mele)
  meleta%<>%gsub('(par\\[)([[:digit:]]{1,10})(\\])', 'par\\2',.)
  meleta%<>%gsub('\\^', '**',.)
  meleta%<>%gsub('\\[\\[', '_',.)%>%gsub('\\]\\]','_', .)%>%gsub(",", "x", .)%>%
    gsub('(nomin\\d)\\[(\\d)\\]', '\\1_\\2_',.)
  #add expressions of cond. means
  mnf1[["meanexpr"]] <- meleta
  
  
  #calc derivatives
  #replace pedef with derivative of conditional mean to get the second part of gradient expression
  expr <- deriv_disc_expr(mnf1, smodu1=partable, attribute='hessian', ffor=fforn,
                          MNtypef = MNtypef)
  
  expr%<>%gsub("par\\[\\[(\\d{1,10})\\]\\]", "par[\\1]",.)
  ple <- paste0('par <- c(par[1:(',fixed,'-1)], 0,par[(',fixed,'):length(par)] )')
  expr <- c('{',ple, expr[-c(1,2)])

  #replace cov and mean expressions in expr
  expr[grepl("covqt([[:digit:]]{1,10}) <\\-", expr)] <- cel
  expr[grepl("meanqt([[:digit:]]{1,10}) <\\-", expr)] <- meannewel
  #add expressions for denominator and nominator
  inde <- expr%>%grep("^meanqt", .)%>%min-1
  iobj <- c(paste0('nomin', 1:kk, ' <- c()'),
            denoml, nominl)
  expr <- insert_obj(inde, expr, iobj)
  #transformed probabilities
  inde <- expr%>%grep("^ife\\[\\[([[:digit:]]{1,10})\\]\\] <\\-", .)%>%max
  transf <- paste0('transf[[', 1:kk, ']] <- qnorm(ife[[', 1:kk, ']]) * avl_', 1:kk, ' * chc_', 1:kk)
  #add transf and err from cont =etas
  etae <- paste0('eta[[',1:neqt, ']] <- ', cheqs0)
  m1 <- paste0('merged <- data.table(',datan,'[, c("PeID", "WeID")],',1:neqt%>%
          paste0('eta', ., ' = eta[[', .,']]')%>%paste0(., collapse=', '),',
               ',1:kk%>%paste0('transf', ., ' = transf[[', .,']]')%>%paste0(., collapse=', '),')')
  m2 <- paste0('merged <- merged[, lapply(.SD, function(x) mean(x, na.rm=FALSE)), by=PeID, .SDcols=c(paste0("eta", 1:',
               neqt,'),paste0("transf", 1:',kk,'))]')
  m2 <- "merged <- merged"
  rhoe <- paste0('rho',1:kk,' <- rho[inde[[',1:kk,']], inde[[',1:kk,']]]')
  sigmae <- paste0('sigma',1:kk,' <- sqrt(diag(sigma[inde[[',1:kk,']], inde[[',1:kk,']]]))')

  iobj <- c('transf <- c()', transf , 'eta <- c()',
          #error from cont eq.
          etae, m1, m2,
          'sigma <- cov(merged[,-c(1,2), with=FALSE],use="pairwise.complete.obs")',
          'rho <- cor(merged[,-c(1,2), with=FALSE],use="pairwise.complete.obs")',
          paste0('inde <- lapply(1:',kk,', function(x) c(1:',neqt,', ',neqt,'+x))'),
          rhoe, sigmae)
  expr <- insert_obj(inde, expr, iobj)

  #meanqt now should be a list
  expr%<>%gsub('meanqt <- c', 'meanqt <- list',.)
  expr%<>%gsub('smodu1', 'partable', .)

  #add gradient values of cont. eq
  inde <- min(grep("derive\\d <-", expr))-1
  iobj <- c('methodo <- methodopt','methodopt <- "BHHH"',
            'separatenmm <- FALSE',
            paste0('try(environment(',cgrad,') <- environment(), silent = TRUE)'),
            #paste0('gradcont <- (',cgrad,'(par)*(1/npaths_cont))[,-',(npar+1),']'),
            paste0('gradcont <- (',cgrad,'(par)*(1/npaths)*(1/npaths_cont))[,-',(npar+1),']'),
            paste0('if(!methodo%in%c("BHHH")) gradcont <- apply(gradcont,  2, function(x)sum(x, na.rm=FALSE))'),
            'methodopt <- methodo')
  expr <- insert_obj(inde, expr, iobj)

  #add hessian of cont
  inde <- min(grep("hessian\\d <-", expr))-1
  iobj <- c('methodo <- methodopt','methodopt <- "BHHH"',
            'separatenmm <- FALSE',
            paste0('try(environment(',hfunc,') <- environment(), silent = TRUE)'),
            paste0('try(environment(',chess,') <- environment(), silent = TRUE)'),
            paste0('hesscont <- ',chess,'(par)*(1/npaths)*(1/npaths_cont)'),
            #paste0('hesscont <- ',chess,'(par)*(1/npaths)*(1/npaths_cont)'),
          'methodopt <- methodo')
  expr <- insert_obj(inde, expr, iobj)
  expr %<>%  gsub("sigma\\[(\\d)\\]", "sqrt(sigma[\\1,\\1])", .)

  
  #function expressions
  funcex <- expr[!grepl(paste0('(((cd)|(bd\\d_\\d{1,10})|(bd\\d)|(gd\\d_\\d{1,10})|(gd\\d)|(edd\\d_\\d{1,10})|(delta)|(edd\\d)|(ad)|(dd)|(gd)|(ed)(r\\d{1,10})|(hesscont)|(gradcont)|(meletah)|(meanhe\\d{1,10})) <-)|(hessian)|(',
                               cgrad,')|(',chess,')|(',gfunc,')|(pedef)|(meange)|(derive)|(term)'),expr)]
  
  #add cont part
  funcex <- funcex[!funcex%in%c("","res","}","{")]
  funcex <- c( funcex,
               'methodo <- methodopt',
             paste0('try(environment(',cfunc,') <- environment(), silent = TRUE)'),
             'methodopt <- "BHHH"',
             'separatenmm <- TRUE',
             #paste0('contpart <- ', cfunc,'(par[1:',npar,'])*(1/npaths_cont)'),
             paste0('contpart <- ', cfunc,'(par[1:',npar,'])*(1/npaths)*(1/npaths_cont)'),
             'methodopt <- methodo'
             )
  discp <- 1:kk%>%paste0('avl_',. , ' * chc_', ., ' * wd_',.,' * log(pife[[',., ']])' )%>%
    paste0(., collapse=' + ')%>%paste0('discpart <- (1/npaths)*(', ., ')')
  discp <- 1:kk%>%paste0("discpart", . ," <-  (1/npaths)*(",'avl_',. , ' * chc_', ., ' * wd_',.,' * log(pife[[',., ']]))' )
  discp <- c(discp,paste0("discpart", 1:kk) %>% paste0(., collapse = ", ") %>% paste0('discpart <- cbind(',., ')'),
             "contpart <- cbind(contpart, discpart)")
  
  funcex <- c(funcex, discp,
              'res <- if(methodopt=="BHHH") rowSums(contpart, na.rm=FALSE)  else sum(contpart%>%rowSums, na.rm=FALSE)',
              ifelse(print_out,
                     'cat(paste0("Continuous LL:", sum(contpart, na.rm=FALSE)%>%round(., 5), ", discrete LL:", sum(discpart, na.rm=FALSE)%>%round(., 5)), "\n")', ''),
              'atrb <- c(contpart%>%colSums(., na.rm=FALSE))')
  
  if(!is.null(hessian)){
    t05 <- c(paste0('if(check_hess==TRUE){',' hess <- try(qr.solve(-',hessian, '(par[-',npar+1,
                    '])), silent=TRUE)'), 'if(any(class(hess)%in%"try-error")) hess <- matrix(-0.01,1,1)',
             '}else{', 'hess <- matrix(0.01,1,1)}')
    t15 <- paste0('if(any(diag(hess)<0)) res <- -Inf')
    funcex <- c(funcex, t05, t15)
  }
  funcex <- c( '{', funcex,'attributes(res) <- c(attributes(res), list(LLv=atrb))','res','}')
  funcex <-c('{', paste0('with(',datan, ','),funcex,')','}')
  funcex <- grep("res <- res", funcex, value=TRUE, invert = TRUE)
  #remove unnecessary lines
  indexrm <- funcex%>%grep("^(methodo|methodopt|separatenmm|methodpt) ", .)
  funcex <- funcex[-indexrm[1:4]]
  strm <- funcex%>%grep("methodo <- methodopt", ., fixed = TRUE)
  enrm <- funcex%>%grep("methodopt <- methodo", ., fixed = TRUE)

  if(bayesian_random==TRUE){
    actbr <- .%>%gsub("par\\[\\[(\\d{1,10})\\]\\]", "par[,\\1]", .)%>%gsub("par\\[(\\d{1,10})", "par[,\\1",.)%>%
      gsub("par\\[\\((\\d{1,10})", "par[,(\\1",.)%>%
      gsub("par <- c(", "par <- cbind(", ., fixed=TRUE)%>%
      gsub("length(par)", "ncol(par)", ., fixed=TRUE)
    funcex%<>%actbr
  }
  
  # get text for function
  indes <- grep("par <- ", funcex, fixed = TRUE) 
  indee <- grep("methodo <- methodopt", funcex, fixed = TRUE)[1]
  text_func <- funcex[(indes+1):(indee-1)]
  funcex1 <- c(funcex[1:indes],'eval(parse(text=text_func))', funcex[indee:length(funcex)])
  ffuncf <- function(par){}
  body(ffuncf) <- parse(text=funcex1)
  
  #grad expressions
  grad <- expr[!grepl(paste0('(((r(|z|p)\\d{1,10})|(cd)|(bd\\d_\\d{1,10})|(bd\\d)|(gd\\d_\\d{1,10})|(gd\\d)|(edd\\d_\\d{1,10})|(delta)|(edd\\d)|(ad)|(ad\\d)|(cd\\d)|(dd)|(dd\\d)|(gd)|(ed)(r\\d{1,10})|(hesscont)|(meletah)|(meanhe\\d{1,10})) <-)|(hessian)|(',
                             chess,')|(',hfunc,')|(',cfunc,')'),expr)]
  xg <- paste0('par', 1:nrow(partable))[-fixed]%>%paste0('"', ., '"')%>%
    paste0(., collapse=", ")%>%paste0('c(', ., ')')
  inde <- which(grad=="")
  iobj <- c(paste0(' if(!methodo%in%c("BHHH")) res <- array(derive[-',fixed,
                   '], c(1,',nrow(partable)-1,'), list(NULL, ',xg,'))%>%unlist'),
            paste0(' if(methodo%in%c("BHHH")) res <- array(derive[ , -',fixed,
                   '], c(nrow(data),',nrow(partable)-1,'), list(NULL, ',xg,'))%>%unlist'))
  grad <- insert_obj(inde, grad, iobj)
  grad <- grad[!grepl("res <- res\\[", grad)]
  grad[grad=='res'] <- paste0('gradcont+res')
  grad%<>%gsub('(par|sigma\\d)(\\[\\[)(\\d{1,10})(\\]\\])', '\\1[\\3]',.)
  grad %<>%c('{', paste0('with(',datan, ','),.,')','}')
  if(bayesian_random==TRUE){
    grad%<>%actbr
  }
  # get text for gradient
  indes <- grep(text_func[length(text_func)], grad, fixed = TRUE) %>% max
  indee <- grep("methodo <- methodopt", grad, fixed = TRUE)[1]
  text_grad <- c("term0 <- c()",grad[(indes+1):(indee-1)])
  indes <- grep("par <- ", grad, fixed = TRUE) 
  grad1 <- c(grad[1:indes],'eval(parse(text=text_func))',
               'eval(parse(text=text_grad))',
             grad[indee:length(grad)])
  
  ffuncg <- function(par){}
  body(ffuncg) <- parse(text=grad1)

  #hessian expressions
  hess <- expr[!grepl(paste0('(derive)|(',gfunc,')|(',cgrad,')|(',cfunc,')|(gradcont)'),expr)]
  hess <- hess[!grepl('res <-', hess)]
  inde <- grep('abind', hess)-1

  iobj <- c('hessian <- lapply(1:(length(par)), function(k) lapply(1:(length(par)),function(j) hessian[[k]][[j]]+ hesscont[,k,j]))',
            paste0('hessian <- lapply(1:length(par), function(k) hessian[[k]][-', fixed, '])'),
            paste0('hessian <- hessian[-', fixed, ']'))
  
  hess <- insert_obj(inde, hess, iobj)
  inde <- grep('abind', hess)
  iobj <- c('if (!methodopt %in%c("BHHH")) hessian <- sapply(1:(length(par)-1), function(k) sapply(1:(length(par)-1), function(j) hessian[, k, j] %>% sum(., na.rm=FALSE)))',
            'res <- hessian')
  hess <- insert_obj(inde, hess, iobj)
  hess%<>%gsub('(par|sigma\\d)(\\[\\[)(\\d{1,10})(\\]\\])', '\\1[\\3]',.)
  hess %<>%c('{', paste0('with(',datan, ','),.,')','}')
  if(bayesian_random==TRUE){
    hess%<>%actbr
  }
  act1 <-  .%>% gsub(".*list\\(", "", .) %>% 
    gsub("\\)", "", .) %>% strsplit(., "\\,") %>% unlist %>% 
    paste0('"', ., '"') %>% paste0(., collapse = ", ")

  
  indf <- (1:length(hess))%in%((grep("^meange", hess)+1):
                                   (grep("cd <- list(", hess, fixed=TRUE)-1))
  hess <- hess[!indf]
  
  # get text for hessian
  indes <- grep(text_grad[length(text_grad)], hess, fixed = TRUE) %>% max
  indff <- grep("methodo <- methodopt", hess, fixed = TRUE)
  if(length(indff)>1){
    indee <- indff[2]
  }else{
    indee <- indff
  }

  text_hess <- hess[(indes+1):(indee-1)]
  text_edd <- text_gd <- text_bd <- NULL
  
  i_bd <- grep("^bd\\d|^bd_b\\d", text_hess)
  text_bd <- text_hess[i_bd]
  text_hess <- c(text_hess[1:(i_bd[1] -1)], 'eval(parse(text=text_bd))',text_hess[(i_bd[length(i_bd)]+1):length(text_hess)])
  
 
  indes <- grep("par <- ", hess, fixed = TRUE) 
  hess1 <- c(hess[1:indes],'eval(parse(text=text_func))',
             'eval(parse(text=text_grad))',
             'eval(parse(text=text_hess))',
             hess[indee:length(hess)])
  i_h <- grep('rm(list=c("hessian', hess1, fixed = TRUE)
  # remove unnecsarry objects
  i_h <- grep('hessian <-\\s{1,3}list', hess1, fixed = FALSE)-1
  hess1 <- c(hess1[1:i_h],
               paste0("rm(list=c(",paste0(paste0('"','term1', 1:kk,'"'), collapse=', ') ,"))"),
               paste0("rm(list=c(",paste0(paste0('"','term2', 1:kk,'"'), collapse=', ') ,"))"),
               paste0("rm(list=c(",paste0(paste0('"','term3', 1:kk,'"'), collapse=', ') ,"))"),
               paste0("rm(list=c(",paste0(paste0('"','meanqt', 1:kk,'"'), collapse=', ') ,"))"),
               paste0("rm(list=c(",paste0(paste0('"','tterm1', 1:kk,'"'), collapse=', ') ,"))"),
               paste0("rm(list=c(",paste0(paste0('"','tterm3', 1:kk,'"'), collapse=', ') ,"))"),
               paste0("rm(list=c(",paste0(paste0('"','bd', 1:kk,'"'), collapse=', ') ,"))"),
               paste0("rm(list=c(",paste0(paste0('"','meanhe', 1:kk,'"'), collapse=', ') ,"))"),
               "rmlist <- c('rp', 'ad', 'cd', 'dd', 'eps', 'eta', 'ife', 'meange','meanqt', 'merged', 'pedef', 'pife', 'pnorme', 'probs', 'const_b', 'dd_b', 'pedef_b' )",
               paste0("rm(list=ls()[ls()%in%rmlist])"),
               "gc()",
               hess1[(i_h+1):length(hess1)])

  
  ffunch <- function(par){}
  body(ffunch) <- parse(text=hess1)
  list(func=ffuncf,grad=ffuncg, hessian=ffunch, 
       text_func=text_func, text_grad=text_grad, text_hess=text_hess,
       text_bd=text_bd, text_gd=text_gd, text_edd=text_edd)
}



#' @keywords internal
get_FIML_pars <- function(x){
  . <- NULL
  neqc <- length(x)
  ndd <- x %>% gsub("(\\[\\d{1,10}\\])", "XX\\1XX", .) %>% strsplit("XX") %>% unlist %>% unique %>% 
    grep("\\[", ., value = TRUE) %>% gsub("\\[|\\]", "", .) %>% as.numeric() %>% max
  sigmas <- paste0("par[", ndd+1:neqc, "]") 
  tffj <-  matrix(rep(1, neqc), ncol=neqc, nrow=neqc)
  nrho <- length(tffj[lower.tri(tffj)])
  ndd <- ndd+neqc
  rhos <- paste0("par[", ndd+1:nrho, "]") 
  list(sigmas=sigmas, rhos=rhos)
}

#' @keywords internal
replace_FIML <- function(cheqs0, mn1){
  . <- NULL
  fimpar <- get_FIML_pars(cheqs0)
  for(cc in names(mn1)){
    oldpars <- data.frame(new=fimpar$sigmas %>% gsub("\\[|\\]","", .), 
                          old=paste0("sigma_", 1:length(fimpar$sigmas), "_"))
    mn1[[cc]] <- replace_par_wrap(oldpars, mn1[[cc]])
    oldpars <- data.frame(new=fimpar$sigmas, 
                          old=paste0("sigma\\[", 1:length(fimpar$sigmas), "\\]"))
    mn1[[cc]] <- replace_par_wrap(oldpars, mn1[[cc]])
    # replacing rhos
    rhosm <- ff_generate_sigma_expr(length(cheqs0), "rho")
    rhosm <- rhosm[lower.tri(rhosm)] 
    oldpars <- data.frame(old=rhosm %>% formula2string, new= fimpar$rhos %>% gsub("\\[|\\]","", .))
    mn1[[cc]] <- replace_par_wrap(oldpars, mn1[[cc]])
    oldpars <- data.frame(old=rhosm %>% gsub("\\[","\\\\[", .) %>%gsub("\\]","\\\\]", .),
                          new= fimpar$rhos )
    mn1[[cc]] <- replace_par_wrap(oldpars, mn1[[cc]])
  }

  mn1
}

#' @keywords internal
add_check_hessian_method <- function(ffj, methodopt = "BHHH"){
  . <- NULL
  inde <- grep("check_hess == TRUE", ffj, fixed = TRUE)
  if(length(inde)>0){
    iobj <- c('methodo <- methodopt', paste0('methodopt <- "BFGS"'),
              'separatenmm <- FALSE',
              paste0('try(environment(joint_hess) <- environment(), silent = TRUE)'))
    ffj <- insert_obj(inde, obj = ffj, insert = iobj)
    iobj <- c('methodopt <- methodo', 'try(environment(joint_hess) <- environment(), silent = TRUE)')
    inde <- grep("hess <- try(qr.solve", ffj, fixed = TRUE)
    ffj <- insert_obj(inde, obj = ffj, insert = iobj)
  }
  ffj
}


#' \code{LL_joint_no_corr} Function for log-likelihood without correlation between
#'  continuous and discrete equations
#' @param ffor Discrete equations.
#' @param cheqs0 Continuous equations.
#' @param datan Character string. Name of data-set.
#' @param fixed_term if \code{TRUE}, includes fixed term in log-likelihood.
#' @param eq_type can be "joint"(both blocks estimated but without correlation), "cont", "disc".
#' @param transform if \code{TRUE}, quantile transformation is applied to discrete equations.
#' @param hessian Name of hessian function.
#' @param bayesian_random If \code{TRUE}, than par[1] is changed to par[,1] to be used for 
#' optimization of random parameters in Bayesian estimation.
#' @param eqsys "sur" or "sem".
#' @param MNtypef "logit" or "dogit".
#' @param add_pars additional parameters to add
#' @return ML estimation function
#' @keywords internal
LL_joint_no_corr <- function(ffor, cheqs0, datan = "data", fixed_term = TRUE, 
                             eq_type = c("joint", "cont", "disc"), transform = FALSE, 
                             hessian = NULL, bayesian_random = FALSE,
                            eqsys = "sur", MNtypef = "logit", add_pars = NULL)
{
  . <- NULL
  #cont. part
  methodopt <- "NA"
  if(eq_type[1]!="disc"){
    #ffmn_p <- if(eqsys=="sur") maxle_p else mlsem
    ffmn <- if(eqsys=="sur") maxle else mlsem
    mn <- maxle_p(cheqs0 = cheqs0, fixed_term = fixed_term) # partitioned LL
    mn1 <- ffmn(cheqs0 = cheqs0, fixed_term = fixed_term) # separatenmm LL for each equation
    parlc <- get_npar(cheqs0, values=TRUE)
    parlc <- c(parlc, add_pars)
    
    npar <- length(parlc)
  }else{
    mn <-list(formula=NULL)
    parlc <- NULL
    npar <- 0
    mn1 <- list(formula=NULL)
    cheqs0 <- NULL
  }

  #disc. part
  if(eq_type[1]!="cont"){
    fforn <- gsubfn('(par\\[)(\\d{1,10})(\\])', function(x, y, z)paste0(x, as.numeric(y) + npar, z),
                    ffor)
    #mnf1 <- MNlogitf(fforn, transform = transform, separatenmm = TRUE,weight_disc = TRUE)
    #mnf1 <- MNdogitf(fforn, transform = transform, separatenmm = TRUE, weight_disc = TRUE)
    if(MNtypef=="logit"){
      mnf1 <- MNlogitf(fforn, transform = transform, separatenmm = TRUE, weight_disc = TRUE)
      ppardogit <- NULL
    }else{
      ndd <- fforn %>% gsub("(\\[\\d{1,10}\\])", "XX\\1XX", .) %>% strsplit("XX") %>% unlist %>% unique %>% 
        grep("\\[", ., value = TRUE) %>% gsub("\\[|\\]", "", .) %>% as.numeric() %>% max
      ppardogit <- c(paste0('par[',(ndd+1):(ndd+length(fforn)),']'))
      print(ppardogit)
      mnf1 <- MNdogitf(fforn, transform = transform, separatenmm = TRUE,
                       weight_disc = TRUE, ppardogit)
    }
    
    parld <- get_npar(fforn, values=TRUE)
    fixed <- npar+1
    summet <- mnf1$sume
  }else{
    mnf1 <- list(formula=NULL)
    parld <- NULL
    fixed <- 0
    ffor <- NULL
    summet <- NULL
    ppardogit <- NULL
  }

  if(length(mn$formula)>0){
     new_f1 <- c(mn$formula%>%paste0('(', ., ')*(1/npaths_cont)'), mnf1$formula)
     new_f <- c(mn$formulap%>%paste0('(', ., ')*(1/npaths_cont)'), mnf1$formula)%>%
       paste0('(', ., ')*(1/npaths)')
  }else{
    new_f1 <- c( mnf1$formula)
    new_f <- c(mnf1$formula)%>%paste0('(', ., ')*(1/npaths)')
  }
  new_f1%<>%paste0(., collapse='+')%>%paste0('(', ., ')*(1/npaths)')

  parll <- c(parlc, parld, ppardogit[ppardogit!=0])
  indexp <- gsub("par\\[|\\]", '', parll)
  parll <- parll[indexp%>%as.numeric%>%order]
  parln <- parll%>%gsub('\\[|\\]', '', .)

  mn1 <- list(expr= formula2string(new_f1), formula=new_f1)

  #combine both
  joint_func <- f_create(new_f, data = datan, cheqs0=cheqs0, transform = FALSE, fixed = fixed,
                         separatenmm = TRUE, hessian = hessian, sume=summet)
  ffj <- deparse(joint_func)
  ffj <- add_check_hessian_method(ffj)
  # add rho and sigmas to the estimation
  
  inde <- grep("if \\(separatenmm == FALSE\\)", ffj)-1
  if(length(inde)==0){
    inde <- min(grep("return\\(res\\)", ffj))-1
    iobj <- c('atrb <- res','attributes(res) <- c(attributes(res), list(LLv=atrb))')
    ffj <- insert_obj(inde, ffj, iobj)
  }else{
    iobj <- 'atrb <- res'
    ffj <- insert_obj(inde, ffj, iobj)

    inde <- min(grep("return\\(res\\)", ffj))-1
    iobj <- 'attributes(res) <- c(attributes(res), list(LLv=atrb))'
    ffj <- insert_obj(inde, ffj, iobj)
  }
  ffj <- ffj[-1]
  if(bayesian_random==TRUE){
    actbr <- .%>%gsub("par\\[(\\d{1,10})", "par[,\\1",.)%>%
      gsub("par\\[\\((\\d{1,10})", "par[,(\\1",.)%>%
      gsub("par <- c(", "par <- cbind(", ., fixed=TRUE)%>%
      gsub("length(par)", "ncol(par)", ., fixed=TRUE)
    ffj%<>%actbr
  }
  joint_func <- function(par){}
  body(joint_func) <- parse(text=ffj)
  
  gf <- grad_hess_eval(mn = mn1, parnl = parln, fixed=fixed, data=datan, cheqs0=cheqs0)
  hf <- grad_hess_eval(mn = mn1, parnl = parln, hessian = TRUE, fixed=fixed, data=datan,
                       cheqs0=cheqs0)
  if(bayesian_random==TRUE){
    gf1 <- gf%>%deparse%>%actbr
    indx <- grep("methodopt", gf1)
    indx <- indx[length(indx)]
    #gf1 <- c(gf1[1:(indx-1)], '    methodopt <-"XXX"', gf1[indx:length(gf1)]) # NEEDED?
    gf1 <- gf1[-c(1)]%>%parse(text=.)
    gf <- function(par){}
    body(gf) <- gf1

    hf1 <- hf%>%deparse%>%actbr
    indx <- grep("methodopt", hf1)
    indx <- indx[length(indx)]
    #hf1 <- c(hf1[1:(indx-1)], '    methodopt <-"XXX"', hf1[indx:length(hf1)]) # NEEDED?
    hf1 <- hf1[-c(1)]%>%parse(text=.)
    hf <- function(par){}
    body(hf) <- hf1
  }
  list(func=joint_func,grad=gf, hessian=hf)
}
