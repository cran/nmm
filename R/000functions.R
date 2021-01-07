#' \code{check_par_av} Checks if all elements of par are available in object.
#' @keywords internal
#' @importFrom dplyr mutate select starts_with filter
#' @importFrom plyr ddply 
#' @param par String vector, what to look for.
#' @param object String vector, where to look.
check_par_av <- function(par, object){
  . <- NULL
  av_par <- sapply(par, function(x)object %>% grepl(x, .) %>% any)
  if(av_par %>% all %>% isFALSE){
    par[!av_par] %>% paste0(., collapse = ", ") %>%
      paste0("Parameter(s) ", ., " not found in the supplied equations!") %>% print
    stop()
  }
}

#' \code{get_par} replaces names of parameters with par[i].
#' @import magrittr
#' @param par Names of parameters, vector of character strings.
#' @param object List consisting formulas. 
#' @return A list object consisting of equations for errors with par[i] (\code{cheqs0}), 
#' the original list with formulas (\code{eqlab}), vector of parameters (\code{parld}, same as par),
#' vector of parameters in form of par[i] (\code{parn}, exogenous variables (\code{exog}), 
#' endogenous variables (\code{endog}, in case of discrete ""), 
#' number of parameters in each equation (\code{neq_par}).
#' @examples
#' # System of Non-linear Regressions
#' eq_c <- c("hg ~ exp( h0 + h1*log(tht) + h2*tht^2 + h3*elev)",
#' "dg ~ exp( d0 + d1*log(dbh) + d2*hg + d3*cr)")
#' par_c <- c(paste0("h", 0:3),paste0("d", 0:3))
#' para_cont <- get_par(par_c, eq_c)
#' # Indirect utility functions for discrete choice:
#' eq_d <- c("a1 + p1 * price_1 + p2 * catch_2", "a2 + p1 * price_2 + p2 * catch_2",
#' "a3 + p1 * price_3 + p2 * catch_3", "a4 + p1 * price_4 + p2 * catch_4")
#' par_d <- c(paste0("a", 1:4), paste0("p", 1:2))
#' disc_par <- get_par(par_d, eq_d)
#' @export
get_par <- function(par, object)
{
  . <- NULL
  par1 <- par
  check_par_av(par1, object)
  type <- "cont"
  # check if cont or discrete, if discrete no ~ available, thus add it
  if(object %>% grepl("~", .) %>% all %>% isFALSE){
    type <- "disc"
    object %<>% paste0("eq", 1:length(.), " ~ ", .)
    par %<>% sort
  }
  eqlab <- object %>% lapply(. %>% as.formula)
  npar <- length(par)
  # in eqlab replace elements of parl with elements of para
  para <- 1:npar %>% paste0('par[', ., ']')
  names(para) <- par
  cheqs <- eqlab %>% as.character
  repdat <- data.frame(old = par, new = para, stringsAsFactors = FALSE)
  #replace parameters with new strings par[x]
  cheqs <- replace_par_wrap(repdat, cheqs)
  tmp <- cheqs %>% as.character %>% strsplit(., '~')
  cheqs0 <- paste0(tmp %>% lapply(., "[[", 1), "-", tmp %>% lapply(., "[[", 2)) %>% as.list
  # how many parameters in each equation
  neq_par <- cheqs %>%lapply(., function(x){
    x %>% gsub("(par\\[\\d{1,10}\\])", "XX\\1XX", .) %>% strsplit(., "XX") %>%
      unlist %>% grep("par\\[",., value = TRUE ) %>% unique %>% length
  })

  # exogenous variables
  p4_ <- "XXZZXX"
  exog <- object %>% gsub("_", p4_ , .) %>% gsub("(\\w)\\(", paste0("\\1", "XXXZZZ "), .) %>%
    gsub('[^[:alnum:]]', ' ', .) %>%  strsplit(., ' ') %>% unlist %>% unique %>%
    grep('^[^[:digit:]]', ., value = TRUE) %>% gsub(p4_, "_", .)

  endog <- object %>% gsub("\\~.*", "", .) %>% trimws
  exog <- exog[(!exog %in% c(par, "", endog, "")) & !grepl("XXXZZZ", exog)]
  if(type == "disc"){
    endog <-  ""
    cheqs0 <- cheqs0 %>% gsub("^eq\\d{1,10} - ", "", .)
    eqlab <- object %>% gsub("^eq\\d{1,10} ~ ", "", .)
  }

  list(cheqs0 = cheqs0, eqlab = eqlab, parld = repdat$old, parn = repdat$new, exog = exog,
       endog = endog, neq_par = neq_par)
}


#' \code{get_npar} Get number of parameters or vector of parameters in supplied equations.
#' Extracts the number of parameters used in equations. Parameters are given as par[1], ..., par[n].
#' @param x List of strings. 
#' @param values if \code{TRUE} returns the character vector with parameters (par[i]).
#' @return Number of parameters or vector with parameters (in form par[i]).
#' @examples
#' eq_d <- c("ASC1 * 1 + B11_dur * dur_1" , "ASC2 * 1 + B12_dur * dur_2",
#' "ASC3 * 1 + B13_dur * dur_3 + B20_cost * cost_3 + B53_parkman * PbAvl_3",
#' "ASC4 * 1 + B14_dur * dur_4 + B20_cost * cost_4 + B34_serv * servIdx_4 + B44_stop * stopUs1R1_4")
#' parl <- c("ASC1", "B11_dur", "ASC2", "B12_dur", "ASC3", "B13_dur", "B20_cost", "B53_parkman",
#' "ASC4", "B14_dur", "B20_cost", "B34_serv", "B44_stop") %>% unique
#' disc_par <- get_par(parl, eq_d)
#' get_npar(disc_par$cheqs0)
#' get_npar(disc_par$cheqs0, values=TRUE)
#'
#' @export
get_npar <- function(x, values=FALSE)
{
  . <- NULL
  xx <- regmatches(x %>% unlist, gregexpr("par\\[.*?\\]", x %>% unlist)) %>% unlist %>% unique
  if(values == TRUE) xx else xx %>% length
}

#' \code{formula2string} removes square brackets from the supplied expressions.
#' Convert par[2] -> par2, sigma[2] -> sigma_2_, sigma[2,3] -> sigma_2x2
#' @param x String with square brackets.
#' @return String without square brackets.
#' @examples
#' xx <- "par[1]*3+par[2]*par+rho1[1,2]+sigma1[2]"
#' formula2string(xx)
#' @export
formula2string <- function(x)
{
  . <- NULL
  x %<>% gsub('(par)\\[(\\d{1,10})\\]','\\1\\2', .)
  x %<>% gsub('([[:alpha:]]|\\d)\\[(\\d{1,10})\\]','\\1_\\2_', .)
  x %<>% gsub('([[:alpha:]]|\\d)\\[(\\d{1,10}),([[:digit:]]{1,10})\\]','\\1_\\2x\\3_', .)
  x
}

#' \code{string2formula} add square brackets to expressions.
#' Reverse of formula2string. Convert par[2] <- par2, sigma[2] <- sigma_2_, sigma[2,3] <- sigma_2x2
#' @param x String without square brackets.
#' @return String with square brackets.
#' @examples
#' xx <- "par[1]*3+par[2]*par+rho1[1,2]+sigma1[2]"
#' xm <- formula2string(xx)
#' string2formula(xm)
#' @export
string2formula <- function(x)
{
  . <- NULL
  x %<>% gsub('(par)(\\d{1,10})','\\1[\\2]', .)
  x %<>% gsub('(sigma)(\\d{1,10}|)_(\\d{1,10})x(\\d{1,10})_', '\\1\\2[\\3,\\4]', .)
  x %<>% gsub('(sigma)(\\d{1,10}|)_(\\d{1,10})_', '\\1\\2[\\3]', .)
  x %<>% gsub('(rho)(\\d{1,10}|)_(\\d{1,10})x(\\d{1,10})_', '\\1\\2[\\3,\\4]', .)
  x
}

#' \code{wxMaxima} does symbolic computation in 'Maxima'
#' Requires installation of Maxima software.
#' @param obj Lines without printed output, that will be evaluated in wxMaxima.
#' @param out Lines with printed output, is end results of the function.
#' @param tex if \code{TRUE} TeX expression of the printed output.
#' @return List of character strings.
#' @examples
#' #components(determinant, Sigma inverse, argument for exponent) for joint normal distribution
#' eq_c <- c("Tw ~ ((((PH) + (tw)) * (ta - Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) -(1 + (PH))) +
#' sqrt((((PH) + (tw)) * (ta - Tc + 2) + (1 +(tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 - 4 * (1 + (PH) +
#' (tw)) *(-(PH) * (ta - Tc + 2) + (1 - (tw) * (ta - Tc + 2)) * (2/w -Ec/w))))/(2 * (1 + (PH) +
#' (tw)))",
#' "Tf1 ~ (th1) * (ta - (((((PH) + (tw)) * (ta - Tc + 2) + (1 + (tw)) *(Ec/w - 2/w) - (1 + (PH))) +
#'  sqrt((((PH) + (tw)) * (ta -Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 - 4 *(1 + (PH) +
#'   (tw)) * (-(PH) * (ta - Tc + 2) + (1 - (tw) *(ta - Tc + 2)) * (2/w - Ec/w))))/(2 * (1 + (PH) +
#'   (tw)))) -Tc + 2) - 1",
#' "Ef1 ~ (ph1)/(PH) * (w * (((((PH) + (tw)) * (ta - Tc + 2) + (1 +(tw)) * (Ec/w - 2/w) -
#' (1 + (PH))) + sqrt((((PH) + (tw)) *(ta - Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 -4 *
#'  (1 + (PH) + (tw)) * (-(PH) * (ta - Tc + 2) + (1 - (tw) *(ta - Tc + 2)) * (2/w - Ec/w))))/(2 *
#'  (1 + (PH) + (tw)))) -Ec + 2) - 1")
#' parl <- c("tw","PH","th1","ph1")
#' para_cont <- get_par(parl, eq_c)
#' cheqs0 <- para_cont$cheqs0
#' npar <- get_npar(cheqs0)
#' neq <- length(cheqs0)
#' sdv <- rep(NA, neq)
#' mv <- rep(NA, neq)
#' sigma <- expand.grid(1:neq, 1:neq)%>%apply(., 1, function(x)paste0(x, collapse=','))%>%
#' paste0('sigma', '[', . ,']')%>%sort
#' sigma <- matrix(sigma, neq, neq, byrow = TRUE)
#' sigma[lower.tri(sigma)] <- sort(sigma[upper.tri(sigma)])
#' #create Y
#' y <- paste0('eps[', 1:neq,']')
#' ypy <- paste0('[', y, ']', collapse = ',')
#' ypy <- paste0('Y : matrix(', ypy, ')')
#' #sigma
#' spy <- paste0(sapply(1:neq, function(i) paste0('[', paste0(sigma[i,], collapse = ', '), ']')),
#' collapse = ',')
#' spy <- paste0('Sigma : matrix(', spy, ')')
#' detS <- 'dets : ratsimp(determinant(Sigma))'
#' invS <- 'invs : ratsimp(invert(Sigma))'
#' expt <- 'expt : ratsimp(transpose(Y) . invs . Y)'
#' obj <- c(ypy, spy, detS, invS, expt)
#' #grind function in Maxima returns an object that can be mathematically evaluated
#' out <- c('print(new)', 'grind(dets)',  'print(new)', 'grind(invs)', 'print(new)', 'grind(expt)')
#' \dontrun{
#' rez <- wxMaxima(obj, out)
#' }
#' @export
wxMaxima <- function(obj, out, tex=FALSE)
{
  . <- NULL
  sysn <- Sys.info()["sysname"] #get system
  #execute lines in wxMaxima
  ngrind <- out %>% grep("grind", ., value=TRUE)

  # for long expressions different treatment, then eq>4 does not function otherwise
  # will not work, as only 8000 characters allowed to cmd..
  curr_dir <- tempdir() %>% gsub("\\\\", "/", .)
  tmpdirxxx <- "tmp_XXXXXX"
  dir.create(paste0(curr_dir, "/", tmpdirxxx), showWarnings = FALSE)

  out[out%in%ngrind] <- (1:length(ngrind)) %>%sapply(.,function(j) gsub("grind\\((.*.)\\)",
                              paste0('stringout("',paste0(curr_dir, "/", tmpdirxxx, "/temp", j,
                                                          ".txt" ), '", \\1)'), ngrind[j]))

  ntex <- out %>% grep("tex \\(", ., value=TRUE)
  out[out%in%ntex] <-
    (1:length(ntex)) %>%sapply(.,function(j) gsub("tex \\((.*.)\\)",
                                          paste0('tex (\\1, "',paste0(curr_dir, "/", tmpdirxxx,
                                                                      "/tex",j, ".txt" ), '")'),
                                           ntex[j]))

  
  if(grepl("windows", sysn, ignore.case = TRUE)){
    sink(paste0(curr_dir, "/", tmpdirxxx, "/main.txt" ))
    c(obj, out)%>%paste0(., " ; ") %>% cat(., "\n")
    sink()
    com <- paste0('batch("',paste0(curr_dir, "/", tmpdirxxx, "/main.txt" ),'");')
    comm <- paste0('echo ', com , ' | maxima')
    a <- shell(comm, intern = TRUE)

  }else{
    # maxima -r 'print(new);'
    #https://linux.die.net/man/1/maxima
    com <-  c(obj, out, 'quit()')%>%paste0(., " ; ")%>%paste0( collapse = "")
    comm <- paste0("maxima -r '", com, "'")
    a <- system(comm, intern = TRUE)
  }
  ## read in the files
  ffd <- list.files(paste0(curr_dir, "/", tmpdirxxx), pattern = "temp")
  aa <- lapply(ffd, function(x)readLines(paste0(curr_dir, "/", tmpdirxxx, "/", x)))
  aa %<>% unlist
  aa <- aa[aa!=""]
  aa <- aa[sapply(aa, function(x)!any(grepl("^\\$\\$|\\$\\$$",x)))] #remove latex blocks, begin w $
  aa <- lapply(aa, function(x)x[!grepl("done|new|false", x)]) #remove some unnec. output
  aa <- lapply(aa, function(x) x%>%paste0(., collapse = "") %>%
                 gsub("(  )|\\$$|\\;","",.))#remove $ " "
  names(aa) <- 1:length(aa)
  # tex files
  ffd <- list.files(paste0(curr_dir, "/", tmpdirxxx), pattern = "tex")
  aa2 <- lapply(ffd, function(x)paste0(curr_dir, "/", tmpdirxxx, "/", x) %>% readLines %>%
                  paste0(., collapse = ""))
  aa2 %<>% unlist
  aa2 <- aa2[aa2!=""]
  aa2 <- lapply(aa2, function(x)x[!grepl("false|done|new", x)])
  aa2 <- lapply(aa2, function(x)paste0(x, collapse = ""))
  unlink(paste0(curr_dir, "/", tmpdirxxx), recursive = TRUE)

  if (length(aa2)>0) names(aa2) <- 1:length(aa2)
  aa2 <- if(tex==FALSE) list(mean="", cov="") else aa2
  list(expr=aa, latex=aa2)
}

#' \code{cond_expr} returns moments of conditional multivariate normal distribution X|Y (last 
#' variable is dependent). Only expression for X|Y. Requires installation of Maxima software. 
#' @param neq Number of equations/variables.
#' @param sdv Vector of standard deviation of normally distributed variables, e.g. c(NA, NA, NA, 1)
#' NA - unknown, any number - know.
#' @param mv Vector of means of normally distributed variables, e.g. rep(0, 4).
#' @param nconteq Number of continuous equations. 
#' @param tex i if \code{TRUE} TeX expressions from wxMaxima are returned.
#' @return List of strings. First element is an expression of conditional mean and covariance. 
#' The second element is a TeX formula.
#' @examples
#' # this means that E[y3|y1,y2] and V[y3|y1,y2] will be returned
#' # all continuous w/ unknown means
#' \dontrun{
#' # To run this, one needs to install Maxima software
#' res <- cond_expr(neq=3)
#' # 3 continuous w/ unknown means and the last one with mean 0 and sd 1, d|c1c2c3
#' res <- cond_expr(neq=4, sdv=c(NA, NA, NA, 1), mv=c(NA, NA, NA, 0))
#' # 2 continuous w/ unknown means and 2 discrete with mean 0 and sd 1, d1|c1c2c3d2
#' res <- cond_expr(neq=4, sdv=c(NA, NA,  1, 1), mv=c(NA, NA,  0, 0), nconteq=2)
#' }
#' @export
cond_expr <- function(neq, sdv, mv, nconteq = neq - 1, tex = FALSE)
{
  . <- NULL
  mf <- match.call()
  m <- match(c("neq", "sdv", "mv"), names(mf), 0L)
  if(m[2]==0){
    sdv <- rep(NA, neq)
  }
  if(m[3]==0){
    mv <- rep(NA, neq)
  }
  # create symmetric sigma matrix
  sigma <- paste0(rep("sigma", neq^2),"[",sort(rep(1:neq, neq)) ,"]*",
                  rep("sigma", neq^2),"[",rep(1:neq, neq) ,"]*",
                  rep("rho[",neq^2),sort(rep(1:neq,neq )),",",rep(1:neq, neq),"]" )
  if(nconteq==neq-2){
    sigma <- gsub(paste0("rho\\[",neq,"\\,",neq-1,"\\]"), "0", sigma)
    sigma <- gsub(paste0("rho\\[",neq-1,"\\,",neq,"\\]"), "0", sigma)
  }
  sigma <- matrix(sigma, neq, neq, byrow = TRUE)

  diag(sigma)%<>%gsub("rho\\[\\d{1,10},\\d{1,10}\\]","",.)
  diag(sigma)%<>%gsub("[*].*","",.)%>%paste0(.,"**2")
  sigma[lower.tri(sigma)] <- sort(sigma[upper.tri(sigma)])

  #dependent(conditional) variable X is the last variable X|Y
  depind <- neq
  #independent variable Y
  indind <- if(neq==1) 1 else 1:(neq-1)
  #if any sd is given
  indsd <- which(!is.na(sdv))
  for(i in indsd){
    sigma <- gsub(paste0("sigma[[]",i,"[]]"),sdv[i], sigma)
  }
  #mean vector
  mvec <- 1:neq%>% paste0("mu[", . ,"]")%>%paste0("[", ., "]")
  indm <- which(!is.na(mv))
  #if any mean is given
  for(i in indm){
    mvec <- gsub(paste0("mu[[]",i,"[]]"),mv[i], mvec)
  }

  #sigmas for wxmaxima
  spy <- apply(sigma, 1, function(x)paste0("[", paste0("",x,"", collapse = ", "), "]"))%>%
    paste0(., collapse = ",")
  spy %<>% paste0("Sigma : matrix(",., ")")

  #X|Y
  sxy <- sigma[depind, indind]%>%paste0("[",.,"]")%>%paste0(., collapse = ",")%>%
    paste0("sigmaxy : matrix(",.,")")
  #Y|Y
  syy <- sapply(indind, function(i) paste0("[",
                                           paste0("",sigma[i,-depind],"", collapse = ", "), "]"))%>%
    paste0(., collapse = ",")
  syy %<>% paste0("sigmayy : matrix(", ., ")")
  #X|X
  sxx <- sigma[depind, depind]%>%paste0("sigmaxx : matrix([",.,"])")
  # if(any(!is.na(sdv))){
  #   sxx <- sdv[depind]%>%paste0("sigmaxx : matrix([",.,"])")
  # }
  

  #means
  mx <- mvec[depind]%>%paste0(., collapse=",")%>%paste0("meanx : matrix(", .,")")
  my <- mvec[indind]%>%paste0(., collapse=",")%>%paste0("meany : matrix(", .,")")

  # create Y, if discrete equations are also conditioned, neq=4, 2 cont and 2 discrete
  if(nconteq < (neq-1)){
    deq <- 1:(neq-nconteq-1)
    y <- c(paste0("eps", 1:nconteq, ""), paste0("disc", deq, ""))
  }else{
    y <- paste0("eps", indind, "")
  }
  ypy <- y%>%paste0("[", . , "]", collapse = ",")%>%paste0("y : matrix(", ., ")")
  # for matrix multimplication maxima uses ".", for scalar "*"
  signp <- if(neq <= 2) "*" else "."

  #conditional mean and cov see under 19 eq Diaz 2008
  mc <- signp%>%paste0("mean : ratsimp(meanx + transpose(sigmaxy)",  .,"invert(sigmayy)", . ,
                       "(y-meany))")
  mc <-  if(neq==1) paste0("mean : meanx") else mc
  cc <- signp%>%paste0("cov : ratsimp(sigmaxx-transpose(sigmaxy)",  .,"invert(sigmayy)", . ,
                       "(sigmaxy))")
  cc <-  if(neq==1) paste0("cov : sigmaxx") else cc
  mct <- "mean[1,1]"
  cct <- "cov[1,1]"
  mce <- "grind(mean[1,1])"
  cce <- "grind(cov[1,1])"

  obj <- c( sxy, syy, sxx, ypy, my, mx)
  out <- c("print(new)", mc, cc, "print(new)",mct, cct,  "print(new)", mce, cce)
  out <- c(out, paste0("tex (''mean[1,1])"),
           paste0("tex (''cov[1,1])"))
  rez <- wxMaxima(obj, out, tex = tex)
  names(rez[[1]]) <-c("mean", "cov")
  names(rez[[2]]) <-c("mean", "cov")
  rez
}

#' \code{cond_mean_cov_expr} returns conditional mean and variance X|Y of conditional multivariate 
#' normal distribution . X is discrete variable, Y are continuous variables. 
#' No correlation between discrete eq.
#' @param neqt Number of continuous variables/equations
#' @param kk Number of discrete variables/equations
#' @return list:
#' \item{meannewel}{conditional mean expressions, w/ "<-"}
#' \item{nomin}{general form of nominator terms, w/o "<-"}
#' \item{nominl}{for each i in 1:kk nominator terms, w/ "<-"}
#' \item{denoml}{denominator in conditional mean, w/ "<-"}
#' \item{denom}{denominator in conditional mean, w/o "<-"}
#' \item{cel}{conditional covariance expressions, w/ "<-"}
#' @keywords internal
cond_mean_cov_expr <- function(neqt, kk)
{
  . <- NULL
  if(neqt>4){
    sdv <- c(rep(NA, neqt), 1)
    mv <- rep(0, neqt+1)
    arg1 <- neqt+1
    arg4 <- arg1-1
    out <- cond_expr(arg1, sdv, mv, arg4)
  }else{
    out <- dat4cond_mean_cov_expr[[neqt]]
  }
  meane <- out$expr$mean
  cove <- out$expr$cov

  #here eps has a meaning from cheqs0, etas in the paper
  mact1 <- .%>%gsub('eps(\\d{1,10})', 'eta[\\1]', .)%>%gsub('disc(\\d{1,10})', 'disc[\\1]', .)
  meane %<>% mact1
  cove %<>% mact1
  meanee <- meane%>%strsplit(.,"/")%>%unlist
  denom <- meanee[[2]]
  nomin <- meanee[[1]]%>%gsub('^\\(|\\)( |)$', '',.)
  #with minus
  minn1 <- (grepl('^-\\(',meanee[[1]])*grepl('( |)$',meanee[[1]]))%>%as.logical
  minn <- meanee[[1]]%>%grepl('^-\\(.*( |)$',.)
  if(minn1!=minn) browser()
  nomin[minn] <- nomin[minn]%>%gsub('^-\\(|\\)( |)$','', .)
  nomin <- strsplit(nomin, 'eta\\[\\d\\]|disc\\[\\d\\]')[[1]]%>%gsub('^( |)\\+|\\*( |)$', '', .)%>%
    paste0('(', ., ')')
  if(minn==TRUE){
    nomin%<>%paste0("(-(", ., "))")
  }
  etase <- meanee[[1]]%>%gregexpr("eta\\[\\d\\]|disc\\[\\d\\]", .)%>%regmatches(meanee[[1]], .)%>%
    unlist

  nomin <- nomin[order(etase)]
  meannewe <- 1:neqt%>%paste0('(nomin[', . , ']*eta[[', .,']])', collapse='+')

  meannewe%<>%paste0('(',., ')','/denom')%>%paste0(., collapse=" + ")
  cel <- c()
  nominl <- c()
  meannewel <- c()
  for(k in 1:kk){
    cel <- c(cel,cove%>%gsub("(rho)|(sigma)", paste0("\\1\\2",k), .)%>%paste0("covqt",k, " <- ",.))
    meannewel <- c(meannewel, meannewe%>%gsub('(nomin)(\\[\\d\\])',paste0('\\1',k,'\\2'), .)%>%
                     paste0("meanqt",k, " <- ",.))
    nominl <- c(nominl, nomin%>%gsub("(rho)|(sigma)", paste0("\\1\\2",k), .)%>%
                  paste0("nomin", k, "[",1:length(nomin),"] <- ",.))
  }
  cel%<>%gsub("\\[","[[", .)%>%gsub("\\]","]]", .)
  denoml <- paste0('denom <- ',denom)
  return(list(meannewel=meannewel, nomin=nomin, nominl=nominl,denoml=denoml, denom=denom, cel=cel))
}



#' \code{extract_attr_deriv} converts attributes(hessian/gradient) of \code{deriv()} into a 
#' matrix of character strings.
#' @param ex Expression of derivative. Results of \code{deriv()}.
#' @param attribute "grad" for gradient or "hessian" for the Hessian matrix.
#' @return Returns a matrix of character strings.
#' @examples
#' eq <- parse(text="2*(log(sin(x)/log(x)))+x^4*log(x)+cos(y+x)")
#' tt <- deriv(eq, c("x", "y"), hessian=TRUE)
#' g <- tt%>%extract_attr_deriv(., attribute = "grad")
#' h <- tt%>%extract_attr_deriv(., attribute = "hessian")
#' @export
extract_attr_deriv <- function(ex, attribute)
{
  . <- NULL
  exx <- ex %>% deparse %>% gsub("structure\\(", "", .)%>%
    gsub("expression\\(\\{", "", .)%>%
    gsub(paste0("  [.]", attribute, "$"), "", .)

  exx <- exx[trimws(exx)!=""]
  #lines ending with " ", means that the expression continuous into the next line
  # we want to collapse
  ind <- grep(" $", exx)
  i <- ind[1]
  exx <- collapse_space(exx, i, ind)
  exx <- gsub("(, srcfile =).*", "", exx)
  exx <- exx[!grepl("\\}\\)$", exx)]

  #expressions
  expre <- grep("[.]expr.*.<-", exx, value=TRUE)
  if(length(expre)>0){
    expre <- strsplit(expre, "<-")
    expre <- data.frame(expr=sapply(expre, "[",1), val=sapply(expre, "[",2),
                        stringsAsFactors = FALSE)
    expre <- matrix(apply(expre, 2, trimws), ncol=2)
    inde <- grep(".expr", expre[,2])

    #replace expressions in expressions
    expre <- replace_expression(expre, inde)
    exx <- grep(attribute, exx, value=TRUE)
    exx <- exx[!grepl("attr", exx)]
    inde <- grep(".expr", exx)

    #replace expressions in expressions
    expre1 <- data.frame(expre, stringsAsFactors = FALSE)
    names(expre1)<- c("old", "new")
    expre1[, "new"] <- paste0("(", expre1[, "new"], ")")
    exx1 <- replace_par_wrap(expre1, exx)
    exx0 <- exx
    #remove after testing
    while(length(inde)>0){
      for(ee in expre[,1]){
        if(length(expre[expre[,1]==ee,2])>0){
          exx[inde] <- sapply(exx[inde], function(x)
            gsub(paste0(ee,"( |$|\\)|\\^|/)"), paste0("(", expre[expre[,1]==ee,2],")\\1"), x) )
        }
      }
      inde <- grep(".expr", exx)
    }
    if(any(exx1!=exx)) browser()
  }

  if(attribute=="grad"){
    exx %<>% grep(attribute, ., value=TRUE)%>%
      grep("attr", ., value=TRUE, invert = TRUE)
  }

  indx <- !grepl("array", exx)&!grepl(paste0(attribute, "*.*", attribute), exx)
  exx[indx] <- gsub("(<- )|$","\\1 \"", exx[indx] )
  da <- attribute%>%paste0(., "*.*", .)%>%grepl(.,exx)
  if(any(da)){
    exx[da] <- gsub("(.*.)(<-)(.*.)(<-)(.*)","\\1\\2\\3\\4 \"\\5\"", exx[da])
  }

  exx %<>% gsub("length[(][.]value[)]", "1", .)%>%gsub('[.]value "$', "", .)
  eval(parse(text=c("{", exx, paste0("    .",attribute),"}")))
}

#' Collapses lines ending with " "
#' @keywords internal
#' @param exx expression
#' @param i first element of ind
#' @param ind vector of lines ending with " "
collapse_space <- function(exx, i, ind)
{
  . <- NULL
  if(length(ind)>0){
    r <- paste0(exx[i], gsub("        ","" ,exx[i+1]), collapse = "")
    si <- if(i==1) 0 else 1
    if((i+2)>length(exx)){
      exx <- c(exx[si:(i-1)], r)
    }else{
      exx <- c(exx[si:(i-1)], r, exx[c((i+2):length(exx))])
    }
    ind <- grep(" $", exx)
    i <- ind[1]
  }
  if(length(ind)>0){
    return(collapse_space(exx, i, ind))
  }else{
    return(exx)
  }
}

#' Replaces expressions ".exprXX" with strings
#' @keywords internal
#' @param exxe expression
#' @param inde vector of lines having expression
replace_expression <- function(expre, inde)
{
  . <- NULL
  if(length(inde)>0){
    for(j in inde){
      tt <- strsplit(expre[j,2], "[[:punct:]]")[[1]]
      tt[grepl("expr", tt)] <- paste0(".", grep("expr", tt, value=TRUE) )
      tt <- grep(".expr", tt, value=TRUE)
      for(t1 in tt){
        t1 <- trimws(t1)
        expre[j,2] <- gsub(paste0(t1,"( |$|\\)|\\^|/)"),
                           paste0("(", expre[expre[,1]==t1,2],")\\1"), expre[j,2])
      }
    }
    inde <- grep(".expr", expre[,2])
  }
  if(length(inde)>0){
    return(replace_expression(expre, inde))
  }else{
    return(expre)
  }
}

#' \code{combine_attr_deriv} combines attributes of two \code{deriv()} results.
#' @keywords internal
#' @param t1 Derivative of first expression. Results of \code{deriv()}.
#' @param t2 Derivative of second expression. Results of \code{deriv()}.
#' @param attribute "grad" for gradient or "hessian" for the Hessian matrix.
#' @return Combined matrices of attributes.
combine_attr_deriv <- function(t1, t2, attribute)
{
  . <- NULL
  t1a <- extract_attr_deriv(t1, attribute)
  t2a <- extract_attr_deriv(t2, attribute)
  varl <- unique(c(colnames(t1a), colnames(t2a)))

  if(attribute=="grad"){
    res <- array(0, c(1,length(varl)), list(NULL, varl ))

    for(i in varl){
      t1ai <- ifelse(grepl(i, colnames(t1a))%>%any, t1a[,i], "0")
      t2ai <- ifelse(grepl(i, colnames(t2a))%>%any, t2a[,i], "0")
      tai<- c(t1ai, t2ai)
      res[,i]<- paste0(tai, collapse=" + ")
    }
  }
  if(attribute=="hessian"){
    res <-   array(0, c(1,length(varl),length(varl)), list(NULL, varl ,varl))
    for(i in varl){
      for(j in varl){
        t1ai <- ifelse(try(t1a[,i, j], silent=TRUE)%>%class%>%
                         grepl(., "try-error")%>%any, "0", t1a[,i, j])
        t2ai <- ifelse(try(t2a[,i, j], silent=TRUE)%>%class%>%
                         grepl(., "try-error")%>%any, "0", t2a[,i, j])
        tai<- c(t1ai, t2ai)
        res[,i, j]<- paste0(tai, collapse=" + ")
      }
    }
  }
  res
}

#' \code{convert_attr2exp} converts symbolic attribute of derivative into expression object.
#' @param obj Symbolic expression of gradient or hessian
#' @return combine expression of derivatives
#' @examples
#' eq1 <- parse(text="2*(log(sin(x)/log(x)))+x^4*log(x)+cos(y+x)")
#' tt1 <- deriv(eq1, c("x", "y"), hessian=TRUE)
#' r1 <- convert_attr2exp(extract_attr_deriv(tt1, "grad"))
#' r2 <- convert_attr2exp(extract_attr_deriv(tt1, "hessian"))
#' @export
convert_attr2exp<-function(obj)
{
  . <- NULL
  dobj <- dim(obj)
  varl <- colnames(obj)
  if(length(dobj)==2){
    res <- paste0("    .grad <- array(0, c(1, ",dobj[2],"), list(NULL, c(", paste0("\"",
                                                          colnames(obj),"\"", collapse=" ,"),")))")
    for(i in varl){
      res<-c(res, paste0("    .grad[, \"",i, "\"] <- ", obj[, i] ))
    }
    res <- c(res, "    .grad")
  }else{
    res <- paste0("    .hessian <- array(0, c(1, ",dobj[2],", ",dobj[3] ,"), list(NULL, c(",
                  paste0("\"",colnames(obj),"\"", collapse=" ,"),"), c(",paste0("\"",
                                                          colnames(obj),"\"", collapse=" ,"),")))")
    for(i in varl){
      for(j in varl){
        res<-c(res, paste0("    .hessian[, \"",i, "\", \"", j ,"\"] <- ", obj[,i, j] ))
      }
    }
    res <- c(res, "    .hessian")
  }
  parse(text=c("{", res,"}"))
}


#' \code{replace_par} replaces text with other text.
#' @param iter Integer, which line of repdat is used, if 1 iteratively all will be replaced,
#' if >1 only the i-th parameter.
#' @param repdat data.frame with columns "old" and "new"
#' @param biogu string which will be modified
#' @return Modified string.
#' @examples
#' eq_c <- c("Tw ~ ((((PH) + (tw)) * (ta - Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) -(1 + (PH))) +
#'  sqrt((((PH) + (tw)) * (ta - Tc + 2) + (1 +(tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 - 4 *
#'   (1 + (PH) + (tw)) *(-(PH) * (ta - Tc + 2) + (1 - (tw) * (ta - Tc + 2)) * (2/w -Ec/w))))/(2 *
#'   (1 + (PH) + (tw)))",
#' "Tf1 ~ (th1) * (ta - (((((PH) + (tw)) * (ta - Tc + 2) + (1 + (tw)) *(Ec/w - 2/w) - (1 + (PH))) +
#'  sqrt((((PH) + (tw)) * (ta -Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 - 4 *(1 + (PH) +
#'   (tw)) * (-(PH) * (ta - Tc + 2) + (1 - (tw) *(ta - Tc + 2)) * (2/w - Ec/w))))/(2 * (1 + (PH) +
#'   (tw)))) -Tc + 2) - 1",
#' "Ef1 ~ (ph1)/(PH) * (w * (((((PH) + (tw)) * (ta - Tc + 2) + (1 +(tw)) * (Ec/w - 2/w) -
#' (1 + (PH))) + sqrt((((PH) + (tw)) *(ta - Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 -4 *
#'  (1 + (PH) + (tw)) * (-(PH) * (ta - Tc + 2) + (1 - (tw) *(ta - Tc + 2)) * (2/w - Ec/w))))/(2 *
#'  (1 + (PH) + (tw)))) -Ec + 2) - 1")
#' parl <- c("tw","PH","th1","ph1")
#' parll <- paste0("par[", 1:length(parl), "]")
#' repdat <- data.frame(old=parl, new=parll)
#' replace_par(1, repdat, eq_c)
#' @export
replace_par <- function(iter, repdat, biogu)
{
  . <- NULL
  cl <- class(biogu)
  if(cl%in%"expression"%>%any){
    biogu <- gsub(paste0("(",repdat[iter, "old"],")","(\\D|$)(| )"), paste0(repdat[iter, "new"],
                                                                            "\\2\\3"), biogu)
    biogu[c(1, length(biogu))] <- gsub("expression[(]|[)]$","",biogu[c(1, length(biogu))])
    biogu <- parse(text=biogu)
  }else{
    biogu <- gsub(paste0("(",repdat[iter, "old"],")","(\\D|$)(| )"),
                  paste0(repdat[iter, "new"],"\\2\\3"), biogu)
  }
  iter <- iter + 1
  if(iter <= nrow(repdat)){
    return(replace_par(iter, repdat, biogu))
  }else{
    return(biogu)
  }
}

#' \code{replace_par_wrap} replace text with other text, wrapper of \code{replace_par}.
#' @param repdat data.frame with columns "old" and "new"
#' @param obj character string or vector of strings that will be modified.
#' @return Modified string/vector.
#' @examples
#' eq_c <- c("Tw ~ ((((PH) + (tw)) * (ta - Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) -(1 + (PH))) +
#' sqrt((((PH) + (tw)) * (ta - Tc + 2) + (1 +(tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 - 4 *
#'  (1 + (PH) + (tw)) *(-(PH) * (ta - Tc + 2) + (1 - (tw) * (ta - Tc + 2)) * (2/w -Ec/w))))/(2 *
#'   (1 + (PH) + (tw)))",
#' "Tf1 ~ (th1) * (ta - (((((PH) + (tw)) * (ta - Tc + 2) + (1 + (tw)) *(Ec/w - 2/w) - (1 + (PH))) +
#'  sqrt((((PH) + (tw)) * (ta -Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 - 4 *(1 + (PH) +
#'   (tw)) * (-(PH) * (ta - Tc + 2) + (1 - (tw) *(ta - Tc + 2)) * (2/w - Ec/w))))/(2 * (1 + (PH) +
#'    (tw)))) -Tc + 2) - 1",
#' "Ef1 ~ (ph1)/(PH) * (w * (((((PH) + (tw)) * (ta - Tc + 2) + (1 +(tw)) * (Ec/w - 2/w) -
#' (1 + (PH))) + sqrt((((PH) + (tw)) *(ta - Tc + 2) + (1 + (tw)) * (Ec/w - 2/w) - (1 + (PH)))^2 -4 *
#'  (1 + (PH) + (tw)) * (-(PH) * (ta - Tc + 2) + (1 - (tw) *(ta - Tc + 2)) * (2/w - Ec/w))))/(2 *
#'  (1 + (PH) + (tw)))) -Ec + 2) - 1")
#' parl <- c("tw","PH","th1","ph1")
#' parll <- paste0("par[", 1:length(parl), "]")
#' repdat <- data.frame(old=parl, new=parll)
#' replace_par(2, repdat, eq_c)
#' replace_par_wrap(repdat, eq_c)
#' @export
replace_par_wrap <- function(repdat,obj)
{
  . <- NULL
  repdat[, "old"] <- as.character(repdat[, "old"])
  repdat[, "new"] <- as.character(repdat[, "new"])
  repdat <- repdat[order(-nchar(repdat[,1])),] # start replacing by longest string, B11, B11B
  obj <- replace_par(1, repdat, obj)
  obj
}
