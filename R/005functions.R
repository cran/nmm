#' Hessian and Gradient expressions
#' 
#' Get expression of Hessian and Gradient for discrete choice model with quantile transformation.
#' 
#' @keywords internal
#' @param mnf1 object containing MNL expressions
#' @param smodu1 parameter matrix
#' @param attribute can be "grad" or "hessian"
#' @param ffor discrete choice equations
#' @param cheqs0 errors from cont. equations
#' @return derivative expression of discrete equations, with dependency on base mode
deriv_disc_expr_corr <-  function(mnf1 = NULL, smodu1, attribute = "grad", ffor = ffor,
                                  cheqs0 = cheqs0)
{
  . <- NULL
  pedef_b <- NULL
  const_b <- NULL
  pedef <- NULL
  kk <- length(ffor)
  k1 <- 1:kk
  # here we form character expressions for the gradient and hessian function
  # if not supplied calculate
  if(is.null(mnf1)){
    mnf1 <- MNlogitf(ffor, transform = TRUE, separatenmm = TRUE, weight_disc = TRUE)
  }
  # MLogit formulation
  ss <- paste0('avl_', k1, ' * eps[', k1, ']' )
  sumss <- ss%>%paste0('(',.,')')%>%paste0( . , collapse = ' + ')%>%paste0('(', . ,')')%>%
    gsub('(eps\\[\\d{1,10}\\])', '(\\1)', .)
  #probabilities
  robsse <- paste0( ss, '/(', paste0(ss, collapse = "+"),')')
  #if prob is 0, qnorm return - Inf, if prob=1, it returns Inf
  ife <- paste0('ifelse((probs[',k1,'])==0, 1e-100, ifelse((probs[',k1,
                '])==1,1-1e-16,(probs[',k1,'])))')
  #distribution of quantile transformed and stand. probabilities
  pnorme <- paste0('pnorm((qnorm(ife[',k1,'])- meanqt[',k1,'])/sqrt(covqt[',k1,']))')
  #check if pnomre is 0 or 1, replace with values close to 0 and 1
  pnormeife <- paste0('ifelse((pnorme[',k1,'])==0, 1e-100, ifelse((pnorme[',k1,
                      '])==1,1-1e-16,(pnorme[',k1,'])))')
  obj <- list(ss=ss, sumss=sumss, robsse=robsse, ife=ife, pnorme=pnorme, pnormeife=pnormeife)

  dl <- deriv_disc_corr_help(mnf1=mnf1, smodu1=smodu1,
                                    attribute=attribute,  obj=obj, ffor=ffor, cheqs0 = cheqs0 )

  #for string expressions
  kko <- kk
  fforexpr <- paste0('eps[', k1, '] <- "exp(', ffor, ')"')
  fforexpr <-c( fforexpr, k1%>%paste0('eps[',.,'] <- replace(eps[',.,
                                        '],eps[',.,']%in%Inf, 999999999999 )'))
  probsse <- paste0('probs[',k1,'] <- "', ss, '/(', paste0(ss, collapse = "+"),')"')
  ifee <- paste0('ife[',k1,'] <- "', ife,'"')


  pnormee <- paste0('pnorme[',k1,'] <- "', lapply(dl, '[[', "pnormen")%>%unlist,'"')
  pnormeife <-paste0('pife[',k1,'] <- "', pnormeife,'"')

  term1 <- paste0('term1', k1, ' <- "',lapply(dl, '[[', 'term1'),'"')
  term2 <- paste0('term2', k1, ' <- "',lapply(dl, '[[', 'term2'),'"')
  term3 <- paste0('term3', k1, ' <- "',lapply(dl, '[[', 'term3'),'"')
  term3_b <- paste0('term3_b', k1, ' <- "',lapply(dl, '[[', 'term3_b'),'"') %>% gsub("probs", "ife", .)

  tterm1 <- paste0('tterm1', k1, ' <- "',lapply(dl, '[[', 'tterm1'),'"')
  tterm3 <- paste0('tterm3', k1, ' <- "',lapply(dl, '[[', 'tterm3'),'"')
  tterm3_b <- paste0('tterm3_b', k1, ' <- "',lapply(dl, '[[', 'tterm3_b'),'"') %>% gsub("probs", "ife", .)

  meanq <- lapply(dl, '[[', 'meanq')%>%unlist%>%paste0
  covq <- lapply(dl, '[[', 'covq')%>%unlist%>%paste0

  bd <- lapply(dl, '[[', "bde")
  bd_b <- lapply(dl, '[[', "bde_b")
  dd <- lapply(dl, '[[', "dd")
  dd_b <- lapply(dl, '[[', "dd_b")
  ad <- lapply(dl, '[[', "ad")
  cd <- lapply(dl, '[[', "cd")
  meange <- lapply(dl, '[[', "meange")
  meanhe <- lapply(dl, '[[', "meanhe")
  envr <- environment()
  a <- lapply( c("pedef", "ad", "dd","cd", "meange", "pedef_b", "dd_b"), function(var){
    vare <-  lapply(dl, '[[', var)%>%lapply(.,function(x)paste0('list(', paste0(paste0('"',x, '"'),
                                                                              collapse=", "), ')'))
    vare %<>% paste0(. , collapse=", ")%>% paste0(var, ' <- list(', . ,')')
    assign(var, vare, envir = envr)
  })
  a <- lapply( c("const_b"), function(var){
    vare <-  lapply(dl, '[[', var)%>%paste0('"',., '"')%>%paste0(., collapse=", ")%>%
      paste0(var, ' <- list(', . ,')')
    assign(var, vare, envir = envr)
  })
  a <- lapply(1:length(bd), function(jj){
    vare <-  bd[[jj]]
    vare <- lapply(1:dim(vare)[2],function(x)paste0('list(', paste0(paste0('"',vare[,x, ], '"'),
                                                                    collapse=", "), ')'))
    vare %<>% paste0(. , collapse=", ")%>% paste0(paste0("bd", jj), ' <- list(', . ,')')
    assign(paste0("bd", jj), vare, envir = envr)
  })
  a <- lapply(1:length(bd_b), function(jj){
    vare <-  bd_b[[jj]]
    vare <- lapply(1:dim(vare)[2],function(x)paste0('list(', paste0(paste0('"',vare[,x, ], '"'),
                                                                    collapse=", "), ')'))
    vare %<>% paste0(. , collapse=", ")%>% paste0(paste0("bd_b", jj), ' <- list(', . ,')')
    assign(paste0("bd_b", jj), vare, envir = envr)
  })
  a <- lapply(1:length(meanhe), function(jj){
    vare <-  meanhe[[jj]]
    vare <- lapply(1:dim(vare)[2],function(x)paste0('list(', paste0(paste0('"',vare[,x, ], '"'),
                                                                    collapse=", "), ')'))
    vare %<>% paste0(. , collapse=", ")%>% paste0(paste0("meanhe", jj), ' <- list(', . ,')')
    assign(paste0("meanhe", jj), vare, envir = envr)
  })

  sumsse <- paste0("sumss <- ",sumss)
  obj <- c( fforexpr , sumsse,probsse, ifee, pnormee, pnormeife, term1, term2, term3,term3_b,
            pedef_b,const_b, pedef, meange)

  obje <- obj%>%gsub('"', "", .)
  obje%<>%string2formula%>% gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3", .)%>%
    gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .)

  #1st derivative
  d0 <- k1%>%paste0('derive',. ,' <-  lapply(1:length(par), function(k) (term1',.,')*(chc_',
                      .,')*(avl_',.,')*(wd_',.,')*(1/npaths)*(term2',.,')*( pedef[[', .,
                      ']][[k]]*(term3', ., ')- meange[[', ., ']][[k]]-const_b[[',.,
                      ']]*pedef_b[[', ., ']][[k]]*(term3_b', ., '))/sqrt(covqt[',.,
                      ']))')
  d1 <- paste0('derive <-  list(',paste0('derive', k1, collapse=', '), ')')
  d2 <- c('derive <- lapply(1:length(par), function(z)lapply(derive,"[[", z)%>%do.call("cbind",.)%>%rowSums(., na.rm=FALSE))%>%do.call("cbind", .)',
          'if(!methodopt%in%c("BHHH")) derive <- apply(derive, 2, function(x)sum(x, na.rm=FALSE))')
  m0 <-  k1%>%paste0('meanqt', ., ' <- mean(qnorm(ife[[',.,']]), na.rm=FALSE)')
  s0 <-  k1%>%paste0('covqt', ., ' <- var(qnorm(ife[[',.,']]), na.rm=FALSE)')


  m1 <- paste0('meanqt <- c(',paste0('meanqt',k1, collapse = ", "),')')
  s1 <- paste0('covqt <- c(',paste0('covqt',k1, collapse = ", "),')')
  xg <- paste0('par', 1:nrow(smodu1))%>%paste0('"', ., '"')%>%
    paste0(., collapse=", ")%>%paste0('c(', ., ')')
  g0 <- paste0('grada <- array(derive[-1], c(1,length(par)-1), list(NULL, ',xg,'))')

  #Hessian
  terms <- c("term2","pedef","tterm1", "tterm3")
  derterms <-c("ad", "bd", "cd", "dd")
  a <- lapply(1:length(terms), function(j){
    tt1 <- terms[-j]
    ttt <- grep('term', tt1, value=TRUE)
    ttp <- grep('term', tt1, value=TRUE, invert=TRUE)
    ttt%<>%paste0(., sort(rep(k1, length(ttt))))%>%matrix(.,ncol=length(ttt), byrow=TRUE)
    ttt <- sapply(1:nrow(ttt), function(z) paste0(ttt[z,], collapse=' * '))
    if(length(ttp)==0){
      var <- paste0('r',j, k1,
                    ' <- lapply(1:length(par), function(k)lapply(1:length(par), function(j) (',
                    derterms[j], k1, '[[k]][[j]]',')  * (',ttt,')))')
    }else{
      var <- paste0('r',j, k1,
                    ' <- lapply(1:length(par), function(k)lapply(1:length(par), function(j) (',
                    derterms[j],'[[',k1,']][[k]]',') * (', ttp,'[[',k1,']][[j]]) * (',ttt,')))')
    }
    assign(paste0('r', j), var, envir =  envr)
  })
  terms <- c("term2","pedef_b","tterm1", "tterm3_b")
  derterms <-c("ad", "bd_b", "cd", "dd_b")
  a <- lapply(1:length(terms), function(j){
    tt1 <- terms[-j]
    ttt <- grep('term', tt1, value=TRUE)
    ttp <- grep('term', tt1, value=TRUE, invert=TRUE)
    ttt%<>%paste0(., sort(rep(k1, length(ttt))))%>%matrix(.,ncol=length(ttt), byrow=TRUE)
    ttt <- sapply(1:nrow(ttt), function(z) paste0(ttt[z,], collapse=' * '))
    if(length(ttp)==0){
      var <- paste0('rb',j, k1,
                    ' <- lapply(1:length(par), function(k)lapply(1:length(par), function(j) (',
                    derterms[j], k1, '[[k]][[j]]',')  * (',ttt,')))')
    }else{
      var <- paste0('rb',j, k1,
                    ' <- lapply(1:length(par), function(k)lapply(1:length(par), function(j) (',
                    derterms[j],'[[',k1,']][[k]]',') * (', ttp,'[[',k1,']][[j]]) * (',ttt,')))')
    }
    assign(paste0('rb', j), var, envir =  envr)
  })

  terms <- c("term2", "meange","tterm1")
  derterms <-c("ad", "meanhe", "cd")
  a <- lapply(1:length(terms), function(j){
    tt1 <- terms[-j]
    ttt <- grep('term', tt1, value=TRUE)
    ttp <- grep('term', tt1, value=TRUE, invert=TRUE)
    ttt%<>%paste0(., sort(rep(k1, length(ttt))))%>%matrix(.,ncol=length(ttt), byrow=TRUE)
    ttt <- sapply(1:nrow(ttt), function(z) paste0(ttt[z,], collapse=' * '))
    if(length(ttp)==0){
      var <- paste0('rp',j, k1,
                    ' <- lapply(1:length(par), function(k)lapply(1:length(par), function(j) (',
                    derterms[j], k1, '[[k]][[j]]',')  * (',ttt,')))')
    }else{
      var <- paste0('rp',j, k1,
                    ' <- lapply(1:length(par), function(k)lapply(1:length(par), function(j) (',
                    derterms[j],'[[',k1,']][[k]]',') * (', ttp,'[[',k1,']][[j]]) * (',ttt,')))')
    }
    assign(paste0('rp', j), var, envir =  envr)
  })


  cd %<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)

  obje <- c("par <- c(0, par)","eps <- list()", "ife <- c()", "probs <- c()","pnorme <- c()",
            "pife <- c()", obje[1:(grep("^pnorme", obje)[1]-1)],m0,m1, s0, s1,
            obje[(grep("^pnorme", obje)[1]):length(obje)], d0, d1, cd, d2)
  addt <- ""

  if(attribute=="grad"){
    obje <- c("{", grep('^cd', obje, value=TRUE, invert = TRUE),g0,"res <- unlist(grada)", addt,
              "res","}")
  }

  if(attribute=="hessian"){
    hh0 <- k1%>%paste0('hessian',.,
                  ' <- lapply(1:length(par), function(k) lapply(1:length(par),  function(j) ((r1',
                         .,'[[k]][[j]] + r2',.,'[[k]][[j]]-r3',.,'[[k]][[j]]-r4',.,
                         '[[k]][[j]])/ifelse((tterm1',.,'*tterm3',.,')^2==0, 1.e-100, tterm1',.,
                         '*tterm3',.,
                         ')^2/sqrt(covqt[[',.,']])-
                         (rb1',.,'[[k]][[j]] + rb2',.,'[[k]][[j]]-rb3',.,'[[k]][[j]]-rb4',.,
                         '[[k]][[j]])/ifelse((tterm1',.,'*tterm3_b',.,')^2==0, 1.e-100, tterm1',.,
                         '*tterm3_b',.,
                         ')^2/sqrt(covqt[[',.,']])-
                         (rp1',.,'[[k]][[j]] + rp2',.,'[[k]][[j]]-rp3',.,
                         '[[k]][[j]])/ifelse((tterm1',.,')^2==0, 1.e-100, tterm1',.,
                         ')^2/sqrt(covqt[[',.,']]))*chc_',.,'*avl_',.,'*wd_',.,'*(1/npaths)))')

    objr <- paste0(rep(c("r", "rb", "rp"), rep(4, 3)), 1:4)
    objr <- objr[objr!="rp4"]
    objri <- paste0(rep(objr, rep(kk, length(objr))), 1:kk) %>% paste0("\\[\\[k\\]\\]\\[\\[j\\]\\]")
    objre <- eval(objr %>% paste0(., collapse = ",") %>% paste0("c(", ., ")") %>% parse(text=.))
    objre %<>% gsub(".*function\\(j\\)","", .) %>% trimws %>% gsub("\\)\\)$", "",.) %>% paste0("(", ., ")")
    repdat <- data.frame(new=objre, old=objri,stringsAsFactors = FALSE)
    hh0 <- replace_par_wrap(repdat, hh0)
    
    hh1 <- paste0('hessian <-  list(',paste0('hessian', k1, collapse=', '), ')')
    hh2 <- c('hessian <- lapply(1:length(par), function(k) lapply(1:length(par),
             function(j) lapply(hessian, "[[", k) %>% lapply(.,
             "[[", j) %>% do.call("cbind", .) %>% rowSums(., na.rm=FALSE)))',
             'hessian <- lapply(hessian, function(x) x%>%do.call("cbind", .))%>%abind(., along=3)',
             'res <- hessian',
             'if (!methodopt %in%c("BHHH")) res <- sapply(1:(length(par)), function(k) sapply(1:(length(par)), function(j) hessian[, k, j] %>% sum(., na.rm=FALSE)))')

    ad %<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)
    dd %<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)
    dd_b %<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)

    tterm1%<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)
    tterm3%<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)
    tterm3_b%<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)
    bdd <- eval(parse(text=paste0("c(",paste0("bd", k1, collapse = ","),")")))%>%
      gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3", .  )
    bdd_b <- eval(parse(text=paste0("c(",paste0("bd_b", k1, collapse = ","),")")))%>%
      gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3", .  )
    meanhed <- eval(parse(text=paste0("c(",paste0("meanhe", k1, collapse = ","),")")))%>%
      gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3", .  )
    obje <- c("{",obje, tterm1, tterm3, tterm3_b, ad, dd,dd_b, bdd,bdd_b, meanhed,
              hh0, hh1, hh2, "res <- res[,-1,-1]",addt, "res", "}")
  }
  obje
}

#' @keywords internal
deriv_disc_corr_help <- function(mnf1, smodu1,  kk,  data="data", attribute="grad", obj,
                                 ffor=ffor, cheqs0=cheqs0, ...)
{
  . <- NULL
  kk <- length(ffor)
  k1 <- 1:kk
  ss <- obj$ss # indirect utility * avl
  sumss <-  obj$sumss # denom in prob. expression
  robsse <-  obj$robsse # probab.
  ife <-  obj$ife #ifelse(probab)
  pnorme <-  obj$pnorme # distribution of transformed probab.
  pnormeife <-  obj$pnormeife # ifelse(pnorme)
  neqt <- length(cheqs0)
  #if expressions of conditional means are supplied
  if(mnf1%>%names%>%grepl('meanexp', .)%>%any){
    expr_meanq <- mnf1$meanexpr
    expr_meanq%<>%string2formula
    expr_mean_nom <- mnf1$nominl%>%string2formula
    expr_mean_denom <- mnf1$denomil%>%string2formula
  }
  namevec <- smodu1[,"parlp"]
  #for each equation get expressions of derivatives
  to_b <- mnf1$expreteso[[kk]] #distribution of transformed prob.
  te_b <- to_b%>%formula2string()%>%gsub("\\_(\\d)\\_", "\\1", .)
  tq_b <- mnf1$probte[[kk]]# quantile transformation
  tqe_b <- tq_b%>%formula2string()%>%gsub("\\_(\\d)\\_", "\\1", .)
  term3_b <- paste0("1/dnorm(", tqe_b, ")") #1/dnorm 1/t3
  tterm3_b <- paste0("dnorm(", tqe_b, ")") #dnorm t3
  po_b <- mnf1$probs[[kk]] #probability
  pe_b <- po_b%>%formula2string() #string of po
  ped_b <- deriv(parse(text=paste0(pe_b)), namevec= namevec, hessian=FALSE) #1st derivatives of prob.
  pede_b <- extract_attr_deriv(ped_b, attribute = "grad") #string expression
  pedef_b <- string2formula(pede_b)%>%gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .) #eval expr t4 from appendix
  bd_b <- deriv(parse(text=paste0(pe_b)), namevec= namevec, hessian=TRUE)#hessian of of prob
  bde_b <- extract_attr_deriv(bd_b, attribute = "hessian")#string expression
  bde_b <- string2formula(bde_b)%>%gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .)#eval expression

  probt_b <- mnf1$probt[[kk]]
  dd_b <- deriv(parse(text='dnorm(Func)'),namevec='Func')%>%
    extract_attr_deriv(., attribute = 'grad')%>%gsub('Func', probt_b, .)
  dd_b%<>%paste0('(',.,') * (', term3_b, ') * (', pedef_b,')')

  dl <- lapply(k1, function(k){
    #conditional mean part(constant) nomin/denomin
    nom_b <- expr_mean_nom%>%grep(paste0("nomin",k,"_",length(cheqs0)+1,"_"), ., value=TRUE,
                                  fixed=TRUE)%>%gsub(".*<-", "", .)
    denom_b <- expr_mean_denom%>%grep("denom0", ., value=TRUE, fixed=TRUE)%>%gsub(".*<-", "", .)
    const_b <- paste0('(', nom_b, ')/(', denom_b,')')
    if(k==kk){
      const_b <- "0"
    }
    to <- mnf1$expreteso[[k]] #distribution of transformed prob.
    te <- to%>%formula2string()%>%gsub("\\_(\\d)\\_", "\\1", .)
    tq <- mnf1$probte[[k]]# quantile transformation
    tqe <- tq%>%formula2string()%>%gsub("\\_(\\d)\\_", "\\1", .)

    #expressionf for conditional mean
    cmeane <- if(mnf1%>%names%>%grepl('meanexp', .)%>%any)mnf1$meanexpr[[k]]%>%
      gsub("\\[(\\d)x(\\d)\\]", "[\\1,\\2]",.)%>%formula2string else 'par1/par2'

    po <- mnf1$probs[[k]] #probability
    pe <- po%>%formula2string() #string of po
    probt <- mnf1$probt[[k]] #qnorm(ifelse(po))
    probtee <- pe%>%paste0("qnorm(",., ")") #string of probt

    term1 <- paste0("1/(", to, ")") #1/pnorm 1/t1 from appendix
    tterm1 <- paste0("(", to, ")") # pnorm t1 from appendix
    term2 <- gsub("pnorm", "dnorm", to) # dnorm t2
    term3 <- paste0("1/dnorm(", tqe, ")") #1/dnorm 1/t3
    tterm3 <- paste0("dnorm(", tqe, ")") #dnorm t3

    #get mean string
    meanq <- paste0("meanqt[", k, "]")
    covq <- paste0("covqt[", k, "]")

    #derivatives of prob.
    ped <- deriv(parse(text=pe), namevec= namevec, hessian=FALSE) #1st derivatives of prob.
    pede <- extract_attr_deriv(ped, attribute = "grad") #string expression
    pedef <- string2formula(pede)%>%gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .) #eval expr t4 from appendix
    bd <- deriv(parse(text=pe), namevec= namevec, hessian=TRUE)#hessian of of prob
    bde <- extract_attr_deriv(bd, attribute = "hessian")#string expression
    bde <- string2formula(bde)%>%gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .)#eval expression

    #cond. mean derivative, we need to divide into 2 parts, one for continuous other for discrete
    if(k==kk){
      contp <- cmeane
    }else{
      indk <- paste0("nomin",k,"_",length(cheqs0)+1,"_")
      parts <- strsplit(cmeane, indk)%>%unlist
      contp <- parts[1]%>%gsub("\\+\\($","", .)%>%paste0(.,")/denom0")
    }

    meang <- deriv(parse(text=contp), namevec= namevec, hessian=TRUE)
    mact <-.%>%string2formula%>%gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .)
    meange <- extract_attr_deriv(meang, attribute = "grad")%>%mact
    meanhe <- extract_attr_deriv(meang, attribute = "hessian")%>%mact
    if(!mnf1%>%names%>%grepl('meanexp', .)%>%any){
      meange[1,] <- ' 0'
      meanhe[,1:2,] <- ' 0'
    }

    pnormen <- pnorme

    #########################second derivative see discexplain.Rmd
    #ad=deriv dnorm()
    #bd=deriv pede
    #cd= 1 deriv
    #dd= deriv dnorm
    argad <- term2 %>% gsub('^dnorm\\(|\\)$', '', .)%>%paste0('(', ., ')')
    ad <- deriv(parse(text='dnorm(Func)'),namevec='Func')%>%
      extract_attr_deriv(., attribute = 'grad')%>%gsub('Func', argad, .)
    ad%<>%paste0('(', . ,') * ((', term3, ') * (', pedef, ') - (', meange ,
                 '))*1/', 'sqrt(covqt[',k,'])' )
    cd <- paste0('(', term2 ,') * ((', term3, ') * (', pedef, ') - (',meange,
                 '))*1/', 'sqrt(covqt[',k,'])' )
    
    dd <- deriv(parse(text='dnorm(Func)'),namevec='Func')%>%
      extract_attr_deriv(., attribute = 'grad')%>%gsub('Func', probt, .)
    dd%<>%paste0('(',.,') * (', term3, ') * (', pedef,')')

    #simplify the expressions
    environment(fff1) <- environment()
    term1 %<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)
    term2 %<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)
    term3 %<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)
    pedef %<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)

    ad%<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)
    cd%<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)
    dd%<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)
    bde%<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)
    tterm1 %<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)
    tterm3 %<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)
    ifeact <- .%>%gsub("probs", "ife", .) 
    term3_b %<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss) %>% ifeact
    tterm3_b %<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)%>% ifeact
    pedef_b%<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)%>% ifeact
    bde_b%<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)%>% ifeact
    dd_b%<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)%>% ifeact
    const_b%<>%fff1(., probtee,  ffor, smodu1, robsse, ife, pnormen, sumss)%>% ifeact

    pnormen <- pnormen[k]
    if(k==kk){
      term3_b <-0
      tterm3_b <-0
      pedef_b <- array(0, dim = dim(pedef_b))
      bde_b <- array(0, dim = dim(bde_b))
      dd_b <-  rep(0, length(dd_b))
    }
    res <- list(term1=term1, term2=term2, term3=term3, pedef=pedef, meanq=meanq,  covq= covq,
                pnormen=pnormen,ad=ad, bde=bde,cd=cd,  dd=dd, tterm3=tterm3, tterm1=tterm1,
                meange=meange, meanhe=meanhe,
                term3_b=term3_b, tterm3_b=tterm3_b, pedef_b=pedef_b, bde_b=bde_b, dd_b=dd_b,
                const_b =const_b)
    res
  })
  return(dl)
}


