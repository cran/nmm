#' Get expression of Hessian and Gradient for discrete choice model with quantile transformation
#' @keywords internal
#' @param mnf1 object containing MNL expressions
#' @param smodu1 parameter matrix
#' @param attribute can be "grad" or "hessian"
#' @param ffor discrete choice equations
#' @param MNtypef logit or dogit
#' @return derivative expression of discrete equations
deriv_disc_expr <-  function(mnf1=NULL, smodu1, attribute="grad", ffor=ffor, 
                             MNtypef="logit")
{
  . <- NULL
  pedef <- NULL
  kk <- length(ffor)
  #here we form character expressions for the gradient and hessian function
  #if not summplied calculate
  if(MNtypef=="logit"){
    ppardogit <- NULL
  }
  if(MNtypef=="dogit"){
    ndd <- ffor %>% gsub("(\\[\\d{1,10}\\])", "XX\\1XX", .) %>% strsplit("XX") %>% unlist %>% unique %>% 
      grep("\\[", ., value = TRUE) %>% gsub("\\[|\\]", "", .) %>% as.numeric() %>% max
    ppardogit <- c(paste0('par[',(ndd+1):(ndd+length(ffor)),']'))
  }
  
  if(is.null(mnf1)){
    if(MNtypef=="logit"){
      mnf1 <- MNlogitf(ffor, transform = TRUE, separatenmm = TRUE, weight_disc = TRUE)
    }
    if(MNtypef=="dogit"){
      mnf1 <- MNdogitf(ffor, transform = TRUE, separatenmm = TRUE, weight_disc = TRUE,
                       ppardogit=ppardogit)
    }
  }

  #MLogit formulation
  ss <- paste0('avl_', 1:kk, ' * eps[', 1:kk, ']' )
  sumss <- ss%>%paste0('(',.,')')%>%paste0( . , collapse = ' + ')%>%paste0('(', . ,')')%>%
    gsub('(eps\\[\\d{1,10}\\])', '(\\1)', .)
  #probabilities
  robsse <- paste0( ss, '/(', paste0(ss, collapse = "+"),')')
  

  #if prob is 0, qnorm return - Inf, if prob=1, it returns Inf
  ife <- paste0('ifelse((probs[',1:kk,'])==0, 1e-100, ifelse((probs[',1:kk,
                '])==1,1-1e-16,(probs[',1:kk,'])))')
  #distribution of quantile transformed and stand. probabilities
  pnorme <- paste0('pnorm((qnorm(ife[',1:kk,'])- meanqt[',1:kk,'])/sqrt(covqt[',1:kk,']))')
  #check if pnomre is 0 or 1, replace with values close to 0 and 1
  pnormeife <- paste0('ifelse((pnorme[',1:kk,'])==0, 1e-100, ifelse((pnorme[',1:kk,
                      '])==1,1-1e-16,(pnorme[',1:kk,'])))')

  if(MNtypef=="logit"){
    ss1 <- NULL
    robsse1 <- NULL
    obj <- list(ss=ss, sumss=sumss, robsse=robsse, ife=ife, pnorme=pnorme, pnormeife=pnormeife)
    
    dl <- deriv_disc_help(mnf1 = mnf1, smodu1 = smodu1,
                          attribute = attribute,  obj = obj, ffor = ffor, MNtypef = MNtypef)
  }
  if(MNtypef=="dogit"){
    sumss1 <- ss%>%paste0( . , collapse = ' + ')%>%paste0('(', . ,')')%>%
      gsub('(eps\\[\\d{1,10}\\])', '(\\1)', .)
    #ss1 <- paste0('eps[', 1:kk, '] + ', ppardogit, " * (", sumss1,")")
    ss1 <- paste0(ss, ' + avl_', 1:kk,' * ', ppardogit, " * (", sumss1,")")
    pjis <- ppardogit %>%
      paste0(., collapse = ' + ') %>% paste0('(1 + ', ., ')')
    robsse1 <- paste0('(', ss1, ")/(", pjis, "*",sumss1, ")")
    obj <- list(ss=ss1, sumss=sumss1, robsse=robsse1, ife=ife, pnorme=pnorme, pnormeife=pnormeife)
    dl <- deriv_disc_help(mnf1=mnf1, smodu1=smodu1,
                          attribute=attribute,  obj=obj, ffor=ffor, MNtypef=MNtypef)
    sumss <- sumss1
    robsse <- robsse1
    ss <- ss1
  }
  
  #for string expressions
  fforexpr <- paste0('eps[', 1:kk, '] <- "exp(', ffor, ')"')
  fforexpr <-c( fforexpr, 1:kk%>%paste0('eps[',.,
                                        '] <- replace(eps[',.,'],eps[',.,']%in%Inf, 999999999999 )'))
  #probsse <- paste0('probs[',1:kk,'] <- "', ss, '/(', paste0(ss, collapse = "+"),')"')
  probsse <- paste0('probs[',1:kk,'] <- "',robsse,'')
  ifee <- paste0('ife[',1:kk,'] <- "', ife,'"')


  pnormee <- paste0('pnorme[',1:kk,'] <- "', lapply(dl, '[[', "pnormen")%>%unlist,'"')
  pnormeife <-paste0('pife[',1:kk,'] <- "', pnormeife,'"')

  term1 <- paste0('term1', 1:kk, ' <- "',lapply(dl, '[[', 'term1'),'"')
  term2 <- paste0('term2', 1:kk, ' <- "',lapply(dl, '[[', 'term2'),'"')
  term3 <- paste0('term3', 1:kk, ' <- "',lapply(dl, '[[', 'term3'),'"')

  tterm1 <- paste0('tterm1', 1:kk, ' <- "',lapply(dl, '[[', 'tterm1'),'"')
  tterm3 <- paste0('tterm3', 1:kk, ' <- "',lapply(dl, '[[', 'tterm3'),'"')

  meanq <- lapply(dl, '[[', 'meanq')%>%unlist%>%paste0
  covq <- lapply(dl, '[[', 'covq')%>%unlist%>%paste0

  #ggr <- expand.grid( 1:nrow(smodu1), 1:kk)
  bd <- lapply(dl, '[[', "bde")
  dd <- lapply(dl, '[[', "dd")
  ad <- lapply(dl, '[[', "ad")
  cd <- lapply(dl, '[[', "cd")
  meange <- lapply(dl, '[[', "meange")
  meanhe <- lapply(dl, '[[', "meanhe")
  envr <- environment()
  a <- lapply( c("pedef", "ad", "dd","cd", "meange"), function(var){
    vare <-  lapply(dl, '[[', var)%>%lapply(.,function(x)paste0('list(', paste0(paste0('"',x,
                                                                        '"'), collapse=", "), ')'))
    vare %<>% paste0(. , collapse=", ")%>% paste0(var, ' <- list(', . ,')')
    assign(var, vare, envir = envr)
  })
  a <- lapply(1:length(bd), function(jj){
    vare <-  bd[[jj]]
    vare <- lapply(1:dim(vare)[2],function(x)paste0('list(', paste0(paste0('"',vare[,x, ], '"'),
                                                                  collapse=", "), ')'))
    vare %<>% paste0(. , collapse=", ")%>% paste0(paste0("bd", jj), ' <- list(', . ,')')
    assign(paste0("bd", jj), vare, envir = envr)
  })
  a <- lapply(1:length(meanhe), function(jj){
    vare <-  meanhe[[jj]]
    vare <- lapply(1:dim(vare)[2],function(x)paste0('list(', paste0(paste0('"',vare[,x, ], '"'),
                                                                    collapse=", "), ')'))
    vare %<>% paste0(. , collapse=", ")%>% paste0(paste0("meanhe", jj), ' <- list(', . ,')')
    assign(paste0("meanhe", jj), vare, envir = envr)
  })

  sumsse <- paste0("sumss <- ",sumss)
  obj <- c( fforexpr , sumsse,probsse, ifee, pnormee, pnormeife, term1, term2, term3,
            pedef, meange)

  obje <- obj%>%gsub('"', "", .)
  obje%<>%string2formula%>% gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3", .)%>%
    gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .)

  #1st derivative
  d0 <- 1:kk%>%paste0('derive',. ,' <-  lapply(1:length(par), function(k) (term1',.,')*(chc_',.,
                      ')*(avl_',.,')*(wd_',.,')*(1/npaths)*(term2',.,')*( pedef[[', .,
                      ']][[k]]*(term3', ., ')- meange[[', ., ']][[k]])/sqrt(covqt[',.,']))')
  d1 <- paste0('derive <-  list(',paste0('derive', 1:kk, collapse=', '), ')')
  d2 <- c('derive <- lapply(1:length(par), function(z)lapply(derive,"[[", z)%>%do.call("cbind",.)%>%rowSums(., na.rm=FALSE))%>%do.call("cbind", .)',
          'if(!methodopt%in%c("BHHH")) derive <- apply(derive, 2, function(x)sum(x, na.rm=FALSE))')
  m0 <-  1:kk%>%paste0('meanqt', ., ' <- mean(qnorm(ife[[',.,']]), na.rm=FALSE)')
  s0 <-  1:kk%>%paste0('covqt', ., ' <- var(qnorm(ife[[',.,']]), na.rm=FALSE)')


  m1 <- paste0('meanqt <- c(',paste0('meanqt',1:kk, collapse = ", "),')')
  s1 <- paste0('covqt <- c(',paste0('covqt',1:kk, collapse = ", "),')')
  xg <- paste0('par', 1:nrow(smodu1))%>%paste0('"', ., '"')%>%paste0(., collapse=", ")%>%
    paste0('c(', ., ')')
  g0 <- paste0('grada <- array(derive[-1], c(1,length(par)-1), list(NULL, ',xg,'))')

  #Hessian
  terms <- c("term2","pedef","tterm1", "tterm3")
  derterms <-c("ad", "bd", "cd", "dd")
  a <- lapply(1:length(terms), function(j){
    tt1 <- terms[-j]
    ttt <- grep('term', tt1, value=TRUE)
    ttp <- grep('term', tt1, value=TRUE, invert=TRUE)
    ttt%<>%paste0(., sort(rep(1:kk, length(ttt))))%>%matrix(.,ncol=length(ttt), byrow=TRUE)
    ttt <- sapply(1:nrow(ttt), function(z) paste0(ttt[z,], collapse=' * '))
    if(length(ttp)==0){
      var <- paste0('r',j, 1:kk,
                    ' <- lapply(1:length(par), function(k)lapply(1:length(par), function(j) (',
                    derterms[j], 1:kk, '[[k]][[j]]',')  * (',ttt,')))')
    }else{
      var <- paste0('r',j, 1:kk,
                    ' <- lapply(1:length(par), function(k)lapply(1:length(par), function(j) (',
                    derterms[j],'[[',1:kk,']][[k]]',') * (', ttp,'[[',1:kk,']][[j]]) * (',ttt,')))')
    }
    assign(paste0('r', j), var, envir =  envr)
  })

  terms <- c("term2", "meange","tterm1")
  derterms <-c("ad", "meanhe", "cd")
  a <- lapply(1:length(terms), function(j){
    tt1 <- terms[-j]
    ttt <- grep('term', tt1, value=TRUE)
    ttp <- grep('term', tt1, value=TRUE, invert=TRUE)
    ttt%<>%paste0(., sort(rep(1:kk, length(ttt))))%>%matrix(.,ncol=length(ttt), byrow=TRUE)
    ttt <- sapply(1:nrow(ttt), function(z) paste0(ttt[z,], collapse=' * '))
    if(length(ttp)==0){
      var <- paste0('rp',j, 1:kk,
                    ' <- lapply(1:length(par), function(k)lapply(1:length(par), function(j) (',
                    derterms[j], 1:kk, '[[k]][[j]]',')  * (',ttt,')))')
    }else{
      var <- paste0('rp',j, 1:kk,
                    ' <- lapply(1:length(par), function(k)lapply(1:length(par), function(j) (',
                    derterms[j],'[[',1:kk,']][[k]]',') * (', ttp,'[[',1:kk,']][[j]]) * (',ttt,')))')
    }
    assign(paste0('rp', j), var, envir =  envr)
  })


  cd %<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)

  obje <- c("par <- c(0, par)","eps <- list()", "ife <- c()", "probs <- c()","pnorme <- c()","pife <- c()",
            "term0 <- c()",
            obje[1:(grep("^pnorme", obje)[1]-1)],m0,m1, s0, s1,
            obje[(grep("^pnorme", obje)[1]):length(obje)],
            d0, d1, cd, d2)
  addt <- ""

  if(attribute=="grad"){
    obje <- c("{", grep('^cd', obje, value=TRUE, invert = TRUE),g0,
              "res <- unlist(grada)", addt, "res","}")
  }

  if(attribute=="hessian"){
    hh0 <- 1:kk%>%paste0('hessian',.,' <- lapply(1:length(par), function(k) lapply(1:length(par),  function(j) ((r1',
                         .,'[[k]][[j]] + r2',.,'[[k]][[j]]-r3',.,'[[k]][[j]]-r4',.,
                         '[[k]][[j]])/ifelse((tterm1',.,'*tterm3',.,')^2==0, 1.e-100, tterm1',.,'*tterm3',.,
                         ')^2/sqrt(covqt[[',.,']])-(rp1',.,'[[k]][[j]] + rp2',.,'[[k]][[j]]-rp3',.,
                         '[[k]][[j]])/ifelse((tterm1',.,')^2==0, 1.e-100, tterm1',.,
                         ')^2/sqrt(covqt[[',.,']]))*chc_',.,'*avl_',.,'*wd_',.,'*(1/npaths)))')
    objr <- paste0(rep(c("r", "rp"), rep(4, 2)), 1:4)
    objr <- objr[objr!="rp4"]
    objri <- paste0(rep(objr, rep(kk, length(objr))), 1:kk) %>% paste0("\\[\\[k\\]\\]\\[\\[j\\]\\]")
    objre <- eval(objr %>% paste0(., collapse = ",") %>% paste0("c(", ., ")") %>% parse(text=.))
    objre %<>% gsub(".*function\\(j\\)","", .) %>% trimws %>% gsub("\\)\\)$", "",.) %>% paste0("(", ., ")")
    repdat <- data.frame(new=objre, old=objri,stringsAsFactors = FALSE)
    hh0 <- replace_par_wrap(repdat, hh0)
    
    
    hh1 <- c(paste0('hessian <-  list(',paste0('hessian', 1:kk, collapse=', '), ')'),
             paste0('rm(list=c(',paste0(paste0('"hessian', 1:kk, '"'),collapse=', '), '))'))
    hh2 <- c('hessian <- lapply(1:length(par), function(k) lapply(1:length(par),
             function(j) lapply(hessian, "[[", k) %>% lapply(.,
             "[[", j) %>% do.call("cbind", .) %>% rowSums(., na.rm=FALSE)))',
             'hessian <- lapply(hessian, function(x) x%>%do.call("cbind", .))%>%abind(., along=3)',
             'res <- hessian',
             'if (!methodopt %in%c("BHHH")) res <- sapply(1:(length(par)), function(k) sapply(1:(length(par)), function(j) hessian[, k, j] %>% sum(., na.rm=FALSE)))')

    ad %<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)
    dd %<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)

    tterm1%<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)
    tterm3%<>%gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3",.)
    bdd <- eval(parse(text=paste0("c(",paste0("bd", 1:kk, collapse = ","),")")))%>%
      gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3", .  )
    meanhed <- eval(parse(text=paste0("c(",paste0("meanhe", 1:kk, collapse = ","),")")))%>%
      gsub('"', "", .)%>%gsub("(\\[)(\\d{1,10})(\\])", "\\1[\\2]\\3", .  )
    obje <- c("{",obje, tterm1, tterm3, ad, dd, bdd, meanhed,
              #eval(parse(text=paste0("c(",paste0("r", 1:4, collapse = ","),")"))),
              #eval(parse(text=paste0("c(",paste0("rp", 1:3, collapse = ","),")"))),
              hh0, hh1, hh2, "res <- res[,-1,-1]",addt, "res", "}")
  }
  obje
}

#' @keywords internal
fff1 <- function(x, probtee,  ffor, smodu1, robsse, ife, pnormen, sumss){
  . <- NULL
  kk <- length(ffor)
  repdat <- unique(smodu1[,c("parlp","parnl")])
  names(repdat) <- c("old", "new")
  x <- replace_par_wrap(repdat, x)

  probtee <- replace_par_wrap(repdat, probtee)

  repl <- gsub("qnorm\\(|\\)$", "", probtee)
  repll <- paste0("qnorm(ifelse((", repl,")==0, 1e-100, ifelse((", repl,")==1,1-1e-16,(", repl,
                  "))))")
  x <- gsub(probtee, repll, x, fixed = TRUE)
  for(j in 1:kk){
    x <- gsub(paste0("exp(",ffor[j],")"), paste0("eps[",j,"]"), x, fixed=TRUE)
  }
  for(j in 1:kk){
    x <- gsub(robsse[j], paste0("probs[",j,"]"), x, fixed=TRUE)
  }
  for(j in 1:kk){
    x <- gsub(ife[j], paste0("ife[",j,"]"), x, fixed=TRUE)
  }
  for(j in 1:kk){
    x <- gsub(pnormen[j], paste0("pnorme[",j,"]"), x, fixed=TRUE)
  }
  for(j in 1:kk){
    x <- gsub(paste0("pnorme[",j,"]"), paste0("pife[",j,"]"), x, fixed=TRUE)
  }
  x <- gsub(sumss, paste0("sumss"), x, fixed=TRUE)
  x
}

#' @keywords internal
deriv_disc_help <- function(mnf1, smodu1,  kk,  data="data",
                            attribute="grad", obj, ffor=ffor, MNtypef="logit", ...)
{
  . <- NULL
  kk <- length(ffor)
  ss <- obj$ss # indirect utility * avl
  sumss <-  obj$sumss # denom in prob. expression
  robsse <-  obj$robsse # probab.
  ife <-  obj$ife #ifelse(probab)
  pnorme <-  obj$pnorme # distribution of transformed probab.
  pnormeife <-  obj$pnormeife # ifelse(pnorme)
  #if expressions of conditional means are supplied
  if(mnf1%>%names%>%grepl('meanexp', .)%>%any){
    expr_meanq <- mnf1$meanexpr
    expr_meanq%<>%string2formula
  }
  #for each equation get expressions of derivatives
  dl <- lapply(1:kk, function(k){
    to <- mnf1$expreteso[[k]] #distribution of transformed prob.
    te <- to%>%formula2string()%>%gsub("\\_(\\d)\\_", "\\1", .)
    tq <- mnf1$probte[[k]]# quantile transformation
    tqe <- tq%>%formula2string()%>%gsub("\\_(\\d)\\_", "\\1", .)

    #expressionf for conditional mean
    #cmeane <- mnf1$meanexpr[[k]]
    cmeane <- if(mnf1%>%names%>%grepl('meanexp', .)%>%any)mnf1$meanexpr[[k]]%>%
      gsub("\\[(\\d)x(\\d)\\]", "[\\1,\\2]",.)%>%formula2string else 'par1/par2'

    po <- mnf1$probs[[k]] #probability
    pe <- po%>%formula2string() #string of po
    probt <- mnf1$probt[[k]] #qnorm(ifelse(po))
    probtee <- pe%>%paste0("qnorm(",., ")") #sting of probt
    namevec <- smodu1[,"parlp"]

    term1 <- paste0("1/(", to, ")") #1/pnorm 1/t1 from appendix
    tterm1 <- paste0("(", to, ")") # pnorm t1 from appendix
    term2 <- gsub("pnorm", "dnorm", to) # dnorm t2
    term3 <- paste0("1/dnorm(", tqe, ")") #1/dnorm 1/t3
    tterm3 <- paste0("dnorm(", tqe, ")") #dnorm t3

    #get mean string
    meanq <- paste0("meanqt[", k, "]")
    #gsub(".*(meanqt\\[)(\\d{1,10})(\\]).*", "\\1\\2\\3", to)
    covq <- paste0("covqt[", k, "]")
    #gsub(".*(covqt\\[)(\\d{1,10})(\\]).*", "\\1\\2\\3", to)

    #derivatives of prob.
    ped <- deriv(parse(text=pe), namevec= namevec, hessian=FALSE) #1st derivatives of prob.
    pede <- extract_attr_deriv(ped, attribute = "grad") #string expression
    pedef <- string2formula(pede)%>%gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .) #eval expression t5 from appendix
    bd <- deriv(parse(text=pe), namevec= namevec, hessian=TRUE)#hessian of of prob
    bde <- extract_attr_deriv(bd, attribute = "hessian")#string expression
    bde <- string2formula(bde)%>%gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .)#eval expression

    #cond. mean derivative
    meang <- deriv(parse(text=cmeane), namevec= namevec, hessian=TRUE)
    meange <- extract_attr_deriv(meang, attribute = "grad")%>%string2formula%>%
      gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .)
    meanhe <- extract_attr_deriv(meang, attribute = "hessian")%>%string2formula%>%
      gsub('(nomin\\d)_(\\d)_', '\\1[\\2]', .)
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
    ad%<>%paste0('(', . ,') * ((', term3, ') * (', pedef, ') - (', meange ,'))*1/',
                 'sqrt(covqt[',k,'])' )
    cd <- paste0('(', term2 ,') * ((', term3, ') * (', pedef, ') - (',meange,'))*1/',
                 'sqrt(covqt[',k,'])' )
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

    pnormen <- pnormen[k]
    res <- list(term1=term1, term2=term2, term3=term3, pedef=pedef, meanq=meanq,  covq= covq,
                pnormen=pnormen,ad=ad, bde=bde,cd=cd,  dd=dd, tterm3=tterm3, tterm1=tterm1, meange=meange, meanhe=meanhe)
    res
  })
  return(dl)
}
