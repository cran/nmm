globalVariables(c("dat4cond_mean_cov_expr", "dat4cond_mean_cov_expr",
                  "datmaxle", "expr_ll_norm", "expr_ll_norm_v2", "datmlsem", "atus_ce"))

#' Help function
#' @keywords internal
ff_generate4maxle_p <- function(outl, neqt, fixed_term){
  . <- NULL
  #phi <- paste0('phi(x, mu, var):= 1/(sqrt(var))*exp((-((x-mu)**2)/(2*var)))')# change to log form
  phi <- paste0('phi(x, mu, var):= - log(sqrt(var))-((x-mu)**2)/(2*var)')
  logdist <- c()
  dist <- c()
  fixedp <- if(fixed_term==TRUE)'log(1/(sqrt(2*%pi)))+' else ''
  ta1 <- .%>%gsub('\\^', '**', .)
  for(j in 1:neqt){
    obj <- outl[[j]]$expr%>%ta1
    logdist <- c(logdist, paste0(fixedp,'(phi(eps',j,', ',obj[1],',', obj[2],'))'))
    dist <- c(dist, paste0('phi(eps',j,', ',obj[1],',', obj[2],')'))
  }

  #normal marginal distributions
  philoglist <- 1:neqt%>%paste0('phi',.,':',logdist)

  #joint distribution, combination of marginals
  jointd  <- dist%>%paste0(., collapse = '*')%>%paste0('joint:',.)
  outobj0 <- c('simp:false ',phi) # do not simplify
  #return marginal and joint distributions
  outobj <- c('print(new)', philoglist , 'print(new)',jointd, 'print(new)',
              paste0('grind(phi',1:neqt,')'), 'print(new)', 'grind(joint)')
  #latex expressions for joint distribution
  outobj <- c(outobj, paste0("tex (''joint)"))

  #evaluate in wxmaxima
  aa <- wxMaxima(outobj0, outobj, tex = FALSE)
  aa
}

#' Help function
#' @keywords internal
ff_generate4maxle_p_v2 <- function(outl, neqt, fixed_term){
  . <- NULL
  phi <- paste0('phi(x, mu, var):= 1/(sqrt(var))*exp((-((x-mu)**2)/(2*var)))')# change to log form
  #phi <- paste0('phi(x, mu, var):= - log(sqrt(var))-((x-mu)**2)/(2*var)')
  logdist <- c()
  dist <- c()
  fixedp <- if(fixed_term==TRUE)'log(1/(sqrt(2*%pi)))+' else ''
  ta1 <- .%>%gsub('\\^', '**', .)
  for(j in 1:neqt){
    obj <- outl[[j]]$expr%>%ta1
    logdist <- c(logdist, paste0(fixedp,'log(phi(eps',j,', ',obj[1],',', obj[2],'))'))
    #logdist <- c(logdist, paste0(fixedp,'(phi(eps',j,', ',obj[1],',', obj[2],'))'))
    dist <- c(dist, paste0('phi(eps',j,', ',obj[1],',', obj[2],')'))
  }
  
  #normal marginal distributions
  philoglist <- 1:neqt%>%paste0('phi',.,':',logdist)
  
  #joint distribution, combination of marginals
  jointd  <- dist%>%paste0(., collapse = '*')%>%paste0('joint:',.)
  outobj0 <- c('simp:false ',phi) # do not simplify
  #return marginal and joint distributions
  outobj <- c('print(new)', philoglist , 'print(new)',jointd, 'print(new)',
              paste0('grind(phi',1:neqt,')'), 'print(new)', 'grind(joint)')
  #latex expressions for joint distribution
  outobj <- c(outobj, paste0("tex (''joint)"))
  
  #evaluate in wxmaxima
  aa <- wxMaxima(outobj0, outobj, tex = FALSE)
  aa
}

#' Log-likelihood expressions for cont. equations
#'
#' @docType data
#'
#' @usage data(expr_ll_norm)
#'
#' @keywords datasets
#'
#' @examples
#' data(expr_ll_norm)
"expr_ll_norm"


#' Another Log-likelihood expressions for cont. equations version 2
#'
#' @docType data
#'
#' @usage data(expr_ll_norm_v2)
#'
#' @keywords datasets
#'
#' @examples
#' data(expr_ll_norm_v2)
"expr_ll_norm_v2"


#' @keywords internal
ff_generate_sigma_expr <- function(neq, word="sigma"){
  . <- NULL
  sigma <- 1:neq%>%expand.grid(., .)%>%apply(., 1, function(x)paste0(x, collapse=','))%>%
    paste0(word, '[', . ,']')%>%sort
  sigma <- matrix(sigma, neq, neq, byrow = TRUE)
  sigma[lower.tri(sigma)] <- sort(sigma[upper.tri(sigma)])
  sigma
}

#' help function
#' @keywords internal
ff_generate4maxle <- function(neq){
  . <- NULL
  sigma <- ff_generate_sigma_expr(neq)

  #create Y
  y <- 1:neq%>%paste0('eps[',.,']')
  ypy <- y%>%paste0('[', ., ']', collapse = ',')%>%
    paste0('Y : matrix(', ., ')')

  #sigma
  spy <- apply(sigma, 1, function(x)paste0("[", paste0("",x,"", collapse = ", "), "]"))%>%paste0(., collapse = ",")
  spy <- paste0('Sigma : matrix(', spy, ')')

  detS <- 'dets : ratsimp(determinant(Sigma))'
  invS <- 'invs : ratsimp(invert(Sigma))'
  expt <- 'expt : ratsimp(transpose(Y) . invs . Y)'

  obj <- c(ypy, spy, detS, invS, expt)
  #grind function in Maxima returns an object that can be mathematically evaluated
  out <- c('print(new)', 'grind(dets)',  'print(new)', 'grind(invs)', 'print(new)', 'grind(expt)')
  rez <- wxMaxima(obj, out, tex = FALSE)
  rez
}

#' help function
#' @keywords internal
ff_generate4mlsem <- function(neq){
  . <- NULL
  sigma <- 1:neq%>%expand.grid(., .)%>%apply(., 1, function(x)paste0(x, collapse=','))%>%
    paste0('sigma', '[', . ,']')%>%sort
  sigma <- matrix(sigma, neq, neq, byrow = TRUE)
  Jt <- sigma%>%gsub("sigma", "deriv", .)
  sigma[lower.tri(sigma)] <- sort(sigma[upper.tri(sigma)])

  #create Y
  y <- 1:neq%>%paste0('eps[',.,']')
  ypy <- y%>%paste0('[', ., ']', collapse = ',')%>%
    paste0('Y : matrix(', ., ')')

  #sigma
  spy <- apply(sigma, 1, function(x)paste0("[", paste0("",x,"", collapse = ", "), "]"))%>%paste0(., collapse = ",")
  spy <- paste0('Sigma : matrix(', spy, ')')

  #Jt
  spj <- apply(Jt, 1, function(x)paste0("[", paste0("",x,"", collapse = ", "), "]"))%>%
    paste0(., collapse = ",")
  spj <- paste0('Jt : matrix(', spj, ')')


  detS <- 'dets : ratsimp(determinant(Sigma))'
  detJ <- 'detJ : ratsimp(determinant(Jt))'
  invS <- 'invs : ratsimp(invert(Sigma))'
  expt <- 'expt : ratsimp(transpose(Y) . invs . Y)'

  obj <- c(ypy, spy, spj, detS, invS, expt, detJ)
  #grind function in Maxima returns an object that can be mathematically evaluated
  out <- c('print(new)', 'grind(dets)',  'print(new)', 'grind(invs)',
           'print(new)', 'grind(expt)', 'print(new)', 'grind(detJ)')
  rez <- wxMaxima(obj, out, tex = FALSE)
  rez
}


#' Log-likelihood expressions for cont. equations
#'
#' @docType data
#'
#' @usage data(datmaxle)
#'
#' @keywords datasets
#'
#' @examples
#' data(datmaxle)
"datmaxle"

#' Log-likelihood expressions for cont. equations sem
#'
#' @docType data
#'
#' @usage data(datmlsem)
#'
#' @keywords datasets
#'
#' @examples
#' data(datmlsem)
"datmlsem"

#' Log-likelihood expressions for cont. equations plus 1 discrete
#'
#' @docType data
#'
#' @usage data(dat4cond_mean_cov_expr)
#'
#' @keywords datasets
#'
#' @examples
#' data(dat4cond_mean_cov_expr)
"dat4cond_mean_cov_expr"


#' Time-use and expenditure dataset 
#' 
#' Data gathered in Austria in 2015 according to Mobility-Activity-Expenditure-Dairy (MAED), which reported all trips, activities (time use) and expenditures of 737 persons over a whole week
#' 
#' Time and expenditure data correspond to weekly totals. Time in hours and expenditure in EUR.
#' 
#' For more on data collection and description see \insertCite{aschauer2018}{nmm} and \insertCite{aschauer2015}{nmm}.
#' 
#' A variant of this dataset was used in:  \insertCite{schmid2017}{nmm},\insertCite{jokubauskaite2018}{nmm} and \insertCite{hoessinger2018}{nmm}.
#' 
#' To get the full dataset please contact r.hoessinger@boku.ac.at.
#' 
#' @docType data
#'
#' @usage data(MAEDtimeExpenditure)
#'
#' @name MAEDtimeExpenditure
#' @format A data frame containing:
#' \describe{
#'   \item{PeID}{individual index}
#'   \item{PeGenF}{gender of the individual}
#'   \item{PeAge}{age in years}
#'   \item{PeEduc}{education level}
#'   \item{PeEmploy}{employment state}
#'   \item{HhCh}{type of household: with children or without children}
#'   \item{w}{hourly wage rate,  EUR/h}
#'   \item{I}{income not realted to work, EUR/week}
#'   \item{Tw}{time spent at work, h/week}
#'   \item{Tf1}{freely chosen activities group 1 (leisure), h/week}
#'   \item{Tf2}{freely chosen activities group 2 (eating, shopping, unspecified), h/week}
#'   \item{Tc}{time spent on committed activities (sleep, domestic work, personal care, travel, education, other), h/week}
#'   \item{Ef1}{freely chosen expenditure group 1 (leisure, accommodation, electronics), EUR/week}
#'   \item{Ef2}{freely chosen expenditure group 2 (clothes), EUR/week}
#'   \item{Ef3}{freely chosen expenditure group 3 (savings), EUR/week}
#'   \item{Ec}{committed expenditures (housing, food, mobility, insurance, other, services, health, furniture, education, financing), EUR/week}
#'   \item{ta}{total time budget = 168 h/week}
#'   \item{Td}{time spent on domestic chores, h/week. Td is part of Tc.}
#' }
#' @importFrom Rdpack reprompt
#' @keywords datasets
#' @references
#' \insertRef{aschauer2015}{nmm}
#' 
#' \insertRef{aschauer2018}{nmm}
#' 
#' \insertRef{hoessinger2018}{nmm}
#' 
#' \insertRef{jokubauskaite2018}{nmm}
#' 
#' \insertRef{schmid2017}{nmm}
#' @examples
#' data(MAEDtimeExpenditure)
"MAEDtimeExpenditure"


#' Trip dataset
#' 
#' Data gathered in Austria in 2015 according to Mobility-Activity-Expenditure-Dairy (MAED), which reported all trips, activities (time use) and expenditures of 737 persons over a whole week.
#' 
#' For more on data collection and description see \insertCite{aschauer2015}{nmm} and \insertCite{aschauer2018}{nmm}.
#' 
#' A variant of this dataset was used in:  \insertCite{schmid2017}{nmm}, \insertCite{jokubauskaite2018}{nmm} and \insertCite{hoessinger2018}{nmm}.
#' 
#' To get the full dataset please contact r.hoessinger@boku.ac.at.
#' 
#' Transport modes available: walk, bike, car, public transport (PT).
#' The inertia variable (int_i) is a dummy, which is equal to one if the mode 
#' chosen by a person for a trip at the start of the current tour is the same as
#' the one chosen in the previous tour made for the same purpose, and zero 
#' otherwise.
#' Variables for trip purpose (leis, work, oth) were created using the effect 
#' coding.
#' @docType data
#'
#' @usage data(MAEDtravel)
#'
#' @name MAEDtravel
#' @format A dataframe containing:
#' \describe{
#'   \item{PeID}{individual index}
#'   \item{PeGenF}{gender of the individual}
#'   \item{PeAge}{age in years}
#'   \item{PeEduc}{education level}
#'   \item{PeEmploy}{employment state}
#'   \item{HhCh}{type of household: with children or without children}
#'   \item{WeID}{trip index}
#'   \item{choice}{chosen mode: 1 - walk, 2 - bike, 3 - car, 4 -public transport (PT)}
#'   \item{dist}{trip distance, km}
#'   \item{avl_1}{availability dummy for mode 1, walk}
#'   \item{avl_2}{availability dummy for mode 2, bike}
#'   \item{avl_3}{availability dummy for mode 3, car}
#'   \item{avl_4}{availability dummy for mode 4, PT}
#'   \item{chc_1}{choice dummy for mode 1, walk}
#'   \item{chc_2}{choice dummy for mode 2, bike}
#'   \item{chc_3}{choice dummy for mode 3, car}
#'   \item{chc_4}{choice dummy for mode 4, PT}
#'   \item{cost_3}{cost of car mode}
#'   \item{cost_4}{cost of PT mode}
#'   \item{dur_1}{trip duration with mode 1, minutes}
#'   \item{dur_2}{trip duration with mode 2, minutes}
#'   \item{dur_3}{trip duration with mode 3, minutes}
#'   \item{vdur_4}{in vehicle time mode 4, minutes}
#'   \item{acc_4}{time to stop or from stop for mode 4, minutes}
#'   \item{HhCarPark}{dummy, car parking at home available}
#'   \item{JobCarPark}{dummy,car parking at workplace available}
#'   \item{PbAvl_3}{dummy, ar parking restrictions (time and/or cost) in-force at the destination of the trip}
#'   \item{servIdx_4}{public transport service interval in minutes}
#'   \item{stopUs1R1_4}{necessary number of changes to reach the destination with public transport}
#'   \item{leis}{trip purpose leisure, effect coding}
#'   \item{work}{trip purpose work, effect coding}
#'   \item{oth}{trip purpose other, effect coding}
#'   \item{int_1}{inertia for mode 1}
#'   \item{int_2}{inertia for mode 2}
#'   \item{int_3}{inertia for mode 3}
#'   \item{int_4}{inertia for mode 4}}
#' @keywords datasets
#' @references
#' \insertRef{aschauer2018}{nmm}
#' 
#' \insertRef{aschauer2015}{nmm}
#' 
#' \insertRef{hoessinger2018}{nmm}
#' 
#' \insertRef{jokubauskaite2018}{nmm}
#' 
#' \insertRef{schmid2017}{nmm}
#' @examples
#' data(MAEDtravel)
"MAEDtravel"




#' Example dataset 
#' 
#' Data "MathPlacement" taken from Stat2Data package.  
#' 
#' Code for data modifications can be found in the example section.
#' 
#' @docType data
#'
#' @usage data(dataM)
#'
#' @name dataM
#' @format A data frame containing:
#' \describe{
#'   \item{Student}{Identification number for each student}
#'   \item{Gender}{0=Female, 1=Male}
#'   \item{PSATM}{PSAT score in Math}
#'   \item{SATM}{SAT score in Math}
#'   \item{ACTM}{ACT Score in Math}
#'   \item{Rank}{Adjusted rank in HS class}
#'   \item{Size}{Number of students in HS class}
#'   \item{GPAadj}{Adjusted GPA}
#'   \item{PlcmtScore}{Score on math placement exam}
#'   \item{Recommends}{Recommended course: R0 R01 R1 R12 R2 R3 R4 R6 R8}
#'   \item{Course}{Actual course taken}
#'   \item{Grade}{Course grade}
#'   \item{RecTaken}{1=recommended course, 0=otherwise}
#'   \item{TooHigh}{1=took course above recommended, 0=otherwise}
#'   \item{TooLow}{1=took course below recommended, 0=otherwise}
#'   \item{CourseSuccess}{1=B or better grade, 0=grade below B}
#'   \item{DR_Course}{according to recommendations, which level of course was taken: alow - lower, bnormal - recommended, chigh - higher}
#' }
#' @importFrom Rdpack reprompt
#' @keywords datasets
#' @examples
#' \donttest{
#' data(dataM)
#' library(magrittr)
#' library(dplyr)
#' if (requireNamespace("recipes", quietly = TRUE)&requireNamespace("Stat2Data", quietly = TRUE)) {
#' data("MathPlacement", package="Stat2Data")
#' head(MathPlacement)
#' library(recipes)
#' # As some of the data is missing, k-nearest neighbors (knn) imputation is 
#' # used to fill the gaps. This is done with recipes package and function 
#' # step_knnimpute.
#' dataM <- recipe(~ ., data = MathPlacement) %>%
#' step_knnimpute(everything()) %>% prep() %>% juice()
#' # Afterwards we create a categorical variable that will show whether a 
#' # student took a course which was too high, too low, the recommended one or
#' # something else happened:
#' dataM %<>% mutate(Student = 1:n(), DR_Course = case_when(
#' TooHigh == 1 ~ "chigh",
#' TooLow == 1 ~ "alow",
#' RecTaken == 1 ~ "bnormal",
#' TRUE ~"dother"
#' )) 
#' # We remove observations with ambiguous course status:
#' dataM %<>% filter(DR_Course!="dother")
#' dataM %>% select(DR_Course) %>% table %>% t 
#' }
#' }
"dataM"
