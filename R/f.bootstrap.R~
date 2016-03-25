#' A bootstrapper that utilizes the same cpp function to retain xts rownames via a hashing procedure
#' Takes a single column input with date rownames , generates a bootstrap sample and returns it
#' Extensible to multiple column input ( via either a separate wrapper, or h.boot.sample )
#' The e parameter guarantees that cases of an integer will map to the right power, assuming no significant numerical error
#' @export

xtshashboot<-function(x,b,type,e=0.00001)
{
    if (type == "stationary") {type = 0}
    else if (type == "circular") {type = 1}
    else {stop("only stationary and circular bootstrap implemented")}

    x = as.matrix(x)

    date<-rownames(x)
    mindate<-min(date)
    date<-as.Date(date)-as.Date(mindate)
    maxelement=max(x)-min(x)
    minelement=min(x)
    x<-x-minelement
    power = ceiling(log10(maxelement+e))
    x<-x+date*(10^(power))
    bootstrap<-f_bootstrap(x,b,type)
    bootstrapdata<-bootstrap%%(10^(power))
    bootstrapdates<-as.Date(mindate)+((bootstrap-bootstrapdata)/(10^(power)))
    bootstrap<-xts(bootstrapdata+minelement,bootstrapdates)
    return (bootstrap) 
} 

#' Wrapper cum Iterator that applies the bootstrap with a generic function 
#' For xts objects only 
#' NON-FUNCTIONAL FOR SHARPE ! As these PerfA metrics do not work with bootstrap samples that repeat elements ( try Sharpe to observe this )
#' Hence, running the function on any sample with repetitions bugs it out
#' Demonstration purposes for now mostly
#' @export

xtsbootfunc<-function (x,b,type,fun,rep=100,e=0.00001,...)
{
    x<-as.matrix(x)
    var<-c()

    for (j in 1:ncol(x))

    { 
    hold<-c()
    for (i in 1:rep) { hold<-c(hold,fun(nse:::xtshashboot(x[,j],b,type,e))) }
    var<-c(var,var(hold))  
    }

    var<-as.matrix(var)
    colnames(var)<-paste("Bootstrapped Variance")
    rownames(var)<-colnames(x) 
    return (t(var))

}

h.boot.sample <-
  function(x, b, type)
  {
    if (type == "stationary") {type = 0}
    else if (type == "circular") {type = 1}
    else {stop("only stationary and circular bootstrap implemented")}
    return(f_bootstrap(x, b, type))
  }

.f.bootstrap <- function(x, nb = 1, statistic = NULL, b = NULL, type, ...)
{
  y <- embed(x, 1)
  if(is.null(statistic)) {
    n = NROW(y)
    boot <- matrix(y, nrow=n, ncol=nb)
    out <- apply(boot, 2, h.boot.sample, b, type)
    return(drop(out))
  }
  else {
    
    yi <- 1:NROW(y)
    orig.statistic <- statistic(y,...)
    l.stat <- length(orig.statistic)
    stat <- matrix(0, nb, l.stat)
    for(i in 1:nb){
      stat[i,] <- statistic(as.matrix(y[h.boot.sample(yi, b, type),]),...)
    }
    out <- list(statistic = stat)
    
    return(out)
  }
}

# Bootstrap samples generation
# @description Generate bootstrap samples
# @details Two bootstrap schemes are available, the stationary bootstrap of Politis and Romano  (1994)
# and the circular bootstrap of Politis and Romano (1992)
#     @param x       A numeric vector or a matrix
#     @param nb   The number of bootstrap replication
#     @param statistic Function to be used on the bootstrap replication
#     @param b  The block length
#     @param type    The bootstrap schemes c("stationary","circular")
#     @param seed    The seed for the bootstrap simulation
#     @param ...    Aditional argument passed to the statistic function
#     @return  The bootstrap samples or the statistic if statistic is non-null
#    @references Politis, Dimitris N., and Joseph P. Romano. "A circular block-resampling procedure for stationary data." Exploring the limits of bootstrap (1992): 263-270.
#    @references Politis, Dimitris N., and Halbert White. "Automatic block-length selection for the dependent bootstrap." Econometric Reviews 23.1 (2004): 53-70.
#    @references Politis, Dimitris N., and Joseph P. Romano. "The stationary bootstrap." Journal of the American Statistical association 89.428 (1994): 1303-1313.
#@import np
f.bootstrap = compiler::cmpfun(.f.bootstrap)
