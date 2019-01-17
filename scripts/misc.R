# return sorted unique values ---------------------------------------------
id <- function(x) unique(sort(as.vector(unlist(x))))


# standardize a vector ----------------------------------------------------
standardize <- function(x)
  (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)


# expit and logit functions -----------------------------------------------
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))


# load and return loaded object -------------------------------------------
load.local <- function(file) {
 v <- load(file)
 stopifnot(length(v) == 1)
 get(v)
}


# nice little pdf function ------------------------------------------------
pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}


# Extraction of the mode of a posterior distribution  ---------------------
################################################
# multiple input (x) possible to this function 
## 1. With the simulations matrix (only one of the chains will be represented here)
# model.summary$BUGSoutput$sims.matrix
## 2. With the parameters extracted from the array 
# param=model.summary$parameters.to.save
# param
# model.summary$BUGSoutput$sims.array[,,param[3]]
## 3. By calling the parameter in the list 
# model.summary$BUGSoutput$sims.list$alpha.0
################################################

post.mode=function (x, adjust = 0.1, ...) {  
  find.mode <- function(x, adjust, ...) {
    dx <- density(x, adjust = adjust, ...)
    dx$x[which.max(dx$y)]
  }
  apply(as.matrix(x), 2, find.mode, adjust = adjust, ...)
}




# Pairs plots -------------------------------------------------------------

# To create panels in pairs function  -------------------------------------

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gold", ...) 
}

# Pairs plot  -------------------------------------------------------------
autopairs<-function(x, ...) {
  pairs(x, 
        lower.panel = panel.smooth, 
        upper.panel = panel.cor, 
        diag.panel = panel.hist, ...)
}


# Stats: mode -------------------------------------------------------------
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# Compute Galton's skewness measure for x ---------------------------------
#  NOTE: this procedure assumes no x values are missing
#

galtonskew.proc <- function(x){
  quarts <- as.numeric(quantile(x, probs = c(0.25, 0.5, 0.75)))
  num <- quarts[1] + quarts[3] - 2*quarts[2]
  denom <- quarts[3] - quarts[1]
  gskew <- num/denom
  gskew
}


# invlogit and logit ------------------------------------------------------
# Function that will transforme the Y values from linear scale to [0,1]
invlogit <- function(x){exp(x)/(exp(x) + 1)}
logi <-     function(x){exp(x)/(1+exp(x))}


# SE ----------------------------------------------------------------------
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))


