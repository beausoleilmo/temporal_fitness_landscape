### ### ### ### ### ### ### ### ### ### 
# R function to do the calculations in Janzen and Stern, Evolution, 1998
# Janzen, F. J., and H. S. Stern. 1998. Logistic Regression for Empirical Studies of Multivariate Selection. Evolution 52:1564–1571.
### ### ### ### ### ### ### ### ### ### 

### ### ### ### ### ### ### ### ### ### 
# Contact 
# Fred Janzen: fjanzen@iastate.edu 
# Hal Stern: sternh@uci.edu

# Adapted to R by Marc-Olivier Beausoleil July 5 2018 
# Commented and explaination added 
# marc-olivier.beausoleil@mail.mcgill.ca 
### ### ### ### ### ### ### ### ### ### 

### ### ### ### ### ### ### ### ### ### 
# Janzen and Stern published a paper in Evolution (1998, 52:1564-1571) describing the use of logistic regression for analyzing selection and a new method for transforming the logistic regression coefficients into selection gradients that can be used in equations describing microevolutionary change. It should be noted that the former is not a novel contribution to the field (others have also used logistic regression for this purpose) but the latter is. 
# The point of this web page (http://www.public.iastate.edu/~fjanzen/software/regress.htm) is simply to make it easier for users to implement these statistical methods in their own work on selection. We provide a SAS program that performs the logistic regression analysis and calculates the transformed selection gradients. We also provide the dataset from Bumpus' classic study of sparrow mortality and the output from the SAS logistic regression analysis of a subset (the males) of those data.

# One important point in performing logistic regression in SAS is to be sure that your model predicts the probability of survival rather than the probability of mortality. The program below achieves this goal by having the data sorted with survivors first and using the order=data option in PROC LOGISTIC. Another option, assuming the traditional 0=dead/1=survived coding, is to use the DESCENDING option. 

# You will likely need to modify the SAS program to analyze your own datasets. It is especially important to note where the survival variable and the predictors appear in the input statement as this arrangement affects some later parts of the program. 
# Please contact Fred Janzen email: fjanzen@iastate.edu if you have thoughts, comments, or suggestions. Note that an S-PLUS program is also available (and considerably simpler). If interested, please contact Hal Stern email: sternh@uci.edu 
### ### ### ### ### ### ### ### ### ### 

### ### ### ### ### ### ### ### ### ### 
# Description of function 

# x is a matrix of predictors
# y is a vector of binary outcomes
# produces a table with:
#    first 3 columns = linear regression slopes, std err, p-values
#    next 3 columns = logistic regression slopes, std err, p-values
#    final column = average gradient transformation of logistic regression slope
#
# Note: Final 3 lines run the analysis on the male Bumpus data 
#       (without the log transform used in the paper)
#       That data is provided below along with the output
# 
# For more information go here: http://www.public.iastate.edu/~fjanzen/software/regress.htm 
### ### ### ### ### ### ### ### ### ### 

lgtsel <- function(x, # trait values
                   y  # Binomial response 
                   )  {
  x = as.data.frame(x)
    if (length(dim(x)) != 2) dim(x) <- c(length(x),1)
    
    # For the linear comparison 
    alin <- lsfit(x,y) # Find the Least Squares Fit
    blin <- ls.diag(alin) # Compute Diagnostics for lsfit Regression Results
    # For the logistic  comparison 
    algt <- glm(y ~ ., family = binomial(link = "logit"), data=x) # logistic regression
    blgt <- summary(algt)
    
    n <- nrow(x)
    np <- ncol(x)
    rw <- n/sum(y) # Relative fitness 
    mn <- apply(x,2,mean) # mean of each trait 
    md <- apply(x,2,median) # median of each trait 
    sd <- sqrt(apply(x,2,var)) # standard deviation, same as apply(x,2,sd) 
    
    # For the linear comparison 
    c1 <- alin$coef[2:(np+1)]*sd*rw # this is rescaling the coefficients with sd of the traits AND relative fitness; (2:(np+1) is removing the intercept)
    c2 <- blin$std.err[2:(np+1)]*sd*rw # this is rescaling the std.err with sd of the traits AND w 
    c3 <- 2*(1-pt(q = abs(c1)/c2,
                  df = n-np-1)) # pt: The Student t Distribution; q is vector of quantiles; n-np-1 are the degrees of freedom; 2* is because it's 2 tailed 
    # For the logistic  comparison 
    c4 <- blgt$coefficients[2:(np+1),1]*sd # Estimates of the GLM (logistic regression), scaled by the sd 
    c5 <- blgt$coefficients[2:(np+1),2]*sd # SE of the GLM (logistic regression), scaled by the sd 
    c6 <- 2*(1-pnorm(q = abs(c4)/c5, 
                     mean = 0, 
                     sd = 1, 
                     lower.tail = TRUE, 
                     log.p = FALSE)) # pnorm: Normal Distribution to find the 2 tailed P-value (2* and abs() function)
    
    # Get the predicted values on the LINEAR SCALE; verify like this range(z)
    z <- blgt$coefficients[1,1] + # intercept 
      as.matrix(x) %*% blgt$coefficients[2:(np+1),1] # α_0 + α_1*x_1 + α_i x_i + ... + α_n + x_n
    # Transform 
    invlogit <- function(x) {exp(x)/(1+exp(x))}
    pr <- invlogit(z) # transormation to the logistic scale [0,1]; verify like this range(pr)
    c7 <- mean(pr*(1-pr))*c4*rw # pr*(1-pr) refers to {\hat{p}*(1-\hat{p})}, c4 contains the estimates from the GLM and rw is the relative fitness. This line is an approximation of the average selection gradient
    c8 <- mean(pr*(1-pr))*c5*rw 
    ag = mean(pr*(1-pr))
    cag= rep(ag, length(c1))
    Relfit = rep(1/rw, length(c1))
    c0 = blgt$coefficients[2:(np+1),1]
    # Create the table 
    tab <- cbind(c0, c1,c2,c3,c4,c5,c6,c7,c8, cag, Relfit,sd)
    dimnames(tab)[[2]] <- c("coeff", 
                            "lin.beta","lin.se","lin.p",
                            "lgt.alph","lgt.se","lgt.p",
                            "bavggrd","bavggrds", "ag",
                            "relfit", "sd")
    # round(
      as.data.frame(tab)
      # ,3)
}

### ### ### ### ### ### ### ### ### ### 
# ## Example 
# bp = read.table("~/Desktop/Bumpus Data.txt", header = TRUE)
## filter data
# sex = 0 # males, 0 are females
# bumpusm = bp[bp$sex== sex ,]
# y <- bumpusm[,2]
# x <- bumpusm[,3:11]
# t1 <- lgtsel(x,y)
# t1
## Look at the difference in terms of gradients between linear and logistic regression
# abs(t1$lin.beta) - abs(t1$bavggrd)
### ### ### ### ### ### ### ### ### ### 
# library(vegan)
# library(rgl)
# library(ggplot2)
# dput(names(bumpusm))
# source("~/Dropbox/finch_recap/v2/data/creating_data/PCA_custom.R")
# pca.out=rda(bumpusm[,c("total.length", "alar.extent", "weight", 
#                "length.beak.head", "length.humerus", "length.femur", "length.tibiotarsus", 
#                "width.skull", "length.keel.sternum")],scale = TRUE)
# custom_pca(pca.out, centered = TRUE)
# x.pca = scores(pca.out,choices = c(1,2))$sites
                          ## PC1         PC2
### total.length        1.175058 -0.59223016 # 
### alar.extent         1.269677 -0.13762706 # 
### weight              1.155034 -0.52266424 # 
### length.beak.head    1.329730  0.11473985 # 
### length.humerus      1.360582  0.25750859 # -
### length.femur        1.262479  0.67516842 # -
### length.tibiotarsus  1.269448  0.64851011 # -
### width.skull         1.183673 -0.02985344 # 
### length.keel.sternum 1.088915 -0.61427424 # +
# lm.out.pc=lm(y~x.pca+I(x.pca^2))
# summary(lm.out.pc)
# coef(lm.out.pc)[4:5]*2
# x.df = cbind(x.pca,x.pca^2)
# colnames(x.df) <- c("pc1","pc2","pc1.2","pc2.2")
# row.names(x.df) <- 1:nrow(x.df)
# lgtsel(x.df[,c(1,3)],y)
# pred <- predict(lm.out.pc)
# plot(y~x.pca[,1])
# lines(x.pca[,1], y=pred, col = "blue")
# df = as.data.frame(cbind(x.df,y))
# lm.out.pc=lm(y~pc1+pc1.2, data = df)
# newdata <- data.frame(pc1  =seq(min(df$pc1), max(df$pc1),     length.out = 1000),
#                       pc2  =seq(min(df$pc2), max(df$pc2),     length.out = 1000),
#                       pc1.2=seq(min(df$pc1.2), max(df$pc1.2), length.out = 1000),
#                       pc2.2=seq(min(df$pc2.2), max(df$pc2.2), length.out = 1000))
# newdata$pred1 <- predict(lm.out.pc, newdata)
# plot(y ~ pc2, data = df)
# lines(newdata$pc2, newdata$pred1, col = "red")
# ggplot(df, aes(y=y, x=c(pc1))) + 
#   geom_point(alpha = .5) + 
#   stat_smooth(method = "lm", formula = y ~ poly(x,2))

# plot3d(x = x.pca[,1],y = x.pca[,2],z = y)
