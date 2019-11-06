
load("./output/logistic_gam_models.selection.RData", verbose = TRUE)
# Model GLM with all data -------------------------------------------------
load("./output/climate.year.month.database.RData", 
     verbose = TRUE) # created from ~/Dropbox/finch_recap/v2/data/creating_data/Climate_all_time_Santa_cruz.R
db.climate.year # Per year precipitation, sum, mean, var, sd, se 

as.factor(df.final.df.new.analysis$y)

# This is to find the rows in each year to get the rain data  
row.precip = match(c(as.numeric(as.character(df.final.df.new.analysis$yea.var)) -2),db.climate.year$Year.marco)
# I decided to look for the rain data in the complete year (wet and dry combined) 
df.final.df.new.analysis$precipitation.yr.before = db.climate.year[row.precip,"sum.prec"]
# I decided to look for the rain data in the dry season only 
df.final.df.new.analysis$prcp.yr.before.dry = dry.rain.our.study[row.precip,"sum.prec"]
df.final.df.new.analysis$prcp.yr.before.wet = wet.rain.our.study[row.precip,"sum.prec"]

df.final.df.new.analysis = df.final.df.new.analysis[df.final.df.new.analysis$yea.var != 2007,]
plot(y~my.x, 
     data = df.final.df.new.analysis)
df.final.df.new.analysis$x = df.final.df.new.analysis$my.x
df.final.df.new.analysis$x2 = df.final.df.new.analysis$my.x^2
df.final.df.new.analysis$yea.var  = as.factor(df.final.df.new.analysis$yea.var)
glm.out1 =  glm(y~ x,                  data = df.final.df.new.analysis,family = "binomial")
glm.out2 =  glm(y~ x + x2,             data = df.final.df.new.analysis,family = "binomial")
glm.out3 =  glm(y~ x + x2 + yea.var,   data = df.final.df.new.analysis,family = "binomial")
glm.out4 =  glm(y~ x + x2 * yea.var,   data = df.final.df.new.analysis,family = "binomial")
glm.out5 =  glm(y~ x + x2 : yea.var,   data = df.final.df.new.analysis,family = "binomial")
glm.out6 =  glm(y~(x + x2) : yea.var,  data = df.final.df.new.analysis,family = "binomial")
glm.out7 =  glm(y~ x * yea.var + x2 * yea.var, data = df.final.df.new.analysis,family = "binomial")
glm.out8 =  glm(y~ x:x2:yea.var, data = df.final.df.new.analysis,family = "binomial")
glm.out9 = glm(y~ x:yea.var+ x2:yea.var, data = df.final.df.new.analysis,family = "binomial")
glm.out10 = glm(y~ x:yea.var+ x2:yea.var + yea.var, data = df.final.df.new.analysis,family = "binomial")
glm.out11 =  glm(y~(x + x2) * yea.var,  data = df.final.df.new.analysis,family = "binomial")
glm.out12 =  glm(y~ x*x2*yea.var, data = df.final.df.new.analysis,family = "binomial")
glm.out13 =  glm(y~ x*x2*yea.var*precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out14 =  glm(y~ x + x2 * precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out15 =  glm(y~ x + x2 : precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out16 =  glm(y~ (x + x2)* precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
# glm.out3 =  glm(y~ x + x2,             data = df.final.df.new.analysis,family = "binomial")
summary(glm.out11)
summary(glm.out16)

glm.out17 =  glm(y~ x + x2 + precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out18 =  glm(y~ x + x2 * precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out19 =  glm(y~ x + x2 : precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out20 =  glm(y~(x + x2) : precipitation.yr.before,  data = df.final.df.new.analysis,family = "binomial")
glm.out21 =  glm(y~ x * precipitation.yr.before + x2 * precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out22 =  glm(y~ x:x2:precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out23 = glm(y~ x:precipitation.yr.before+ x2:precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out24 = glm(y~ x:precipitation.yr.before+ x2:precipitation.yr.before + precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out25 =  glm(y~(x + x2) * precipitation.yr.before,  data = df.final.df.new.analysis,family = "binomial")
glm.out26 =  glm(y~ x*x2*precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out27 =  glm(y~ x*x2*precipitation.yr.before*precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out28 =  glm(y~ x + x2 * precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out29 =  glm(y~ x + x2 : precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
table(df.final.df.new.analysis$yea.var)
glm.out30 =  glm(y~ (x + x2)* precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")

glm.out31 =  glm(y~ (x + x2)* prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")
glm.out32 =  glm(y~ x + (x2)* prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")
glm.out33 =  glm(y~ x + x2 + prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")
glm.out34 =  glm(y~ x + x2: prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")


glm.out28 =  glm(y ~  x + x2  * precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out32 =  glm(y ~  x + (x2)* prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")
glm.out35.1 =  glm(y~ x + (x2)* prcp.yr.before.wet,   data = df.final.df.new.analysis,family = "binomial")
summary(glm.out28)
summary(glm.out32)
summary(glm.out35.1)


summary(glm.out31)
summary(glm.out32)
summary(glm.out33)
summary(glm.out34)
summary(glm.out35.1)


df.final.df.new.analysis$prc.bfr = df.final.df.new.analysis$precipitation.yr.before
glm.out102=summary(glm(y~x2 * prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out103=summary(glm(y~x + x2 * prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out104=summary(glm(y~(x + x2) * prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out105=summary(glm(y~ x2+ prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out106=summary(glm(y~ x + x2,data=df.final.df.new.analysis,family="binomial"))
glm.out107=summary(glm(y~ x2,data=df.final.df.new.analysis,family="binomial"))
glm.out108=summary(glm(y~ x+x2+prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out109=summary(glm(y~ x*prc.bfr+x2,data=df.final.df.new.analysis,family="binomial"))
glm.out110=summary(glm(y~ prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out111=summary(glm(y~ 1,data=df.final.df.new.analysis,family="binomial"))
glm.out112=summary(glm(y~ x*prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out113=summary(glm(y~ x+prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out114=summary(glm(y~ x,data=df.final.df.new.analysis,family="binomial"))

mod.l=list(glm.out102,glm.out103,glm.out104,
           glm.out105,glm.out106,glm.out107,glm.out108,
           glm.out109,glm.out110,glm.out111,
           glm.out112,glm.out113,glm.out114)
numb.modl = unlist(lapply(lapply(mod.l,coefficients),nrow))

test.all.var = do.call("rbind", (lapply(mod.l,coefficients)))
rown= row.names(test.all.var)

test.all.var = data.frame(test.all.var)
test.all.var$var= rown
mlistt = NULL
for (i in 1:length(numb.modl)) {
  mlistt = c(mlistt,rep(i,numb.modl[i]))
}
test.all.var$model = (mlistt)
write.csv(test.all.var,"./output/glm.model.p.values.coefficients.csv",na = "")

# model with just interaction 
# model GLMM precipitation (fixed), year (random)
summary(glm.out1)
summary(glm.out2)
summary(glm.out3)
summary(glm.out4) 
summary(glm.out5) 
summary(glm.out6) 
summary(glm.out7) 
summary(glm.out8) 
summary(glm.out9) 
summary(glm.out10) 
summary(glm.out11) 
summary(glm.out12) 
summary(glm.out13) 
summary(glm.out14) 
summary(glm.out15) 
summary(glm.out16) 

summary(glm.out17)
summary(glm.out18)
# summary(glm.out19)
# summary(glm.out20) 
summary(glm.out21) 
# summary(glm.out22) 
# summary(glm.out23)
# summary(glm.out24) 
summary(glm.out25) 
# summary(glm.out26) 
# summary(glm.out27) 
summary(glm.out28) 
# summary(glm.out29) 
summary(glm.out30) 

old.options=options()
options(na.action = na.fail)
# see https://rcompanion.org/rcompanion/d_04.html 
options(contrasts = c("contr.sum", "contr.poly"))
### needed for type III tests

glm.out28 =  glm(y ~  (x + x2) * precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out32 =  glm(y ~  (x + x2) * prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")
glm.out35.1 =  glm(y~ (x + x2) * prcp.yr.before.wet,   data = df.final.df.new.analysis,family = "binomial")
model.sel28 = dredge(glm.out28)
model.sel32 = dredge(glm.out32)
model.sel35 = dredge(glm.out35.1)
# Automated model selection
model.sel0 = dredge(glm.out7)
model.sel1 = dredge(glm.out30)
model.sel2 = dredge(glm.out31)
model.sel2 = dredge(glm.out31)
# model.sel2 = dredge(glm.out13) #too much to look at! 
model.avg(model.sel0,subset = delta < 4)
model.avg(model.sel,subset = delta < 4)
att.mod0 =attr(model.sel0,"model.calls")
att.mod1 =attr(model.sel1,"model.calls")
att.mod2 =attr(model.sel2,"model.calls")
mod.0 = gsub(",     ", replacement = "", paste(format(att.mod0)))
mod.0 = gsub('\"', replacement = "", mod.0)
mod.1 = gsub(",     ", replacement = "", paste(format(att.mod1)))
mod.1 = gsub('\"', replacement = "", mod.1)
mod.2 = gsub(",     ", replacement = "", paste(format(att.mod2)))
mod.2 = gsub('\"', replacement = "", mod.2)

model.sel0$models = mod.0
model.sel1$models = mod.1
model.sel2$models = mod.2

write.csv(model.sel0,"./output/glm.all.data_year.csv",na = "")
write.csv(model.sel1,"./output/glm.all.data_precip.csv",na = "")

write.csv(model.sel2,"./output/glm.all.data_year_precip.csv")
model.sel1
model.sel2
str(model.sel1)
coeffs(model.sel1)
attr(model.sel1,"coefTable")
row.names(model.sel0[model.sel0$delta <=4,])
attr(model.sel0,"coefTable")
attr(model.sel0,"model.calls")$`8`
attr(model.sel0,"model.calls")$`4`
attr(model.sel0,"model.calls")$`3`
attr(model.sel0,"model.calls")$`16`
attr(model.sel0,"model.calls")$`24`

attr(model.sel1,"model.calls")$`22`
summary(glm(formula = y ~ x + x2 + yea.var + 1, family = "binomial", 
            data = df.final.df.new.analysis))
summary(glm(formula = y ~ x + x2 + 1, family = "binomial", data = df.final.df.new.analysis))
summary(glm(formula = y ~ x2 + 1, family = "binomial", data = df.final.df.new.analysis))
summary(glm(formula = y ~ x + x2 + yea.var + x:yea.var + 1, family = "binomial", 
            data = df.final.df.new.analysis))
summary(glm(formula = y ~ x + x2 + yea.var + x2:yea.var + 1, family = "binomial", 
            data = df.final.df.new.analysis))

Anova(glm.out1,type = "III", test.statistic = "Wald")
Anova(glm.out2,type = "III", test.statistic = "Wald")
Anova(glm.out3,type = "III", test.statistic = "Wald")
Anova(glm.out4,type = "III", test.statistic = "Wald")
Anova(glm.out5,type = "III", test.statistic = "Wald")
Anova(glm.out6,type = "III", test.statistic = "Wald")
Anova(glm.out7,type = "III", test.statistic = "Wald")
Anova(glm.out8,type = "III", test.statistic = "Wald")
Anova(glm.out9,type = "III", test.statistic = "Wald")
Anova(glm.out10,type = "II", test.statistic = "Wald")
Anova(glm.out10,type = "III", test.statistic = "Wald")
Anova(glm.out11,type = "III", test.statistic = "Wald")
Anova(glm.out12,type = "III", test.statistic = "Wald")
Anova(glm.out13,type = "III", test.statistic = "Wald")
Anova(glm.out14,type = "III", test.statistic = "Wald")
Anova(glm.out15,type = "III", test.statistic = "Wald")
Anova(glm.out16,type = "III", test.statistic = "Wald")

Anova(glm.out17,type = "III", test.statistic = "Wald")
Anova(glm.out18,type = "III", test.statistic = "Wald")
Anova(glm.out19,type = "III", test.statistic = "Wald")
Anova(glm.out20,type = "III", test.statistic = "Wald")
Anova(glm.out21,type = "III", test.statistic = "Wald")
Anova(glm.out22,type = "III", test.statistic = "Wald")
Anova(glm.out23,type = "III", test.statistic = "Wald")
Anova(glm.out24,type = "III", test.statistic = "Wald")
Anova(glm.out25,type = "III", test.statistic = "Wald") #***
Anova(glm.out26,type = "II", test.statistic = "Wald")
Anova(glm.out26,type = "III", test.statistic = "Wald")
Anova(glm.out27,type = "III", test.statistic = "Wald")
Anova(glm.out28,type = "III", test.statistic = "Wald")
Anova(glm.out29,type = "III", test.statistic = "Wald")
Anova(glm.out30,type = "III", test.statistic = "Wald")

attr(model.sel1,"model.calls")$`32`

modmod =glm(formula = y ~ precipitation.yr.before + x2 + precipitation.yr.before:x2 + 
              1, family = "binomial", data = df.final.df.new.analysis)
Anova(modmod, type = "III", test.statistic = "Wald")

modmod =glm(formula = y ~ precipitation.yr.before + x + x2 + precipitation.yr.before:x2 + 
              1, family = "binomial", data = df.final.df.new.analysis)
Anova(modmod, type = "III", test.statistic = "Wald")

modmod =glm(formula = y ~ precipitation.yr.before + x + x2 + precipitation.yr.before:x + 
              precipitation.yr.before:x2 + 1, family = "binomial", data = df.final.df.new.analysis)
Anova(modmod, type = "III", test.statistic = "Wald")

Anova(glm.out1,type = "II", test.statistic = "Wald")
Anova(glm.out2,type = "II", test.statistic = "Wald")
Anova(glm.out3,type = "II", test.statistic = "Wald")
Anova(glm.out4,type = "II", test.statistic = "Wald")
Anova(glm.out5,type = "II", test.statistic = "Wald")
Anova(glm.out6,type = "II", test.statistic = "Wald")
Anova(glm.out7,type = "II", test.statistic = "Wald")
Anova(glm.out8,type = "II", test.statistic = "Wald")
Anova(glm.out9,type = "II", test.statistic = "Wald")
Anova(glm.out10,type = "II", test.statistic = "Wald")
Anova(glm.out11,type = "II", test.statistic = "Wald")
Anova(glm.out12,type = "II", test.statistic = "Wald")

anova(glm.out1,glm.out2,
      glm.out3,glm.out4,glm.out5,glm.out5,glm.out6,glm.out7,glm.out8,glm.out9,glm.out10)

df.final.df.new.analysis$pcr.scale = scale(df.final.df.new.analysis$prc.bfr)
rd.modl = glmer(formula = y ~ prc.bfr + 
                  x + x2 + 
                  prc.bfr:x + prc.bfr:x2 + 
                  (1 | my.band),
                family = "binomial", data = df.final.df.new.analysis)
rd.modl1 = glmer(formula = y ~ x + x2 + prc.bfr+
                   (1 | my.band), family = "binomial", 
                 data = df.final.df.new.analysis)
rd.modl1 = glmer(formula = y ~ x + x2* prc.bfr+
                   (1 | my.band), family = "binomial", 
                 data = df.final.df.new.analysis)
summary(rd.modl)
summary(rd.modl1)
dredge(rd.modl) # Takes a long time 
sjp.glmer(rd.modl, type = "fe")

par(mfrow=c(1,1))

# plot(glm.out)
newx <- seq(from = min(df.final.df.new.analysis$my.x), 
            to = max(df.final.df.new.analysis$my.x), 
            length.out = 2000)

z1 <- predict(glm.out1, 
              newdata = data.frame(x = newx,  
                                   x2 = newx^2, 
                                   yea.var = 2009), 
              se.fit = TRUE)
# Function that will transforme the Y values from linear scale to [0,1]
# This is the actual transformed data 
yhat <-  invlogit(z1$fit)
upper <- invlogit(z1$fit + z1$se.fit)
lower <- invlogit(z1$fit - z1$se.fit)

# satisfied...= "n"
# while(satisfied... != "y"){
# Here is the plot of the fitness function using the evenly spaced data (newx) and the response to it (yhat)
plot(newx, yhat, type="l", 
     ylim = c(0,1), 
     xlab = "PC1",
     main = "GLM all data")
# Adding error 
lines(newx, upper, lty = 2)
lines(newx, lower, lty = 2)
points(df.final.df.new.analysis$my.x,df.final.df.new.analysis$y)
# par(new = TRUE)



# LDA ---------------------------------------------------------------------
mld = lda(y~my.mbl +my.mbw + my.mbd,data = df.final.df.new.analysis)
mqd = qda(y~my.mbl +my.mbw + my.mbd,data = df.final.df.new.analysis)
mld
mqd
plot(mld)
# par(mfrow=c(4,2))
# for (j  in 1:length(unique(df.final.df.new.analysis$yea.var))) {
#   all.y = unique(df.final.df.new.analysis$yea.var)
#   my.df.for = df.final.df.new.analysis[df.final.df.new.analysis$yea.var == all.y[j],]
# partimat(as.factor(my.df.for$y) ~ my.mbl +my.mbw + my.mbd,
#          data = my.df.for, method = "lda", plot.matrix = TRUE)
# }
if (pdf.output) {
  pdf("./output/quadratic_discriminant_function.pdf")
}
# partimat(as.factor(df.final.df.new.analysis$y) ~ my.mbl +my.mbw + my.mbd,
#          data = df.final.df.new.analysis, method = "qda", prec = 50, 
#          col.correct = "black",
#          col.wrong = "red", 
#          imageplot = TRUE,
#          # gs = as.numeric((df.final.df.new.analysis$y)),
#          pch.mean = 1:length(levels(as.factor(df.final.df.new.analysis$y))),
#          cex.mean = 1.4,
#          # image.colors = couleur,
#          col.mean = "yellow",
#          plot.matrix = TRUE)
if (pdf.output) {
  dev.off()
}


# Relation with environment -----------------------------------------------
my.eco.evo.df$preci.yr1

# Extract only linear 
if (!linear.only) {
  lin = dt[seq(1,nrow(dt),by = 2),]
} else {lin = dt}

# lin = lin[lin$lgt.p.1.t<= 0.05,]
plot(lin$bavggrd~my.eco.evo.df$preci.yr2)
plot(lin$bavggrd.mean.st~my.eco.evo.df$preci.yr2)

# Extract only quadratic 
quad = dt[seq(2,nrow(dt),by = 2),]

# Here I wanted to know what was the sample size for males and females all at once. 
yr = c(2004:2006,2007, 2008:2011)
for(i in 1:c(length(yr)-1)){
  yr.list <- c(yr[i],yr[i+jump])
  mdat <- prep.data(sp.keep = sp.list,
                    center = FALSE,yr.keep = yr.list,
                    site.keep = site.list,flip.pc1 = TRUE,
                    adults.only = FALSE,gr.sel = NULL,
                    valley.fortis = FALSE, valley.fortis.fuli = FALSE,
                    peaks.fortis.fuli = FALSE,two.univariate.traits = FALSE,
                    keep.resid = FALSE, splines = splines,
                    keep.last.year.data = FALSE,gam.analysis = TRUE,
                    K.nots = 5, # Manually select the number of knots  
                    traits = c('PC1'))#, 'PC2')) # This is only for the spline 
  tab.int = data.frame(surv = mdat$X[,2],
                       sex = droplevels(mdat$ind.vars$age),
                       pc1 = mdat$ind.vars$pc1)
  yearly.number.of.id = c(yearly.number.of.id,list(tab.int))
} 
new.yearly.list = NULL
for (i in 1:7) {
  vector.find.cut.each.year = yearly.number.of.id[[i]]$pc1>= all.ranges.x[[i]][1]  & 
    yearly.number.of.id[[i]]$pc1<= all.ranges.x[[i]][2]
  new.yearly.list = c(new.yearly.list,list(yearly.number.of.id[[i]] [vector.find.cut.each.year,]))
}

# Without cutting anything from the data 
lapply(yearly.number.of.id, function(x) table(x$sex))
lapply(yearly.number.of.id, function(x) table(x$surv))
lapply(new.yearly.list, function(x) table(x$sex))
lapply(new.yearly.list, function(x) table(x$surv))

lapply(yearly.number.of.id, nrow)
par(mfrow=c(2,2))
if (!linear.only) {
  lm.o1 = lm(quad$bavggrd~my.eco.evo.df$preci.yr1)
  summary(lm.o1)
  plot(quad$bavggrd~my.eco.evo.df$preci.yr1);abline(lm.o1)
  
  lm.o2 = lm(quad$bavggrd.mean.st~my.eco.evo.df$preci.yr1)
  summary(lm.o2)
  plot(quad$bavggrd.mean.st~my.eco.evo.df$preci.yr1);abline(lm.o2)
  
  lm.o3 = lm(quad$bavggrd~my.eco.evo.df$preci.yr2)
  summary(lm.o3)
  plot(quad$bavggrd~my.eco.evo.df$preci.yr2);abline(lm.o3)
  
  lm.o4 = lm(quad$bavggrd.mean.st~my.eco.evo.df$preci.yr2)
  summary(lm.o4)
  plot(quad$bavggrd.mean.st~my.eco.evo.df$preci.yr2);abline(lm.o4)
}
if(pdf.output){
  dev.off()
  dev.off()
  dev.off()
  dev.off()
}